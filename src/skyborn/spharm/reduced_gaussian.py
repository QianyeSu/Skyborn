"""Experimental reduced-Gaussian backend scaffolding.

This module defines the Python-side data model and method signatures for a
future ECMWF ectrans-backed spherical-harmonic backend that can operate on
packed reduced Gaussian grids directly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.polynomial.legendre import leggauss

from . import _spherepack

try:
    from .ectrans_backend_api import (
        create_setup,
        describe_setup,
        destroy_setup,
        gradient_synthesis_stub,
        scalar_analysis_stub,
        scalar_block_solve_stub,
        scalar_fourier_stub,
        scalar_synthesis_stub,
        uv_synthesis_stub,
        validate_nloen,
        vrtdiv_analysis_stub,
        weighted_block_solve_stub,
    )

    _ECTRANS_BACKEND_AVAILABLE = True
except ImportError:
    _ECTRANS_BACKEND_AVAILABLE = False

DEFAULT_EARTH_RADIUS = 6.3712e6


class ReducedGaussianBackendError(RuntimeError):
    """Base error for the experimental reduced-Gaussian backend."""


class ReducedGaussianValidationError(ValueError):
    """Raised when reduced-grid metadata or data shapes are invalid."""


@dataclass(frozen=True)
class ReducedGaussianGrid:
    """Describe one reduced Gaussian horizontal grid.

    Parameters
    ----------
    nloen:
        One-dimensional integer array giving the number of longitude points on
        each latitude circle, ordered north-to-south.
    rsphere:
        Sphere radius in meters.
    """

    nloen: np.ndarray
    rsphere: float = DEFAULT_EARTH_RADIUS

    def __post_init__(self) -> None:
        nloen = np.asarray(self.nloen, dtype=np.int32)
        if nloen.ndim != 1:
            raise ReducedGaussianValidationError(
                f"nloen must be rank 1, got rank {nloen.ndim}"
            )
        if nloen.size < 3:
            raise ReducedGaussianValidationError(
                f"nloen must describe at least 3 latitudes, got {nloen.size}"
            )
        if np.any(nloen <= 0):
            raise ReducedGaussianValidationError("nloen entries must be positive")
        object.__setattr__(self, "rsphere", float(self.rsphere))
        if _ECTRANS_BACKEND_AVAILABLE:
            native_ndgl, native_ngptot = validate_nloen(nloen)
            if native_ndgl != int(nloen.size):
                raise ReducedGaussianValidationError(
                    "native reduced-grid validator returned inconsistent ndgl"
                )
            if native_ngptot != int(np.sum(nloen, dtype=np.int64)):
                raise ReducedGaussianValidationError(
                    "native reduced-grid validator returned inconsistent ngptot"
                )
        object.__setattr__(self, "nloen", nloen)

    @property
    def ndgl(self) -> int:
        """Return the number of Gaussian latitude circles."""
        return int(self.nloen.size)

    @property
    def ngptot(self) -> int:
        """Return the total number of packed grid points."""
        return int(np.sum(self.nloen, dtype=np.int64))

    @property
    def max_nlon(self) -> int:
        """Return the maximum longitude count on any latitude circle."""
        return int(np.max(self.nloen))


class ReducedGaussianSpharmt:
    """Experimental reduced-Gaussian spherical-harmonic transform interface.

    This class mirrors the shape of the public `Spharmt` API where practical,
    but it operates on packed reduced-Gaussian arrays of shape
    `(ngptot,)` or `(ngptot, *extra_dims)` instead of rectangular
    `(nlat, nlon, *extra_dims)` arrays.

    The actual native ectrans transform implementation is not wired yet; this
    class currently provides validation, shape handling, and method contracts.
    """

    def __init__(
        self,
        nloen: np.ndarray,
        rsphere: float = DEFAULT_EARTH_RADIUS,
        backend: str = "ectrans",
    ) -> None:
        self.grid = ReducedGaussianGrid(nloen=nloen, rsphere=rsphere)
        if backend != "ectrans":
            raise ReducedGaussianValidationError(
                f"backend must be 'ectrans', got {backend!r}"
            )
        self.backend = backend
        self._setup_handle = None
        self._gaussian_mu = None
        self._gaussian_weights = None
        self._gaussian_lats = None
        self._lat_offsets = None
        self._bridge_spharmt = None
        self._analysis_basis_cache = {}
        self._vector_analysis_basis_cache = {}
        self._spectral_degree_index_cache = {}
        self._native_uv_synthesis_available = None
        self._native_vrtdiv_analysis_available = None
        self._native_gradient_synthesis_available = None

    @property
    def ndgl(self) -> int:
        """Return the number of Gaussian latitude circles."""
        return self.grid.ndgl

    @property
    def nloen(self) -> np.ndarray:
        """Return longitude counts per latitude circle."""
        return self.grid.nloen

    @property
    def ngptot(self) -> int:
        """Return total packed grid-point count."""
        return self.grid.ngptot

    @property
    def rsphere(self) -> float:
        """Return sphere radius in meters."""
        return self.grid.rsphere

    def __repr__(self) -> str:
        return (
            "ReducedGaussianSpharmt("
            f"ndgl={self.ndgl}, "
            f"ngptot={self.ngptot}, "
            f"max_nlon={self.grid.max_nlon}, "
            f"rsphere={self.rsphere:e})"
        )

    def _validate_ntrunc(self, ntrunc: Optional[int]) -> int:
        """Validate and normalize triangular truncation."""
        max_allowed = self.ndgl - 1
        if ntrunc is None:
            return max_allowed
        if ntrunc < 0 or ntrunc > max_allowed:
            raise ReducedGaussianValidationError(
                f"ntrunc must be between 0 and {max_allowed}, got {ntrunc}"
            )
        return ntrunc

    @staticmethod
    def _ncoeff_from_ntrunc(ntrunc: int) -> int:
        """Return the triangular complex coefficient count for one truncation."""
        return ((ntrunc + 1) * (ntrunc + 2)) // 2

    def _infer_ntrunc_from_ncoeff(self, ncoeff: int, operation_name: str) -> int:
        """Infer triangular truncation from public complex coefficient count."""
        ntrunc = int(-1.5 + 0.5 * math.sqrt(1.0 + 8.0 * ncoeff))
        if self._ncoeff_from_ntrunc(ntrunc) != ncoeff:
            raise ReducedGaussianValidationError(
                f"{operation_name} got inconsistent spectral leading size {ncoeff}"
            )
        if ntrunc > self.ndgl - 1:
            raise ReducedGaussianValidationError(
                f"{operation_name} inferred ntrunc {ntrunc}, but max allowed is {self.ndgl - 1}"
            )
        return ntrunc

    @classmethod
    def _ectrans_scalar_nspec2_from_ntrunc(cls, ntrunc: int) -> int:
        """Return ECMWF-style packed scalar spectral size (real-imag interleaved)."""
        return 2 * cls._ncoeff_from_ntrunc(ntrunc)

    def _validate_packed_grid_data(
        self, data: np.ndarray, operation_name: str
    ) -> Tuple[int, np.ndarray, Tuple[int, ...]]:
        """Validate packed reduced-grid data and flatten extra dimensions."""
        array = np.asarray(data)
        if array.ndim < 1:
            raise ReducedGaussianValidationError(
                f"{operation_name} needs at least a rank 1 array, got rank {array.ndim}"
            )
        if array.shape[0] != self.ngptot:
            raise ReducedGaussianValidationError(
                f"{operation_name} needs leading dimension {self.ngptot}, got {array.shape[0]}"
            )
        if array.ndim == 1:
            nt = 1
            normalized = np.expand_dims(array, 1)
            extra_shape: Tuple[int, ...] = ()
        else:
            extra_shape = tuple(array.shape[1:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized = array.reshape(self.ngptot, nt)
        return nt, normalized, extra_shape

    def _validate_spectral_data(
        self, data: np.ndarray, operation_name: str
    ) -> Tuple[int, int, np.ndarray, Tuple[int, ...]]:
        """Validate spectral data and flatten trailing dimensions."""
        array = np.asarray(data)
        if array.ndim < 1:
            raise ReducedGaussianValidationError(
                f"{operation_name} needs at least a rank 1 array, got rank {array.ndim}"
            )
        ncoeff = array.shape[0]
        ntrunc = self._infer_ntrunc_from_ncoeff(ncoeff, operation_name)
        if array.ndim == 1:
            nt = 1
            normalized = np.expand_dims(array, 1)
            extra_shape: Tuple[int, ...] = ()
        else:
            extra_shape = tuple(array.shape[1:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized = array.reshape(ncoeff, nt)
        return nt, ntrunc, normalized, extra_shape

    def _validate_ectrans_scalar_packed(
        self, data: np.ndarray, operation_name: str
    ) -> Tuple[int, int, np.ndarray, Tuple[int, ...]]:
        """Validate ECMWF-style packed scalar spectra and flatten trailing dimensions."""
        array = np.asarray(data)
        if array.ndim < 1:
            raise ReducedGaussianValidationError(
                f"{operation_name} needs at least a rank 1 array, got rank {array.ndim}"
            )
        nspec2 = array.shape[0]
        if nspec2 % 2 != 0:
            raise ReducedGaussianValidationError(
                f"{operation_name} needs even packed spectral leading size, got {nspec2}"
            )
        ncoeff = nspec2 // 2
        ntrunc = self._infer_ntrunc_from_ncoeff(ncoeff, operation_name)
        if array.ndim == 1:
            nt = 1
            normalized = np.expand_dims(array, 1)
            extra_shape: Tuple[int, ...] = ()
        else:
            extra_shape = tuple(array.shape[1:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized = array.reshape(nspec2, nt)
        return nt, ntrunc, normalized, extra_shape

    @staticmethod
    def _restore_grid_shape(
        datagrid: np.ndarray, extra_shape: Tuple[int, ...]
    ) -> np.ndarray:
        """Restore packed grid fields to their original extra dimensions."""
        if not extra_shape:
            return datagrid[:, 0]
        return datagrid.reshape((datagrid.shape[0],) + extra_shape)

    @staticmethod
    def _restore_spectral_shape(
        dataspec: np.ndarray, extra_shape: Tuple[int, ...]
    ) -> np.ndarray:
        """Restore spectral fields to their original extra dimensions."""
        if not extra_shape:
            return dataspec[:, 0]
        return dataspec.reshape((dataspec.shape[0],) + extra_shape)

    def _spectogrd_pair(
        self,
        spec_a: np.ndarray,
        spec_b: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Synthesize two scalar spectra together through one batched call."""
        if np.asarray(spec_a).shape != np.asarray(spec_b).shape:
            raise ReducedGaussianValidationError(
                "paired spectra must have the same shape"
            )
        paired_specs = np.stack((spec_a, spec_b), axis=-1)
        paired_grids = self.spectogrd(paired_specs)
        return paired_grids[..., 0], paired_grids[..., 1]

    def _pack_scalar_spectrum_for_ectrans(
        self, dataspec: np.ndarray, operation_name: str
    ) -> Tuple[np.ndarray, int, Tuple[int, ...]]:
        """Convert public complex triangular spectra to ECMWF packed real layout."""
        _, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, operation_name
        )
        packed = np.empty(
            (self._ectrans_scalar_nspec2_from_ntrunc(ntrunc), normalized_spec.shape[1]),
            dtype=np.float64,
        )
        packed[0::2, :] = normalized_spec.real
        packed[1::2, :] = normalized_spec.imag
        return self._restore_spectral_shape(packed, extra_shape), ntrunc, extra_shape

    def _unpack_scalar_spectrum_from_ectrans(
        self,
        packedspec: np.ndarray,
        operation_name: str,
    ) -> np.ndarray:
        """Convert ECMWF packed real spectra back to public complex triangular layout."""
        _, _, normalized_packed, extra_shape = self._validate_ectrans_scalar_packed(
            packedspec, operation_name
        )
        complex_spec = normalized_packed[0::2, :] + 1j * normalized_packed[1::2, :]
        return self._restore_spectral_shape(complex_spec, extra_shape)

    def _not_implemented(self, method_name: str) -> None:
        """Raise a consistent error for not-yet-wired native transforms."""
        raise ReducedGaussianBackendError(
            f"{method_name} is not implemented yet for the experimental "
            "ectrans reduced-Gaussian backend scaffold."
        )

    def _require_setup_handle(self):
        """Create and cache a native setup handle on first use."""
        if not _ECTRANS_BACKEND_AVAILABLE:
            self._not_implemented("setup")
        if self._setup_handle is None:
            mu, weights, _ = self._get_gaussian_metadata()
            self._setup_handle = create_setup(self.nloen, mu, weights, self.rsphere)
        return self._setup_handle

    def _compute_gaussian_mu_weights(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute Gaussian sin(latitude) roots and quadrature weights."""
        mu, weights = leggauss(self.ndgl)
        mu = np.asarray(mu[::-1], dtype=np.float64)
        weights = np.asarray(weights[::-1], dtype=np.float64)
        return mu, weights

    def _get_gaussian_metadata(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return cached Gaussian mu, weights, and latitudes."""
        if self._gaussian_mu is None or self._gaussian_weights is None:
            mu, weights = self._compute_gaussian_mu_weights()
            self._gaussian_mu = mu
            self._gaussian_weights = weights
            self._gaussian_lats = np.degrees(np.arcsin(mu))
        return self._gaussian_mu, self._gaussian_weights, self._gaussian_lats

    def _get_lat_offsets(self) -> np.ndarray:
        """Return cached packed-grid latitude starting offsets."""
        if self._lat_offsets is None:
            lat_offsets = np.empty(self.ndgl, dtype=np.int64)
            lat_offsets[0] = 0
            if self.ndgl > 1:
                lat_offsets[1:] = np.cumsum(self.nloen[:-1], dtype=np.int64)
            self._lat_offsets = lat_offsets
        return self._lat_offsets

    @staticmethod
    def _resample_longitude_periodic(data: np.ndarray, target_nlon: int) -> np.ndarray:
        """Resample periodic longitude data using truncated/zero-padded real Fourier modes."""
        array = np.asarray(data, dtype=np.float64)
        squeeze = False
        if array.ndim == 1:
            array = array[:, None]
            squeeze = True
        source_nlon = array.shape[0]
        if source_nlon == target_nlon:
            result = array.copy()
            return result[:, 0] if squeeze else result

        coeff = np.fft.rfft(array, axis=0)
        target_nfreq = target_nlon // 2 + 1
        resized = np.zeros((target_nfreq, array.shape[1]), dtype=np.complex128)
        copy_nfreq = min(coeff.shape[0], target_nfreq)
        resized[:copy_nfreq, :] = coeff[:copy_nfreq, :]
        resized *= float(target_nlon) / float(source_nlon)
        result = np.fft.irfft(resized, n=target_nlon, axis=0)
        return result[:, 0] if squeeze else result

    def _get_bridge_spharmt(self):
        """Return a cached full-Gaussian bridge transform on the maximum longitude count."""
        if self._bridge_spharmt is None:
            from .spherical_harmonics import Spharmt

            self._bridge_spharmt = Spharmt(
                int(self.grid.max_nlon),
                self.ndgl,
                rsphere=self.rsphere,
                gridtype="gaussian",
                legfunc="stored",
            )
        return self._bridge_spharmt

    def _packed_to_full_gaussian_grid(self, normalized_grid: np.ndarray) -> np.ndarray:
        """Resample packed reduced-grid data onto a full Gaussian rectangular grid."""
        nt = normalized_grid.shape[1]
        target_nlon = int(self.grid.max_nlon)
        full_grid = np.empty((self.ndgl, target_nlon, nt), dtype=np.float64)
        lat_offsets = self._get_lat_offsets()
        for ilat, nlon in enumerate(self.nloen):
            offset = int(lat_offsets[ilat])
            full_grid[ilat, :, :] = self._resample_longitude_periodic(
                normalized_grid[offset : offset + int(nlon), :],
                target_nlon,
            )
        return full_grid

    def _spectral_column_slices(self, ntrunc: int) -> list[slice]:
        """Return triangular coefficient slices grouped by zonal wavenumber."""
        slices: list[slice] = []
        start = 0
        for m in range(ntrunc + 1):
            count = ntrunc - m + 1
            slices.append(slice(start, start + count))
            start += count
        return slices

    def _get_spectral_degree_index(self, ntrunc: int) -> np.ndarray:
        """Return cached total spherical wavenumber index per packed coefficient."""
        cached = self._spectral_degree_index_cache.get(ntrunc)
        if cached is not None:
            return cached
        degree_index = np.empty((self._ncoeff_from_ntrunc(ntrunc),), dtype=np.int32)
        start = 0
        for m in range(ntrunc + 1):
            count = ntrunc - m + 1
            degree_index[start : start + count] = np.arange(
                m,
                ntrunc + 1,
                dtype=np.int32,
            )
            start += count
        self._spectral_degree_index_cache[ntrunc] = degree_index
        return degree_index

    def _scalar_fourier(
        self,
        normalized_grid: np.ndarray,
        mmax: int,
    ) -> np.ndarray:
        """Return native reduced-grid Fourier coefficients shaped as (ndgl, mmax + 1, nt)."""
        setup_handle = self._require_setup_handle()
        fourier, ierror = scalar_fourier_stub(setup_handle, normalized_grid, mmax)
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"native scalar_fourier_stub failed with ierror={ierror}"
            )
        return np.asarray(fourier, dtype=np.complex128)

    @staticmethod
    def _split_columnwise(
        array: np.ndarray,
        counts: tuple[int, ...],
    ) -> tuple[np.ndarray, ...]:
        """Split an array by consecutive blocks along its last axis."""
        start = 0
        chunks = []
        for count in counts:
            stop = start + count
            chunks.append(array[..., start:stop])
            start = stop
        return tuple(chunks)

    def _scalar_block_solve(
        self,
        basis_block: np.ndarray,
        observed: np.ndarray,
    ) -> np.ndarray:
        """Return one native weighted least-squares block solve result."""
        setup_handle = self._require_setup_handle()
        solution, ierror = scalar_block_solve_stub(setup_handle, basis_block, observed)
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"native scalar_block_solve_stub failed with ierror={ierror}"
            )
        return np.asarray(solution, dtype=np.complex128)

    @staticmethod
    def _weighted_block_solve(
        weights: np.ndarray,
        basis_block: np.ndarray,
        observed: np.ndarray,
    ) -> np.ndarray:
        """Return one native general weighted least-squares solve result."""
        solution, ierror = weighted_block_solve_stub(
            np.asarray(weights, dtype=np.float64),
            basis_block,
            observed,
        )
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"native weighted_block_solve_stub failed with ierror={ierror}"
            )
        return np.asarray(solution, dtype=np.complex128)

    def _get_analysis_basis(self, ntrunc: int) -> tuple[list[slice], list[np.ndarray]]:
        """Return cached reduced-grid Fourier/latitude basis blocks from native synthesis."""
        cached = self._analysis_basis_cache.get(ntrunc)
        if cached is not None:
            return cached

        slices = self._spectral_column_slices(ntrunc)
        basis_blocks: list[np.ndarray] = []
        ncoeff = self._ncoeff_from_ntrunc(ntrunc)
        basis_spec = np.eye(ncoeff, dtype=np.complex128)
        basis_grid = np.asarray(self.spectogrd(basis_spec), dtype=np.float64)
        basis_fourier = self._scalar_fourier(basis_grid, ntrunc)

        for m, spectral_slice in enumerate(slices):
            basis_blocks.append(
                np.asarray(
                    basis_fourier[:, m, spectral_slice],
                    dtype=np.complex128,
                )
            )

        cached = (slices, basis_blocks)
        self._analysis_basis_cache[ntrunc] = cached
        return cached

    @staticmethod
    def _solve_weighted_complex_lstsq(
        basis: np.ndarray,
        observed: np.ndarray,
        weights: np.ndarray,
    ) -> np.ndarray:
        """Solve one weighted complex least-squares system."""
        sqrt_weights = np.sqrt(np.asarray(weights, dtype=np.float64))[:, None]
        weighted_basis = sqrt_weights * np.asarray(basis, dtype=np.complex128)
        weighted_observed = sqrt_weights * np.asarray(observed, dtype=np.complex128)
        solution, _, _, _ = np.linalg.lstsq(
            weighted_basis,
            weighted_observed,
            rcond=None,
        )
        return np.asarray(solution, dtype=np.complex128)

    def _get_vector_analysis_basis(
        self,
        ntrunc: int,
    ) -> tuple[list[slice], list[np.ndarray]]:
        """Return cached reduced-grid vector-analysis basis blocks."""
        use_native_uv = self._supports_native_uv_synthesis()
        cache_key = (ntrunc, use_native_uv)
        cached = self._vector_analysis_basis_cache.get(cache_key)
        if cached is not None:
            return cached

        spectral_slices = self._spectral_column_slices(ntrunc)
        basis_blocks: list[np.ndarray] = []
        ncoeff = self._ncoeff_from_ntrunc(ntrunc)
        basis_spec = np.eye(ncoeff, dtype=np.complex128)
        zero_spec = np.zeros_like(basis_spec)
        joint_vrtspec = np.concatenate((basis_spec, zero_spec), axis=1)
        joint_divspec = np.concatenate((zero_spec, basis_spec), axis=1)
        joint_u, joint_v = self.getuv(joint_vrtspec, joint_divspec)
        if use_native_uv and not self._supports_native_uv_synthesis():
            return self._get_vector_analysis_basis(ntrunc)
        joint_fourier = self._scalar_fourier(
            np.concatenate(
                (
                    np.asarray(joint_u, dtype=np.float64),
                    np.asarray(joint_v, dtype=np.float64),
                ),
                axis=1,
            ),
            ntrunc,
        )
        u_fourier_all, v_fourier_all = self._split_columnwise(
            np.asarray(joint_fourier, dtype=np.complex128),
            (2 * ncoeff, 2 * ncoeff),
        )

        for m, spectral_slice in enumerate(spectral_slices):
            vort_slice = slice(spectral_slice.start, spectral_slice.stop)
            div_slice = slice(
                ncoeff + spectral_slice.start,
                ncoeff + spectral_slice.stop,
            )
            u_basis = np.concatenate(
                (
                    u_fourier_all[:, m, vort_slice],
                    u_fourier_all[:, m, div_slice],
                ),
                axis=1,
            )
            v_basis = np.concatenate(
                (
                    v_fourier_all[:, m, vort_slice],
                    v_fourier_all[:, m, div_slice],
                ),
                axis=1,
            )
            basis_blocks.append(np.concatenate((u_basis, v_basis), axis=0))

        cached = (spectral_slices, basis_blocks)
        self._vector_analysis_basis_cache[cache_key] = cached
        return cached

    def _lapspec(
        self,
        dataspec: np.ndarray,
        operation_name: str = "_lapspec",
    ) -> np.ndarray:
        """Apply the Laplacian to a scalar spectrum and preserve shape."""
        _, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec,
            operation_name,
        )
        lap_spec = np.zeros_like(normalized_spec, dtype=np.complex64)
        rsphere_inv_sq = np.float32(1.0) / (
            np.float32(self.rsphere) * np.float32(self.rsphere)
        )
        start = 0
        for m in range(ntrunc + 1):
            n_values = np.arange(m, ntrunc + 1, dtype=np.float32)
            count = ntrunc - m + 1
            factors = -(n_values * (n_values + np.float32(1.0))) * rsphere_inv_sq
            lap_spec[start : start + count, :] = (
                normalized_spec[start : start + count, :].astype(np.complex64)
                * factors[:, None]
            )
            start += count
        return self._restore_spectral_shape(lap_spec, extra_shape)

    def _invlapspec(
        self,
        dataspec: np.ndarray,
        operation_name: str = "_invlapspec",
    ) -> np.ndarray:
        """Apply the inverse Laplacian to a scalar spectrum and preserve shape."""
        _, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec,
            operation_name,
        )
        invlap_spec = np.zeros_like(normalized_spec, dtype=np.complex64)
        rsphere_sq = np.float32(self.rsphere) * np.float32(self.rsphere)
        start = 0
        for m in range(ntrunc + 1):
            n_values = np.arange(m, ntrunc + 1, dtype=np.float32)
            count = ntrunc - m + 1
            factors = np.zeros((count,), dtype=np.float32)
            nonzero = n_values > np.float32(0.0)
            factors[nonzero] = -rsphere_sq / (
                n_values[nonzero] * (n_values[nonzero] + np.float32(1.0))
            )
            invlap_spec[start : start + count, :] = (
                normalized_spec[start : start + count, :].astype(np.complex64)
                * factors[:, None]
            )
            start += count
        return self._restore_spectral_shape(invlap_spec, extra_shape)

    def _apply_spectral_smoothing(
        self,
        dataspec: np.ndarray,
        smooth: np.ndarray,
        operation_name: str = "_apply_spectral_smoothing",
        expected_size: Optional[int] = None,
    ) -> np.ndarray:
        """Apply isotropic smoothing factors by total spherical wavenumber."""
        _, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec,
            operation_name,
        )
        smooth_array = np.asarray(smooth, dtype=np.float32)
        if expected_size is None:
            expected_size = self.ndgl
        if smooth_array.ndim != 1 or smooth_array.shape[0] != expected_size:
            raise ReducedGaussianValidationError(
                f"smooth must be rank 1 with size {expected_size}, got shape {smooth_array.shape}"
            )
        factors = smooth_array[self._get_spectral_degree_index(ntrunc)]
        smoothed_spec = normalized_spec.astype(np.complex64) * factors[:, None]
        return self._restore_spectral_shape(smoothed_spec, extra_shape)

    def _full_gaussian_to_packed_grid(self, full_grid: np.ndarray) -> np.ndarray:
        """Resample a full Gaussian rectangular grid down to the packed reduced layout."""
        nt = full_grid.shape[2]
        packed_grid = np.empty((self.ngptot, nt), dtype=np.float64)
        lat_offsets = self._get_lat_offsets()
        for ilat, nlon in enumerate(self.nloen):
            offset = int(lat_offsets[ilat])
            packed_grid[offset : offset + int(nlon), :] = (
                self._resample_longitude_periodic(
                    full_grid[ilat, :, :],
                    int(nlon),
                )
            )
        return packed_grid

    def _packed_to_bridge_grid(
        self,
        normalized_grid: np.ndarray,
        extra_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Convert packed reduced-grid data to the bridge full-Gaussian grid shape."""
        bridge = self._get_bridge_spharmt()
        full_grid = self._packed_to_full_gaussian_grid(normalized_grid)
        return bridge._restore_grid_shape(full_grid, extra_shape)

    def _packed_pair_to_bridge_grids(
        self,
        normalized_a: np.ndarray,
        normalized_b: np.ndarray,
        extra_shape: Tuple[int, ...],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Convert two packed reduced-grid fields to bridge-grid shape in one resampling pass."""
        if normalized_a.shape != normalized_b.shape:
            raise ReducedGaussianValidationError(
                "paired packed fields must have the same normalized shape"
            )
        bridge = self._get_bridge_spharmt()
        full_pair = self._packed_to_full_gaussian_grid(
            np.concatenate((normalized_a, normalized_b), axis=1)
        )
        restored_pair = bridge._restore_grid_shape(full_pair, extra_shape + (2,))
        return restored_pair[..., 0], restored_pair[..., 1]

    def _bridge_to_packed_grid(
        self,
        full_grid: np.ndarray,
        operation_name: str,
    ) -> np.ndarray:
        """Convert a bridge full-Gaussian grid field back to packed reduced layout."""
        bridge = self._get_bridge_spharmt()
        _, normalized_full_grid, extra_shape = bridge._validate_grid_data(
            full_grid,
            operation_name,
        )
        packed_grid = self._full_gaussian_to_packed_grid(normalized_full_grid)
        return self._restore_grid_shape(packed_grid, extra_shape)

    def _bridge_pair_to_packed_grids(
        self,
        full_a: np.ndarray,
        full_b: np.ndarray,
        operation_name: str,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Convert two bridge-grid fields back to packed reduced layout in one resampling pass."""
        if np.asarray(full_a).shape != np.asarray(full_b).shape:
            raise ReducedGaussianValidationError(
                "paired bridge-grid fields must have the same shape"
            )
        bridge = self._get_bridge_spharmt()
        _, normalized_a, extra_shape = bridge._validate_grid_data(
            full_a, operation_name
        )
        _, normalized_b, extra_shape_b = bridge._validate_grid_data(
            full_b,
            operation_name,
        )
        if extra_shape != extra_shape_b:
            raise ReducedGaussianValidationError(
                "paired bridge-grid fields must have consistent dimensions"
            )
        packed_pair = self._full_gaussian_to_packed_grid(
            np.concatenate((normalized_a, normalized_b), axis=2)
        )
        restored_pair = self._restore_grid_shape(packed_pair, extra_shape + (2,))
        return restored_pair[..., 0], restored_pair[..., 1]

    def _getvrtdivspec_bridge(
        self,
        normalized_u: np.ndarray,
        normalized_v: np.ndarray,
        extra_shape: Tuple[int, ...],
        ntrunc: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Analyze winds through the full-Gaussian bridge fallback path."""
        bridge_ugrid, bridge_vgrid = self._packed_pair_to_bridge_grids(
            normalized_u,
            normalized_v,
            extra_shape,
        )
        bridge = self._get_bridge_spharmt()
        return bridge.getvrtdivspec(
            bridge_ugrid,
            bridge_vgrid,
            ntrunc=ntrunc,
        )

    def _getuv_bridge(
        self,
        vrtspec: np.ndarray,
        divspec: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Synthesize winds through the full-Gaussian bridge fallback path."""
        bridge = self._get_bridge_spharmt()
        full_ugrid, full_vgrid = bridge.getuv(vrtspec, divspec)
        return self._bridge_pair_to_packed_grids(
            full_ugrid,
            full_vgrid,
            "reduced getuv bridge",
        )

    def _getgrad_bridge(
        self,
        chispec: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Synthesize gradient winds through the full-Gaussian bridge fallback path."""
        bridge = self._get_bridge_spharmt()
        full_ugrad, full_vgrad = bridge.getgrad(chispec)
        return self._bridge_pair_to_packed_grids(
            full_ugrad,
            full_vgrad,
            "reduced getgrad bridge",
        )

    def _try_native_vrtdiv_analysis(
        self,
        ugrid: np.ndarray,
        vgrid: np.ndarray,
        ntrunc: int,
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Try native reduced-grid vrt/div analysis and return None on stub fallback."""
        if not self._supports_native_vrtdiv_analysis():
            return None
        setup_handle = self._require_setup_handle()
        vrtspec, divspec, ierror = vrtdiv_analysis_stub(
            setup_handle,
            np.asarray(ugrid, dtype=np.float64),
            np.asarray(vgrid, dtype=np.float64),
            ntrunc,
        )
        if ierror == -1:
            self._native_vrtdiv_analysis_available = False
            return None
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"getvrtdivspec native stub returned ierror={ierror}"
            )
        self._native_vrtdiv_analysis_available = True
        return vrtspec, divspec

    def _try_native_uv_synthesis(
        self,
        vrtspec: np.ndarray,
        divspec: np.ndarray,
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Try native reduced-grid uv synthesis and return None on stub fallback."""
        if not self._supports_native_uv_synthesis():
            return None
        setup_handle = self._require_setup_handle()
        ugrid, vgrid, ierror = uv_synthesis_stub(
            setup_handle,
            np.asarray(vrtspec, dtype=np.complex128),
            np.asarray(divspec, dtype=np.complex128),
        )
        if ierror == -1:
            self._native_uv_synthesis_available = False
            return None
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"getuv native stub returned ierror={ierror}"
            )
        self._native_uv_synthesis_available = True
        return ugrid, vgrid

    def _try_native_gradient_synthesis(
        self,
        chispec: np.ndarray,
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Try native reduced-grid gradient synthesis and return None on stub fallback."""
        if not self._supports_native_gradient_synthesis():
            return None
        setup_handle = self._require_setup_handle()
        ugrad, vgrad, ierror = gradient_synthesis_stub(
            setup_handle,
            np.asarray(chispec, dtype=np.complex128),
        )
        if ierror == -1:
            self._native_gradient_synthesis_available = False
            return None
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"getgrad native stub returned ierror={ierror}"
            )
        self._native_gradient_synthesis_available = True
        return ugrad, vgrad

    def _try_native_scalar_analysis(
        self,
        datagrid: np.ndarray,
        ntrunc: int,
    ) -> Optional[np.ndarray]:
        """Try native reduced-grid scalar analysis and return None on stub fallback."""
        setup_handle = self._require_setup_handle()
        dataspec_packed, ierror = scalar_analysis_stub(
            setup_handle,
            np.asarray(datagrid, dtype=np.float64),
            ntrunc,
        )
        if ierror == -1:
            return None
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"grdtospec native stub returned ierror={ierror}"
            )
        return dataspec_packed

    def _try_native_scalar_synthesis(
        self,
        packedspec: np.ndarray,
    ) -> Optional[np.ndarray]:
        """Try native reduced-grid scalar synthesis and return None on stub fallback."""
        setup_handle = self._require_setup_handle()
        datagrid, ierror = scalar_synthesis_stub(
            setup_handle,
            np.asarray(packedspec, dtype=np.float64),
        )
        if ierror == -1:
            return None
        if ierror != 0:
            raise ReducedGaussianBackendError(
                f"spectogrd native stub returned ierror={ierror}"
            )
        return datagrid

    def _spectogrd_specintrp_prototype(
        self,
        normalized_spec: np.ndarray,
        ntrunc: int,
        extra_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Prototype reduced-grid synthesis using a full Gaussian bridge grid."""
        bridge = self._get_bridge_spharmt()
        bridge_spec = bridge._restore_spectral_shape(normalized_spec, extra_shape)
        full_grid = bridge.spectogrd(bridge_spec)
        _, normalized_full_grid, _ = bridge._validate_grid_data(
            full_grid, "reduced spectogrd bridge"
        )
        packed_grid = self._full_gaussian_to_packed_grid(normalized_full_grid)
        return self._restore_grid_shape(packed_grid, extra_shape)

    def _grdtospec_prototype(
        self,
        normalized_grid: np.ndarray,
        ntrunc: int,
        extra_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Prototype reduced-grid analysis using native-synthesis-consistent basis blocks."""
        nt = normalized_grid.shape[1]
        ncoeff = self._ncoeff_from_ntrunc(ntrunc)
        dataspec = np.zeros((ncoeff, nt), dtype=np.complex128)
        spectral_slices, basis_blocks = self._get_analysis_basis(ntrunc)
        observed_fourier = self._scalar_fourier(normalized_grid, ntrunc)

        for m, (spectral_slice, basis_block) in enumerate(
            zip(spectral_slices, basis_blocks)
        ):
            observed = observed_fourier[:, m, :]
            try:
                solution = self._scalar_block_solve(basis_block, observed)
            except ReducedGaussianBackendError:
                _, weights, _ = self._get_gaussian_metadata()
                sqrt_weights = np.sqrt(weights)[:, None]
                weighted_basis = sqrt_weights * basis_block
                weighted_observed = sqrt_weights * observed
                solution, _, _, _ = np.linalg.lstsq(
                    weighted_basis,
                    weighted_observed,
                    rcond=None,
                )
            dataspec[spectral_slice, :] = solution

        return self._restore_spectral_shape(dataspec, extra_shape)

    def _getvrtdivspec_prototype(
        self,
        normalized_u: np.ndarray,
        normalized_v: np.ndarray,
        ntrunc: int,
        extra_shape: Tuple[int, ...],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Prototype reduced-grid vector analysis using synthesis-consistent basis blocks."""
        nt = normalized_u.shape[1]
        ncoeff = self._ncoeff_from_ntrunc(ntrunc)
        vrtspec = np.zeros((ncoeff, nt), dtype=np.complex128)
        divspec = np.zeros((ncoeff, nt), dtype=np.complex128)
        spectral_slices, basis_blocks = self._get_vector_analysis_basis(ntrunc)
        combined_fourier = self._scalar_fourier(
            np.concatenate((normalized_u, normalized_v), axis=1),
            ntrunc,
        )
        u_fourier, v_fourier = self._split_columnwise(
            np.asarray(combined_fourier, dtype=np.complex128),
            (nt, nt),
        )
        _, weights, _ = self._get_gaussian_metadata()
        vector_weights = np.concatenate((weights, weights))

        for m, (spectral_slice, basis_block) in enumerate(
            zip(spectral_slices, basis_blocks)
        ):
            observed = np.concatenate(
                (
                    np.asarray(u_fourier[:, m, :], dtype=np.complex128),
                    np.asarray(v_fourier[:, m, :], dtype=np.complex128),
                ),
                axis=0,
            )
            try:
                solution = self._weighted_block_solve(
                    vector_weights,
                    basis_block,
                    observed,
                )
            except ReducedGaussianBackendError:
                solution = self._solve_weighted_complex_lstsq(
                    basis_block,
                    observed,
                    vector_weights,
                )
            nblock = spectral_slice.stop - spectral_slice.start
            vrtspec[spectral_slice, :] = solution[:nblock, :]
            divspec[spectral_slice, :] = solution[nblock:, :]

        return (
            self._restore_spectral_shape(vrtspec, extra_shape),
            self._restore_spectral_shape(divspec, extra_shape),
        )

    def close(self) -> None:
        """Release the native reduced-grid setup handle if it exists."""
        if self._setup_handle is not None:
            destroy_setup(self._setup_handle)
            self._setup_handle = None

    def _describe_setup(self) -> dict:
        """Return internal native setup metadata for focused testing/debugging."""
        setup_handle = self._require_setup_handle()
        return describe_setup(setup_handle)

    def _supports_native_uv_synthesis(self) -> bool:
        """Return whether native reduced-grid uv synthesis appears to be implemented."""
        if self._native_uv_synthesis_available is None:
            self._native_uv_synthesis_available = _ECTRANS_BACKEND_AVAILABLE
        return self._native_uv_synthesis_available

    def _supports_native_vrtdiv_analysis(self) -> bool:
        """Return whether a bridge-free reduced-grid vrt/div analysis path is enabled."""
        if self._native_vrtdiv_analysis_available is None:
            self._native_vrtdiv_analysis_available = _ECTRANS_BACKEND_AVAILABLE
        return self._native_vrtdiv_analysis_available

    def _supports_native_gradient_synthesis(self) -> bool:
        """Return whether native reduced-grid gradient synthesis appears to be implemented."""
        if self._native_gradient_synthesis_available is None:
            self._native_gradient_synthesis_available = _ECTRANS_BACKEND_AVAILABLE
        return self._native_gradient_synthesis_available

    def __del__(self) -> None:
        """Best-effort native setup cleanup."""
        try:
            self.close()
        except Exception:
            pass

    def grdtospec(
        self, datagrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Analyze packed reduced-grid scalar data to spherical-harmonic space."""
        _, normalized_grid, extra_shape = self._validate_packed_grid_data(
            datagrid, "grdtospec"
        )
        ntrunc_value = self._validate_ntrunc(ntrunc)
        native_result = self._try_native_scalar_analysis(datagrid, ntrunc_value)
        if native_result is not None:
            return self._unpack_scalar_spectrum_from_ectrans(
                native_result,
                "grdtospec",
            )
        return self._grdtospec_prototype(
            normalized_grid,
            ntrunc_value,
            extra_shape,
        )

    def spectogrd(self, dataspec: np.ndarray) -> np.ndarray:
        """Synthesize scalar spectra back to a packed reduced Gaussian grid."""
        nt, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, "spectogrd"
        )
        packed_spec, _, _ = self._pack_scalar_spectrum_for_ectrans(
            dataspec, "spectogrd"
        )
        native_result = self._try_native_scalar_synthesis(packed_spec)
        if native_result is not None:
            return native_result
        return self._spectogrd_specintrp_prototype(
            normalized_spec,
            ntrunc,
            extra_shape,
        )

    def getvrtdivspec(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Analyze packed reduced-grid winds to vorticity/divergence spectra."""
        ugrid_array = np.asarray(ugrid)
        vgrid_array = np.asarray(vgrid)
        if ugrid_array.shape != vgrid_array.shape:
            raise ReducedGaussianValidationError(
                "getvrtdivspec needs ugrid and vgrid with the same shape"
            )
        _, normalized_u, extra_shape = self._validate_packed_grid_data(
            ugrid_array, "getvrtdivspec"
        )
        _, normalized_v, v_extra_shape = self._validate_packed_grid_data(
            vgrid_array, "getvrtdivspec"
        )
        if extra_shape != v_extra_shape:
            raise ReducedGaussianValidationError(
                "getvrtdivspec needs ugrid and vgrid with consistent dimensions"
            )
        ntrunc_value = self._validate_ntrunc(ntrunc)
        native_result = self._try_native_vrtdiv_analysis(
            ugrid_array,
            vgrid_array,
            ntrunc_value,
        )
        if native_result is not None:
            return native_result
        return self._getvrtdivspec_bridge(
            normalized_u,
            normalized_v,
            extra_shape,
            ntrunc_value,
        )

    def getvrtspec(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Analyze packed reduced-grid winds to vorticity spectrum only."""
        vrtspec, _ = self.getvrtdivspec(ugrid, vgrid, ntrunc=ntrunc)
        return vrtspec

    def getdivspec(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Analyze packed reduced-grid winds to divergence spectrum only."""
        _, divspec = self.getvrtdivspec(ugrid, vgrid, ntrunc=ntrunc)
        return divspec

    def getuv(
        self, vrtspec: np.ndarray, divspec: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Synthesize packed reduced-grid winds from vorticity/divergence spectra."""
        vrtspec_array = np.asarray(vrtspec)
        divspec_array = np.asarray(divspec)
        if vrtspec_array.shape != divspec_array.shape:
            raise ReducedGaussianValidationError(
                "getuv needs vrtspec and divspec with the same shape"
            )
        nt_vrt, ntrunc_vrt, _, extra_shape = self._validate_spectral_data(
            vrtspec_array, "getuv"
        )
        nt_div, ntrunc_div, _, div_extra_shape = self._validate_spectral_data(
            divspec_array, "getuv"
        )
        if (
            nt_vrt != nt_div
            or ntrunc_vrt != ntrunc_div
            or extra_shape != div_extra_shape
        ):
            raise ReducedGaussianValidationError(
                "getuv needs vrtspec and divspec with consistent dimensions"
            )
        native_result = self._try_native_uv_synthesis(
            vrtspec_array,
            divspec_array,
        )
        if native_result is not None:
            return native_result
        return self._getuv_bridge(vrtspec_array, divspec_array)

    def getgrad(self, chispec: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Synthesize the packed reduced-grid horizontal gradient of a scalar spectrum."""
        chispec_array = np.asarray(chispec)
        self._validate_spectral_data(chispec_array, "getgrad")
        native_result = self._try_native_gradient_synthesis(chispec_array)
        if native_result is not None:
            return native_result
        return self._getgrad_bridge(chispec_array)

    def getpsichi_spec_component(
        self,
        ugrid: np.ndarray,
        vgrid: np.ndarray,
        ntrunc: Optional[int],
        component: str,
        operation_name: str,
    ) -> np.ndarray:
        """Compute one streamfunction/velocity-potential spectrum."""
        source_spec = (
            self.getvrtspec(ugrid, vgrid, ntrunc=ntrunc)
            if component == "psi"
            else self.getdivspec(ugrid, vgrid, ntrunc=ntrunc)
        )
        return self._invlapspec(source_spec, operation_name)

    def getpsispec(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Compute streamfunction spectral coefficients from packed reduced-grid wind."""
        return self.getpsichi_spec_component(
            ugrid,
            vgrid,
            ntrunc=ntrunc,
            component="psi",
            operation_name="getpsispec",
        )

    def getchispec(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Compute velocity-potential spectral coefficients from packed reduced-grid wind."""
        return self.getpsichi_spec_component(
            ugrid,
            vgrid,
            ntrunc=ntrunc,
            component="chi",
            operation_name="getchispec",
        )

    def getpsichispec(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute streamfunction and velocity-potential spectra."""
        vrtspec, divspec = self.getvrtdivspec(ugrid, vgrid, ntrunc=ntrunc)
        return (
            self._invlapspec(vrtspec, "getpsichispec"),
            self._invlapspec(divspec, "getpsichispec"),
        )

    def getpsichi_component(
        self,
        ugrid: np.ndarray,
        vgrid: np.ndarray,
        ntrunc: Optional[int],
        component: str,
        operation_name: str,
    ) -> np.ndarray:
        """Compute one streamfunction/velocity-potential packed grid field."""
        fieldspec = self.getpsichi_spec_component(
            ugrid,
            vgrid,
            ntrunc=ntrunc,
            component=component,
            operation_name=operation_name,
        )
        return self.spectogrd(fieldspec)

    def getpsi(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Compute streamfunction from packed reduced-grid wind."""
        return self.getpsichi_component(
            ugrid,
            vgrid,
            ntrunc=ntrunc,
            component="psi",
            operation_name="getpsi",
        )

    def getchi(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> np.ndarray:
        """Compute velocity potential from packed reduced-grid wind."""
        return self.getpsichi_component(
            ugrid,
            vgrid,
            ntrunc=ntrunc,
            component="chi",
            operation_name="getchi",
        )

    def getpsichi(
        self, ugrid: np.ndarray, vgrid: np.ndarray, ntrunc: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute streamfunction and velocity potential from packed reduced-grid wind."""
        psispec, chispec = self.getpsichispec(ugrid, vgrid, ntrunc=ntrunc)
        return self._spectogrd_pair(psispec, chispec)

    def specsmooth(self, datagrid: np.ndarray, smooth: np.ndarray) -> np.ndarray:
        """Apply isotropic spectral smoothing on the reduced-grid spectral chain."""
        self._validate_packed_grid_data(
            datagrid,
            "specsmooth",
        )
        dataspec = self.grdtospec(datagrid, self.ndgl - 1)
        smoothed_spec = self._apply_spectral_smoothing(
            dataspec,
            smooth,
            operation_name="specsmooth",
        )
        return self.spectogrd(smoothed_spec)


def regrid(
    grdin: ReducedGaussianSpharmt,
    grdout: ReducedGaussianSpharmt,
    datagrid: np.ndarray,
    ntrunc: Optional[int] = None,
    smooth: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Regrid packed reduced-Gaussian scalar data through spectral space."""
    grdin._validate_packed_grid_data(datagrid, "regrid")

    max_ntrunc = min(grdin.ndgl - 1, grdout.ndgl - 1)
    if ntrunc is None:
        ntrunc = max_ntrunc
    elif ntrunc < 0 or ntrunc > max_ntrunc:
        raise ReducedGaussianValidationError(
            f"regrid ntrunc must be between 0 and {max_ntrunc}, got {ntrunc}"
        )

    dataspec = grdin.grdtospec(datagrid, ntrunc=ntrunc)
    if smooth is not None:
        dataspec = grdin._apply_spectral_smoothing(
            dataspec,
            smooth,
            operation_name="regrid",
            expected_size=grdout.ndgl,
        )
    return grdout.spectogrd(dataspec)


def regriduv(
    grdin: ReducedGaussianSpharmt,
    grdout: ReducedGaussianSpharmt,
    ugrid: np.ndarray,
    vgrid: np.ndarray,
    ntrunc: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Regrid packed reduced-Gaussian wind data through vrt/div spectra."""
    ugrid_array = np.asarray(ugrid)
    vgrid_array = np.asarray(vgrid)
    if ugrid_array.shape != vgrid_array.shape:
        raise ReducedGaussianValidationError(
            "regriduv needs ugrid and vgrid with the same shape"
        )
    grdin._validate_packed_grid_data(ugrid_array, "regriduv")

    max_ntrunc = min(grdin.ndgl - 1, grdout.ndgl - 1)
    if ntrunc is None:
        ntrunc = max_ntrunc
    elif ntrunc < 0 or ntrunc > max_ntrunc:
        raise ReducedGaussianValidationError(
            f"regriduv ntrunc must be between 0 and {max_ntrunc}, got {ntrunc}"
        )

    vrtspec, divspec = grdin.getvrtdivspec(ugrid_array, vgrid_array, ntrunc=ntrunc)
    return grdout.getuv(vrtspec, divspec)
