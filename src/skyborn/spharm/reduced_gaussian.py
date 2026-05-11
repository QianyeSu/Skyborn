"""Packed reduced-Gaussian spherical harmonic transforms."""

from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from . import _spherepack
from .spherical_harmonics import (
    DEFAULT_EARTH_RADIUS,
    ComplexArray,
    FloatArray,
    SpheremackError,
    ValidationError,
    gaussian_lats_wts,
)

IntArray = NDArray[np.integer]

_private_vars = (
    "pl",
    "nlat",
    "npoints",
    "rsphere",
    "gridtype",
    "legfunc",
    "precision",
)
_GLOBAL_BASIS_CACHE: Dict[Tuple[int, int], FloatArray] = {}
_GLOBAL_DBASIS_CACHE: Dict[Tuple[int, int], FloatArray] = {}
_GLOBAL_BASIS_TRANSPOSED_CACHE: Dict[Tuple[int, int], FloatArray] = {}
_GLOBAL_DBASIS_TRANSPOSED_CACHE: Dict[Tuple[int, int], FloatArray] = {}
_MAX_GLOBAL_BASIS_ENTRIES = 2
_MAX_GLOBAL_TRANSPOSED_BASIS_ENTRIES = 1


def _cache_basis_array(
    cache: Dict[Tuple[int, int], FloatArray],
    key: Tuple[int, int],
    value: FloatArray,
    max_entries: int = _MAX_GLOBAL_BASIS_ENTRIES,
) -> None:
    if key in cache:
        return
    if len(cache) >= max_entries:
        cache.pop(next(iter(cache)))
    cache[key] = value


class ReducedGaussianSpharmt:
    """
    Scalar spherical harmonic transforms for packed reduced Gaussian grids.

    The horizontal grid is represented by a ``pl`` vector and row-packed values
    with shape ``(sum(pl), ...)``.  Spectral coefficients use the same public
    triangular layout as :class:`skyborn.spharm.Spharmt`, so scalar spectra can
    be passed through the existing ``spharm`` spectral utilities without first
    interpolating the reduced grid to a rectangular Gaussian grid.

    Attributes:
        pl: Reduced Gaussian row lengths ordered north-to-south.
        nlat: Number of latitude rows.
        npoints: Total packed grid-point count, equal to ``sum(pl)``.
        rsphere: Sphere radius in meters.
        gridtype: Fixed grid type label ``"reduced_gaussian"``.
        legfunc: Legendre function handling mode - ``"stored"`` or ``"computed"``.
        precision: Public output precision mode - ``"auto"``, ``"single"``, or
            ``"double"``.
    """

    def __init__(
        self,
        pl: IntArray,
        rsphere: float = DEFAULT_EARTH_RADIUS,
        legfunc: str = "stored",
        precision: str = "auto",
    ) -> None:
        """
        Create a reduced-Gaussian scalar transform backend.

        Args:
            pl: Reduced Gaussian row lengths ordered north-to-south. Each entry
                gives the number of longitude points on that latitude row.
            rsphere: Sphere radius in meters.
            legfunc: Legendre function handling mode - ``"stored"`` or
                ``"computed"``.
            precision: Public output precision mode. ``"auto"`` follows the
                highest precision seen in relevant inputs, ``"single"`` returns
                float32/complex64 outputs, and ``"double"`` returns
                float64/complex128 outputs.

        Raises:
            ValidationError: If ``pl``, ``rsphere``, ``legfunc``, or
                ``precision`` is invalid.
        """
        self.pl = self._validate_pl(pl)
        self.nlat = int(self.pl.size)
        self.npoints = int(self.pl.sum())
        self.rsphere = float(rsphere)
        self.gridtype = "reduced_gaussian"
        self.legfunc = self._validate_legfunc(legfunc)
        self.precision = self._validate_precision(precision)

        if self.rsphere <= 0.0:
            raise ValidationError(f"rsphere must be positive, got {rsphere}")

        latitudes, weights = gaussian_lats_wts(self.nlat)
        self._weights = np.asfortranarray(weights.astype(np.float32, copy=False))
        self._sin_theta = np.asfortranarray(
            np.cos(np.deg2rad(latitudes)).astype(np.float32, copy=False)
        )
        self._basis_cache: Dict[int, FloatArray] = {}
        self._dbasis_cache: Dict[int, FloatArray] = {}
        self._basis_transposed_cache: Dict[int, FloatArray] = {}
        self._dbasis_transposed_cache: Dict[int, FloatArray] = {}

    @staticmethod
    def _validate_pl(pl: IntArray) -> IntArray:
        arr = np.asarray(pl, dtype=np.int32)
        if arr.ndim != 1:
            raise ValidationError(f"pl must be rank 1, got shape {arr.shape}")
        if arr.size < 3:
            raise ValidationError(f"pl must contain at least 3 rows, got {arr.size}")
        if np.any(arr < 4):
            raise ValidationError("each reduced Gaussian row must contain >= 4 points")
        return np.ascontiguousarray(arr)

    @staticmethod
    def _validate_legfunc(legfunc: str) -> str:
        if legfunc not in ("stored", "computed"):
            raise ValidationError(
                f'legfunc must be "stored" or "computed", got "{legfunc}"'
            )
        return legfunc

    @staticmethod
    def _validate_precision(precision: str) -> str:
        if precision not in ("auto", "single", "double"):
            raise ValidationError(
                f'precision must be "auto", "single", or "double", got "{precision}"'
            )
        return precision

    def _public_real_dtype(self, *arrays: np.ndarray) -> np.dtype:
        if self.precision == "single":
            return np.dtype(np.float32)
        if self.precision == "double":
            return np.dtype(np.float64)
        for array in arrays:
            dtype = np.dtype(np.asarray(array).dtype)
            if dtype in (np.dtype(np.float64), np.dtype(np.complex128)):
                return np.dtype(np.float64)
        return np.dtype(np.float32)

    def _public_complex_dtype(self, *arrays: np.ndarray) -> np.dtype:
        return (
            np.dtype(np.complex128)
            if self._public_real_dtype(*arrays) == np.dtype(np.float64)
            else np.dtype(np.complex64)
        )

    def __setattr__(self, key: str, val) -> None:
        if key in self.__dict__ and key in _private_vars:
            raise AttributeError(f"Attempt to rebind read-only instance variable {key}")
        self.__dict__[key] = val

    def __delattr__(self, key: str) -> None:
        if key in self.__dict__ and key in _private_vars:
            raise AttributeError(f"Attempt to unbind read-only instance variable {key}")
        del self.__dict__[key]

    def __repr__(self) -> str:
        return (
            f"ReducedGaussianSpharmt(nlat={self.nlat:d}, "
            f"npoints={self.npoints:d}, rsphere={self.rsphere:e}, "
            f"precision='{self.precision}')"
        )

    def close(self) -> None:
        """Compatibility no-op for callers that manage transform lifetimes."""

    def _validate_ntrunc(self, ntrunc: Optional[int]) -> int:
        if ntrunc is None:
            return self.nlat - 1
        ntrunc = int(ntrunc)
        if ntrunc < 0 or ntrunc > self.nlat - 1:
            raise ValidationError(
                f"ntrunc must be between 0 and {self.nlat - 1}, got {ntrunc}"
            )
        return ntrunc

    def _basis(self, ntrunc: int) -> FloatArray:
        basis = self._basis_cache.get(ntrunc)
        if basis is not None:
            return basis

        key = (self.nlat, ntrunc)
        basis = _GLOBAL_BASIS_CACHE.get(key)
        if basis is not None:
            self._basis_cache[ntrunc] = basis
            return basis

        try:
            basis, ierror = _spherepack.reduced_gaussian_legendre_basis(
                self.nlat, ntrunc
            )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian Legendre setup failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                f"reduced Gaussian Legendre setup failed with error code {ierror}"
            )

        basis = np.asfortranarray(np.asarray(basis, dtype=np.float32))
        _cache_basis_array(_GLOBAL_BASIS_CACHE, key, basis)
        self._basis_cache[ntrunc] = basis
        return basis

    def _dbasis(self, ntrunc: int) -> FloatArray:
        dbasis = self._dbasis_cache.get(ntrunc)
        if dbasis is not None:
            return dbasis

        key = (self.nlat, ntrunc)
        dbasis = _GLOBAL_DBASIS_CACHE.get(key)
        if dbasis is not None:
            self._dbasis_cache[ntrunc] = dbasis
            return dbasis

        try:
            basis = self._basis(ntrunc)
            if hasattr(_spherepack, "reduced_gaussian_legendre_derivative_from_basis"):
                dbasis, ierror = (
                    _spherepack.reduced_gaussian_legendre_derivative_from_basis(
                        basis, ntrunc
                    )
                )
            else:
                dbasis, ierror = _spherepack.reduced_gaussian_legendre_derivative_basis(
                    self.nlat, ntrunc
                )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian Legendre derivative setup failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                "reduced Gaussian Legendre derivative setup failed with "
                f"error code {ierror}"
            )

        dbasis = np.asfortranarray(np.asarray(dbasis, dtype=np.float32))
        _cache_basis_array(_GLOBAL_DBASIS_CACHE, key, dbasis)
        self._dbasis_cache[ntrunc] = dbasis
        return dbasis

    def _basis_transposed(self, ntrunc: int) -> FloatArray:
        basis = self._basis_transposed_cache.get(ntrunc)
        if basis is not None:
            return basis

        key = (self.nlat, ntrunc)
        basis = _GLOBAL_BASIS_TRANSPOSED_CACHE.get(key)
        if basis is None:
            basis = np.asfortranarray(self._basis(ntrunc).T)
            _cache_basis_array(
                _GLOBAL_BASIS_TRANSPOSED_CACHE,
                key,
                basis,
                max_entries=_MAX_GLOBAL_TRANSPOSED_BASIS_ENTRIES,
            )
        self._basis_transposed_cache[ntrunc] = basis
        return basis

    def _dbasis_transposed(self, ntrunc: int) -> FloatArray:
        dbasis = self._dbasis_transposed_cache.get(ntrunc)
        if dbasis is not None:
            return dbasis

        key = (self.nlat, ntrunc)
        dbasis = _GLOBAL_DBASIS_TRANSPOSED_CACHE.get(key)
        if dbasis is None:
            dbasis = np.asfortranarray(self._dbasis(ntrunc).T)
            _cache_basis_array(
                _GLOBAL_DBASIS_TRANSPOSED_CACHE,
                key,
                dbasis,
                max_entries=_MAX_GLOBAL_TRANSPOSED_BASIS_ENTRIES,
            )
        self._dbasis_transposed_cache[ntrunc] = dbasis
        return dbasis

    def _validate_grid_data(
        self, data: FloatArray, operation_name: str
    ) -> Tuple[int, FloatArray, Tuple[int, ...]]:
        data = np.asarray(data)
        if data.ndim < 1:
            raise ValidationError(
                f"{operation_name} needs at least a rank 1 packed array, "
                f"got rank {data.ndim}"
            )
        if data.shape[0] != self.npoints:
            raise ValidationError(
                f"{operation_name} needs a packed first dimension of "
                f"{self.npoints}, got {data.shape[0]}"
            )

        if data.ndim == 1:
            nt = 1
            extra_shape: Tuple[int, ...] = ()
            normalized = data.reshape(self.npoints, 1)
        else:
            extra_shape = tuple(data.shape[1:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized = data.reshape(self.npoints, nt)

        return nt, np.asfortranarray(normalized, dtype=np.float32), extra_shape

    def _validate_spectral_data(
        self, data: ComplexArray, operation_name: str
    ) -> Tuple[int, int, ComplexArray, Tuple[int, ...]]:
        data = np.asarray(data)
        if data.ndim < 1:
            raise ValidationError(
                f"{operation_name} needs at least a rank 1 spectrum, "
                f"got rank {data.ndim}"
            )

        ntrunc = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - data.shape[0])))
        expected = (ntrunc + 1) * (ntrunc + 2) // 2
        if expected != data.shape[0]:
            raise ValidationError(
                f"{operation_name} got invalid triangular spectral size "
                f"{data.shape[0]}"
            )
        if ntrunc > self.nlat - 1:
            raise ValidationError(
                f"ntrunc too large - can be max of {self.nlat - 1}, got {ntrunc}"
            )

        if data.ndim == 1:
            nt = 1
            extra_shape: Tuple[int, ...] = ()
            normalized = data.reshape(data.shape[0], 1)
        else:
            extra_shape = tuple(data.shape[1:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized = data.reshape(data.shape[0], nt)

        return (
            nt,
            ntrunc,
            np.asfortranarray(normalized, dtype=np.complex64),
            extra_shape,
        )

    def _restore_spectral_shape(
        self,
        dataspec: ComplexArray,
        extra_shape: Tuple[int, ...],
        output_dtype: Optional[np.dtype] = None,
    ) -> ComplexArray:
        if not extra_shape:
            result = dataspec[:, 0]
        else:
            result = dataspec.reshape((dataspec.shape[0],) + extra_shape)
        if output_dtype is None:
            if self.precision == "double":
                output_dtype = np.complex128
            elif self.precision == "single":
                output_dtype = np.complex64
        if output_dtype is None:
            return result
        return np.asarray(result).astype(np.dtype(output_dtype), copy=False)

    def _restore_grid_shape(
        self,
        datagrid: FloatArray,
        extra_shape: Tuple[int, ...],
        output_dtype: Optional[np.dtype] = None,
    ) -> FloatArray:
        if not extra_shape:
            result = datagrid[:, 0]
        else:
            result = datagrid.reshape((self.npoints,) + extra_shape)
        if output_dtype is None:
            if self.precision == "double":
                output_dtype = np.float64
            elif self.precision == "single":
                output_dtype = np.float32
        if output_dtype is None:
            return result
        return np.asarray(result).astype(np.dtype(output_dtype), copy=False)

    def grdtospec(
        self, datagrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """
        Analyze packed reduced Gaussian scalar data into ``spharm`` spectra.

        Parameters
        ----------
        datagrid:
            Packed reduced-grid values with shape ``(sum(pl), ...)``.
        ntrunc:
            Optional triangular truncation.  Defaults to ``nlat - 1``.
        """
        _, normalized, extra_shape = self._validate_grid_data(datagrid, "grdtospec")
        ntrunc = self._validate_ntrunc(ntrunc)
        basis = self._basis(ntrunc)

        try:
            dataspec, ierror = _spherepack.reduced_gaussian_grdtospec(
                normalized,
                self.pl,
                self._weights,
                basis,
                ntrunc,
            )
        except Exception as exc:
            raise SpheremackError(f"reduced Gaussian grdtospec failed: {exc}") from exc

        if ierror != 0:
            raise SpheremackError(
                f"reduced Gaussian grdtospec failed with error code {ierror}"
            )

        return self._restore_spectral_shape(
            dataspec, extra_shape, self._public_complex_dtype(datagrid)
        )

    def spectogrd(self, dataspec: ComplexArray) -> FloatArray:
        """
        Synthesize ``spharm`` spectra onto the packed reduced Gaussian grid.
        """
        _, ntrunc, normalized, extra_shape = self._validate_spectral_data(
            dataspec, "spectogrd"
        )
        basis = self._basis_transposed(ntrunc)

        try:
            datagrid, ierror = _spherepack.reduced_gaussian_spectogrd(
                normalized,
                self.pl,
                basis,
                ntrunc,
                self.npoints,
            )
        except Exception as exc:
            raise SpheremackError(f"reduced Gaussian spectogrd failed: {exc}") from exc

        if ierror != 0:
            raise SpheremackError(
                f"reduced Gaussian spectogrd failed with error code {ierror}"
            )

        return self._restore_grid_shape(
            datagrid, extra_shape, self._public_real_dtype(dataspec)
        )

    def _spectogrd_pair(
        self, spec_a: ComplexArray, spec_b: ComplexArray
    ) -> Tuple[FloatArray, FloatArray]:
        """Synthesize two scalar spectra together through one batched call."""
        _, ntrunc, normalized_a, normalized_b, extra_shape = (
            self._validate_paired_spectral_data(spec_a, spec_b, "_spectogrd_pair")
        )

        # On the user's O320 single-layer workload (nt=1), the current native
        # paired scalar synthesis kernel is slower than two independent
        # single-field syntheses. Keep the paired kernel for multi-field
        # workloads, but prefer the simpler path for a single slice until the
        # Fortran paired kernel is reworked.
        if normalized_a.shape[1] == 1:
            return self.spectogrd(spec_a), self.spectogrd(spec_b)

        if hasattr(_spherepack, "reduced_gaussian_spectogrd_pair"):
            basis = self._basis_transposed(ntrunc)
            try:
                grid_a, grid_b, ierror = _spherepack.reduced_gaussian_spectogrd_pair(
                    normalized_a,
                    normalized_b,
                    self.pl,
                    basis,
                    ntrunc,
                    self.npoints,
                )
            except Exception as exc:
                raise SpheremackError(
                    f"reduced Gaussian paired spectogrd failed: {exc}"
                ) from exc

            if ierror != 0:
                raise SpheremackError(
                    "reduced Gaussian paired spectogrd failed with "
                    f"error code {ierror}"
                )

            return (
                self._restore_grid_shape(
                    grid_a, extra_shape, self._public_real_dtype(spec_a, spec_b)
                ),
                self._restore_grid_shape(
                    grid_b, extra_shape, self._public_real_dtype(spec_a, spec_b)
                ),
            )

        paired_specs = np.stack((spec_a, spec_b), axis=-1)
        paired_grids = self.spectogrd(paired_specs)
        return paired_grids[..., 0], paired_grids[..., 1]

    def _validate_paired_grid_data(
        self, ugrid: FloatArray, vgrid: FloatArray, operation_name: str
    ) -> Tuple[int, FloatArray, FloatArray, Tuple[int, ...]]:
        if np.shape(ugrid) != np.shape(vgrid):
            raise ValidationError("ugrid and vgrid must have the same shape")

        nt, normalized_u, extra_shape = self._validate_grid_data(ugrid, operation_name)
        _, normalized_v, v_extra_shape = self._validate_grid_data(vgrid, operation_name)
        if extra_shape != v_extra_shape:
            raise ValidationError("ugrid and vgrid must have consistent dimensions")
        return nt, normalized_u, normalized_v, extra_shape

    def _validate_paired_spectral_data(
        self, spec_a: ComplexArray, spec_b: ComplexArray, operation_name: str
    ) -> Tuple[int, int, ComplexArray, ComplexArray, Tuple[int, ...]]:
        if np.shape(spec_a) != np.shape(spec_b):
            raise ValidationError("paired spectra must have the same shape")

        nt_a, ntrunc_a, normalized_a, extra_shape = self._validate_spectral_data(
            spec_a, operation_name
        )
        nt_b, ntrunc_b, normalized_b, b_extra_shape = self._validate_spectral_data(
            spec_b, operation_name
        )
        if nt_a != nt_b or ntrunc_a != ntrunc_b or extra_shape != b_extra_shape:
            raise ValidationError("paired spectra must have consistent dimensions")
        return nt_a, ntrunc_a, normalized_a, normalized_b, extra_shape

    def getvrtdivspec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> Tuple[ComplexArray, ComplexArray]:
        """
        Compute vorticity and divergence spectra from packed reduced-grid winds.
        """
        _, normalized_u, normalized_v, extra_shape = self._validate_paired_grid_data(
            ugrid, vgrid, "getvrtdivspec"
        )
        ntrunc = self._validate_ntrunc(ntrunc)
        basis = self._basis_transposed(ntrunc)
        dbasis = self._dbasis_transposed(ntrunc)

        try:
            vrtspec, divspec, ierror = _spherepack.reduced_gaussian_getvrtdivspec(
                normalized_u,
                normalized_v,
                self.pl,
                self._weights,
                basis,
                dbasis,
                self._sin_theta,
                ntrunc,
                self.rsphere,
            )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian vector analysis failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                "reduced Gaussian vector analysis failed with " f"error code {ierror}"
            )

        return (
            self._restore_spectral_shape(vrtspec, extra_shape),
            self._restore_spectral_shape(divspec, extra_shape),
        )

    def getvrtspec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """Compute only vorticity spectral coefficients from packed winds."""
        _, normalized_u, normalized_v, extra_shape = self._validate_paired_grid_data(
            ugrid, vgrid, "getvrtspec"
        )
        ntrunc = self._validate_ntrunc(ntrunc)
        basis = self._basis_transposed(ntrunc)
        dbasis = self._dbasis_transposed(ntrunc)

        try:
            vrtspec, ierror = _spherepack.reduced_gaussian_getvrtspec(
                normalized_u,
                normalized_v,
                self.pl,
                self._weights,
                basis,
                dbasis,
                self._sin_theta,
                ntrunc,
                self.rsphere,
            )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian vorticity analysis failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                "reduced Gaussian vorticity analysis failed with "
                f"error code {ierror}"
            )

        return self._restore_spectral_shape(
            vrtspec, extra_shape, self._public_complex_dtype(ugrid, vgrid)
        )

    def getdivspec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """Compute only divergence spectral coefficients from packed winds."""
        _, normalized_u, normalized_v, extra_shape = self._validate_paired_grid_data(
            ugrid, vgrid, "getdivspec"
        )
        ntrunc = self._validate_ntrunc(ntrunc)
        basis = self._basis_transposed(ntrunc)
        dbasis = self._dbasis_transposed(ntrunc)

        try:
            divspec, ierror = _spherepack.reduced_gaussian_getdivspec(
                normalized_u,
                normalized_v,
                self.pl,
                self._weights,
                basis,
                dbasis,
                self._sin_theta,
                ntrunc,
                self.rsphere,
            )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian divergence analysis failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                "reduced Gaussian divergence analysis failed with "
                f"error code {ierror}"
            )

        return self._restore_spectral_shape(
            divspec, extra_shape, self._public_complex_dtype(ugrid, vgrid)
        )

    def getuv(
        self, vrtspec: ComplexArray, divspec: ComplexArray
    ) -> Tuple[FloatArray, FloatArray]:
        """
        Synthesize packed reduced-grid winds from vorticity/divergence spectra.
        """
        _, _, normalized_vrt, normalized_div, extra_shape = (
            self._validate_paired_spectral_data(vrtspec, divspec, "getuv")
        )

        psispec = self._restore_spectral_shape(
            _spherepack.invlap(normalized_vrt, self.rsphere), extra_shape
        )
        chispec = self._restore_spectral_shape(
            _spherepack.invlap(normalized_div, self.rsphere), extra_shape
        )
        u_chi, v_chi, v_psi, u_psi = self.getgrad_pair(chispec, psispec)

        return (
            u_chi - u_psi,
            v_chi + v_psi,
        )

    def getgrad(self, dataspec: ComplexArray) -> Tuple[FloatArray, FloatArray]:
        """Compute the packed reduced-grid vector gradient of a scalar spectrum."""
        _, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, "getgrad"
        )
        basis = self._basis_transposed(ntrunc)
        dbasis = self._dbasis_transposed(ntrunc)

        try:
            if hasattr(_spherepack, "reduced_gaussian_getgrad"):
                ugrad, vgrad, ierror = _spherepack.reduced_gaussian_getgrad(
                    normalized_spec,
                    self.pl,
                    basis,
                    dbasis,
                    self._sin_theta,
                    ntrunc,
                    self.npoints,
                    self.rsphere,
                )
            else:
                zero_spec = np.zeros_like(normalized_spec)
                ugrad, vgrad, _, _, ierror = _spherepack.reduced_gaussian_getgrad_pair(
                    normalized_spec,
                    zero_spec,
                    self.pl,
                    basis,
                    dbasis,
                    self._sin_theta,
                    ntrunc,
                    self.npoints,
                    self.rsphere,
                )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian gradient synthesis failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                "reduced Gaussian gradient synthesis failed with "
                f"error code {ierror}"
            )

        return (
            self._restore_grid_shape(ugrad, extra_shape),
            self._restore_grid_shape(vgrad, extra_shape),
        )

    def getgrad_pair(
        self, spec_a: ComplexArray, spec_b: ComplexArray
    ) -> Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
        """
        Compute packed reduced-grid gradients for two scalar spectra at once.

        The paired path shares Legendre-basis reads and longitude synthesis
        recurrences, which is faster than two independent ``getgrad`` calls for
        Helmholtz-style workflows.
        """
        _, ntrunc, normalized_a, normalized_b, extra_shape = (
            self._validate_paired_spectral_data(spec_a, spec_b, "getgrad_pair")
        )
        basis = self._basis_transposed(ntrunc)
        dbasis = self._dbasis_transposed(ntrunc)

        try:
            a_ugrad, a_vgrad, b_ugrad, b_vgrad, ierror = (
                _spherepack.reduced_gaussian_getgrad_pair(
                    normalized_a,
                    normalized_b,
                    self.pl,
                    basis,
                    dbasis,
                    self._sin_theta,
                    ntrunc,
                    self.npoints,
                    self.rsphere,
                )
            )
        except Exception as exc:
            raise SpheremackError(
                f"reduced Gaussian paired gradient synthesis failed: {exc}"
            ) from exc

        if ierror != 0:
            raise SpheremackError(
                "reduced Gaussian paired gradient synthesis failed with "
                f"error code {ierror}"
            )

        public_dtype = self._public_real_dtype(spec_a, spec_b)
        return (
            self._restore_grid_shape(a_ugrad, extra_shape, public_dtype),
            self._restore_grid_shape(a_vgrad, extra_shape, public_dtype),
            self._restore_grid_shape(b_ugrad, extra_shape, public_dtype),
            self._restore_grid_shape(b_vgrad, extra_shape, public_dtype),
        )

    def _invlapspec(
        self, dataspec: ComplexArray, operation_name: str = "_invlapspec"
    ) -> ComplexArray:
        """Apply the inverse Laplacian to a scalar spectrum and preserve shape."""
        _, _, normalized, extra_shape = self._validate_spectral_data(
            dataspec, operation_name
        )
        invlap_spec = _spherepack.invlap(normalized, self.rsphere)
        return self._restore_spectral_shape(invlap_spec, extra_shape)

    def specsmooth(self, datagrid: FloatArray, smooth: FloatArray) -> FloatArray:
        """Apply latitude-dependent spectral smoothing to a packed scalar field."""
        smooth = np.asarray(smooth)
        if smooth.ndim != 1 or smooth.shape[0] != self.nlat:
            raise ValidationError(
                f"smooth must be rank 1 with size {self.nlat}, got shape {smooth.shape}"
            )

        _, _, normalized_spec, extra_shape = self._validate_spectral_data(
            self.grdtospec(datagrid), "specsmooth"
        )
        smoothed = _spherepack.multsmoothfact(
            normalized_spec, smooth.astype(np.float32)
        )
        return self.spectogrd(self._restore_spectral_shape(smoothed, extra_shape))
