"""Thin Python wrapper around the experimental native ectrans backend."""

from __future__ import annotations

import importlib
from collections import OrderedDict

import numpy as np

_backend = importlib.import_module("skyborn.spharm._ectrans_backend")


class UvToVordivBlockSolver:
    """Thin handle-owning wrapper around the cached native block solver."""

    def __init__(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
        rsphere: float,
    ) -> None:
        ntrunc, ncoeff, _ = _uv_to_vordiv_geometry(uspec, vspec, rsphere)
        self._ntrunc = ntrunc
        self._ncoeff = ncoeff
        self._rsphere = float(rsphere)
        self._handle = create_uv_to_vordiv_block_setup(uspec, vspec, self._rsphere)
        self._closed = False

    @property
    def rsphere(self) -> float:
        """Return the sphere radius used to build the cached solver."""
        return self._rsphere

    @property
    def ntrunc(self) -> int:
        """Return the spectral truncation this solver was built for."""
        return self._ntrunc

    @property
    def ncoeff(self) -> int:
        """Return the public triangular coefficient count this solver was built for."""
        return self._ncoeff

    @property
    def closed(self) -> bool:
        """Whether the native setup handle has been closed."""
        return self._closed

    def solve(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Apply the cached native block solver to one public ``u/v`` spectrum pair."""
        if self._closed:
            raise RuntimeError("uv_to_vordiv block solver is closed")
        ntrunc, ncoeff, _ = _uv_to_vordiv_geometry(uspec, vspec, self._rsphere)
        if ntrunc != self._ntrunc or ncoeff != self._ncoeff:
            raise ValueError(
                "uv_to_vordiv block solver geometry mismatch: "
                f"expected ntrunc={self._ntrunc}, ncoeff={self._ncoeff}; "
                f"got ntrunc={ntrunc}, ncoeff={ncoeff}"
            )
        return uv_to_vordiv_block_native_with_setup(self._handle, uspec, vspec)

    def close(self) -> None:
        """Close the native setup handle."""
        if not self._closed:
            destroy_uv_to_vordiv_block_setup(self._handle)
            self._closed = True
            self._handle = None

    def __enter__(self) -> "UvToVordivBlockSolver":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


class UvToVordivBlockSolverCache:
    """Small internal cache for reusing block solvers by spectral geometry."""

    def __init__(self, max_entries: int | None = None) -> None:
        if max_entries is not None and max_entries < 1:
            raise ValueError(f"max_entries must be positive or None, got {max_entries}")
        self._solvers: OrderedDict[tuple[int, float], UvToVordivBlockSolver] = (
            OrderedDict()
        )
        self._max_entries = max_entries
        self._closed = False

    @property
    def closed(self) -> bool:
        """Whether this cache has been closed."""
        return self._closed

    @property
    def size(self) -> int:
        """Return the number of live cached solver entries."""
        return len(self._solvers)

    @property
    def max_entries(self) -> int | None:
        """Return the configured cache entry limit."""
        return self._max_entries

    @staticmethod
    def _cache_key(
        uspec: np.ndarray, vspec: np.ndarray, rsphere: float
    ) -> tuple[int, float]:
        ntrunc, _, normalized_rsphere = _uv_to_vordiv_geometry(uspec, vspec, rsphere)
        return (ntrunc, normalized_rsphere)

    def get_solver(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
        rsphere: float,
    ) -> UvToVordivBlockSolver:
        """Return a cached solver for one spectral geometry, creating it if needed."""
        if self._closed:
            raise RuntimeError("uv_to_vordiv block solver cache is closed")
        key = self._cache_key(uspec, vspec, rsphere)
        solver = self._solvers.get(key)
        if solver is None or solver.closed:
            if solver is not None and solver.closed:
                self._solvers.pop(key, None)
            solver = UvToVordivBlockSolver(uspec, vspec, rsphere)
            self._solvers[key] = solver
            if self._max_entries is not None and len(self._solvers) > self._max_entries:
                evict_key, evict_solver = self._solvers.popitem(last=False)
                if evict_key != key:
                    evict_solver.close()
                else:
                    self._solvers[key] = solver
        else:
            self._solvers.move_to_end(key)
        return solver

    def solve(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
        rsphere: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve through a cached solver keyed by spectral shape and ``rsphere``."""
        solver = self.get_solver(uspec, vspec, rsphere)
        return solver.solve(uspec, vspec)

    def clear(self) -> None:
        """Close and drop all cached solver entries while keeping the cache reusable."""
        for solver in self._solvers.values():
            solver.close()
        self._solvers.clear()

    def close(self) -> None:
        """Close the cache and all cached solver entries."""
        if not self._closed:
            self.clear()
            self._closed = True

    def __enter__(self) -> "UvToVordivBlockSolverCache":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


_UV_TO_VORDIV_BLOCK_NATIVE_PILOT_CACHE = UvToVordivBlockSolverCache(max_entries=4)


def _get_uv_to_vordiv_block_native_pilot_cache() -> UvToVordivBlockSolverCache:
    """Return the module-level pilot cache, recreating it if it was closed."""
    global _UV_TO_VORDIV_BLOCK_NATIVE_PILOT_CACHE
    if _UV_TO_VORDIV_BLOCK_NATIVE_PILOT_CACHE.closed:
        _UV_TO_VORDIV_BLOCK_NATIVE_PILOT_CACHE = UvToVordivBlockSolverCache(
            max_entries=4
        )
    return _UV_TO_VORDIV_BLOCK_NATIVE_PILOT_CACHE


def _clear_uv_to_vordiv_block_native_pilot_cache() -> None:
    """Clear the module-level pilot cache without closing the cache object."""
    _get_uv_to_vordiv_block_native_pilot_cache().clear()


class UvToVordivPoa1BlockSolver:
    """Cached Python-side inverse for the narrowed PRFI/POA1 boundary pilot."""

    def __init__(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
        rsphere: float,
    ) -> None:
        ntrunc, ncoeff, _ = _uv_to_vordiv_geometry(uspec, vspec, rsphere)
        self._ntrunc = ntrunc
        self._ncoeff = ncoeff
        self._rsphere = float(rsphere)
        self._nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2)
        self._blocks: list[tuple[np.ndarray, np.ndarray]] = []
        self._estimated_bytes = 0
        self._closed = False
        self._build_cached_blocks()

    @property
    def rsphere(self) -> float:
        """Return the sphere radius used to build this solver."""
        return self._rsphere

    @property
    def ntrunc(self) -> int:
        """Return the spectral truncation this solver was built for."""
        return self._ntrunc

    @property
    def ncoeff(self) -> int:
        """Return the public triangular coefficient count this solver was built for."""
        return self._ncoeff

    @property
    def estimated_bytes(self) -> int:
        """Approximate live cache footprint for the stored per-block operators."""
        return self._estimated_bytes

    @property
    def closed(self) -> bool:
        """Whether this solver has been closed."""
        return self._closed

    def _build_cached_blocks(self) -> None:
        for m in range(self._ntrunc + 1):
            operator = np.zeros((self._nlei1 * 4, self._nlei1 * 4), dtype=np.float64)

            for basis_idx in range(self._nlei1 * 4):
                poa_basis = np.zeros((self._nlei1, 4), dtype=np.float64)
                poa_basis.reshape(-1)[basis_idx] = 1.0
                vrt_basis, div_basis, ierr = poa1_to_vordiv(
                    self._ntrunc,
                    m,
                    poa_basis,
                    self._rsphere,
                )
                if ierr != 0:
                    raise RuntimeError(
                        f"poa1_to_vordiv failed while building m={m} POA1 basis "
                        f"at flattened index {basis_idx}"
                    )
                uspec_basis, vspec_basis, ierr = vordiv_to_uv(
                    vrt_basis,
                    div_basis,
                    self._rsphere,
                )
                if ierr != 0:
                    raise RuntimeError(
                        f"vordiv_to_uv failed while building m={m} POA1 basis "
                        f"at flattened index {basis_idx}"
                    )
                obs_basis, ierr = prfi1b_uv_block(
                    m,
                    uspec_basis / self._rsphere,
                    vspec_basis / self._rsphere,
                    self._rsphere,
                )
                if ierr != 0:
                    raise RuntimeError(
                        f"prfi1b_uv_block failed while building m={m} POA1 basis "
                        f"at flattened index {basis_idx}"
                    )
                operator[:, basis_idx] = obs_basis.reshape(-1)

            col_norms = np.linalg.norm(operator, axis=0)
            active = col_norms > 1e-18
            reduced_operator = np.ascontiguousarray(
                operator[:, active], dtype=np.float64
            )
            self._blocks.append((reduced_operator, active))
            self._estimated_bytes += reduced_operator.nbytes + active.nbytes

    def solve(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Apply the cached per-``m`` POA1 pilot inverse to one spectral batch."""
        if self._closed:
            raise RuntimeError("uv_to_vordiv POA1 block solver is closed")

        uspec_arr = np.ascontiguousarray(uspec, dtype=np.complex128)
        vspec_arr = np.ascontiguousarray(vspec, dtype=np.complex128)
        ntrunc, ncoeff, _ = _uv_to_vordiv_geometry(uspec_arr, vspec_arr, self._rsphere)
        if ntrunc != self._ntrunc or ncoeff != self._ncoeff:
            raise ValueError(
                "uv_to_vordiv POA1 block solver geometry mismatch: "
                f"expected ntrunc={self._ntrunc}, ncoeff={self._ncoeff}; "
                f"got ntrunc={ntrunc}, ncoeff={ncoeff}"
            )

        extra_shape = uspec_arr.shape[1:]
        nt = int(np.prod(extra_shape, dtype=int)) if extra_shape else 1
        uspec_flat = uspec_arr.reshape(self._ncoeff, nt)
        vspec_flat = vspec_arr.reshape(self._ncoeff, nt)
        vrt_flat = np.zeros((self._ncoeff, nt), dtype=np.complex128)
        div_flat = np.zeros((self._ncoeff, nt), dtype=np.complex128)

        for col in range(nt):
            vrt_col = np.zeros((self._ncoeff,), dtype=np.complex128)
            div_col = np.zeros((self._ncoeff,), dtype=np.complex128)
            uspec_scaled = uspec_flat[:, col] / self._rsphere
            vspec_scaled = vspec_flat[:, col] / self._rsphere

            for m, (reduced_operator, active) in enumerate(self._blocks):
                rhs_obs, ierr = prfi1b_uv_block(
                    m,
                    uspec_scaled,
                    vspec_scaled,
                    self._rsphere,
                )
                if ierr != 0:
                    raise RuntimeError(
                        f"prfi1b_uv_block failed while building m={m} rhs"
                    )
                reduced_solution, _, _, _ = np.linalg.lstsq(
                    reduced_operator,
                    rhs_obs.reshape(-1),
                    rcond=None,
                )
                poa_flat = np.zeros((self._nlei1 * 4,), dtype=np.float64)
                poa_flat[active] = reduced_solution
                vrt_block, div_block, ierr = poa1_to_vordiv(
                    self._ntrunc,
                    m,
                    poa_flat.reshape(self._nlei1, 4),
                    self._rsphere,
                )
                if ierr != 0:
                    raise RuntimeError(
                        f"poa1_to_vordiv failed while solving m={m} block"
                    )
                vrt_col += vrt_block
                div_col += div_block

            vrt_flat[:, col] = vrt_col
            div_flat[:, col] = div_col

        if extra_shape:
            return vrt_flat.reshape((self._ncoeff, *extra_shape)), div_flat.reshape(
                (self._ncoeff, *extra_shape)
            )
        return vrt_flat[:, 0], div_flat[:, 0]

    def close(self) -> None:
        """Release cached block operators."""
        if not self._closed:
            self._blocks.clear()
            self._estimated_bytes = 0
            self._closed = True

    def __enter__(self) -> "UvToVordivPoa1BlockSolver":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


class UvToVordivPoa1BlockSolverCache:
    """Small cache for the Python-side POA1 pilot solvers."""

    def __init__(self, max_entries: int | None = None) -> None:
        if max_entries is not None and max_entries < 1:
            raise ValueError(f"max_entries must be positive or None, got {max_entries}")
        self._solvers: OrderedDict[tuple[int, float], UvToVordivPoa1BlockSolver] = (
            OrderedDict()
        )
        self._max_entries = max_entries
        self._closed = False

    @property
    def closed(self) -> bool:
        """Whether this cache has been closed."""
        return self._closed

    @property
    def size(self) -> int:
        """Return the number of live cached solver entries."""
        return len(self._solvers)

    @property
    def max_entries(self) -> int | None:
        """Return the configured cache entry limit."""
        return self._max_entries

    @staticmethod
    def _cache_key(
        uspec: np.ndarray, vspec: np.ndarray, rsphere: float
    ) -> tuple[int, float]:
        ntrunc, _, normalized_rsphere = _uv_to_vordiv_geometry(uspec, vspec, rsphere)
        return (ntrunc, normalized_rsphere)

    def get_solver(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
        rsphere: float,
    ) -> UvToVordivPoa1BlockSolver:
        """Return a cached POA1 solver for one spectral geometry."""
        if self._closed:
            raise RuntimeError("uv_to_vordiv POA1 block solver cache is closed")
        key = self._cache_key(uspec, vspec, rsphere)
        solver = self._solvers.get(key)
        if solver is None or solver.closed:
            if solver is not None and solver.closed:
                self._solvers.pop(key, None)
            solver = UvToVordivPoa1BlockSolver(uspec, vspec, rsphere)
            self._solvers[key] = solver
            if self._max_entries is not None and len(self._solvers) > self._max_entries:
                evict_key, evict_solver = self._solvers.popitem(last=False)
                if evict_key != key:
                    evict_solver.close()
                else:
                    self._solvers[key] = solver
        else:
            self._solvers.move_to_end(key)
        return solver

    def solve(
        self,
        uspec: np.ndarray,
        vspec: np.ndarray,
        rsphere: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve through a cached POA1 solver keyed by spectral shape and ``rsphere``."""
        solver = self.get_solver(uspec, vspec, rsphere)
        return solver.solve(uspec, vspec)

    def clear(self) -> None:
        """Close and drop all cached solver entries while keeping the cache reusable."""
        for solver in self._solvers.values():
            solver.close()
        self._solvers.clear()

    def close(self) -> None:
        """Close the cache and all cached solver entries."""
        if not self._closed:
            self.clear()
            self._closed = True

    def __enter__(self) -> "UvToVordivPoa1BlockSolverCache":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


_UV_TO_VORDIV_POA1_BLOCK_PILOT_CACHE = UvToVordivPoa1BlockSolverCache(max_entries=1)


def _get_uv_to_vordiv_poa1_block_pilot_cache() -> UvToVordivPoa1BlockSolverCache:
    """Return the module-level POA1 pilot cache, recreating it if needed."""
    global _UV_TO_VORDIV_POA1_BLOCK_PILOT_CACHE
    if _UV_TO_VORDIV_POA1_BLOCK_PILOT_CACHE.closed:
        _UV_TO_VORDIV_POA1_BLOCK_PILOT_CACHE = UvToVordivPoa1BlockSolverCache(
            max_entries=1
        )
    return _UV_TO_VORDIV_POA1_BLOCK_PILOT_CACHE


def _clear_uv_to_vordiv_poa1_block_pilot_cache() -> None:
    """Clear the module-level POA1 pilot cache without closing the cache object."""
    _get_uv_to_vordiv_poa1_block_pilot_cache().clear()


def _infer_ntrunc_from_ncoeff(ncoeff: int) -> int:
    """Infer triangular truncation from the public spectral coefficient count."""
    ntrunc = int(-1.5 + 0.5 * np.sqrt(1.0 + 8.0 * ncoeff))
    if ((ntrunc + 1) * (ntrunc + 2)) // 2 != ncoeff:
        raise ValueError(f"invalid triangular coefficient count: {ncoeff}")
    return ntrunc


def _uv_to_vordiv_geometry(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[int, int, float]:
    """Validate one public ``u/v`` spectral pair and return its geometry key."""
    uspec_arr = np.asarray(uspec)
    vspec_arr = np.asarray(vspec)
    if uspec_arr.ndim < 1 or vspec_arr.ndim < 1:
        raise ValueError("uspec and vspec must be at least rank-1 arrays")
    if uspec_arr.shape != vspec_arr.shape:
        raise ValueError(
            f"uspec and vspec must have the same shape, got {uspec_arr.shape} and {vspec_arr.shape}"
        )
    ncoeff = int(uspec_arr.shape[0])
    ntrunc = _infer_ntrunc_from_ncoeff(ncoeff)
    return ntrunc, ncoeff, float(rsphere)


def _public_spectral_index(ntrunc: int, m: int, n: int) -> int:
    """Return the public triangular coefficient index for one ``(m, n)``."""
    idx = 0
    for mm in range(ntrunc + 1):
        for nn in range(mm, ntrunc + 1):
            if mm == m and nn == n:
                return idx
            idx += 1
    raise ValueError(f"invalid (m, n)=({m}, {n}) for ntrunc={ntrunc}")


def _public_indices_for_zonal_wavenumber(ntrunc: int, m: int) -> np.ndarray:
    """Return the public triangular coefficient indices for one zonal wavenumber."""
    indices: list[int] = []
    idx = 0
    for mm in range(ntrunc + 1):
        for _ in range(mm, ntrunc + 1):
            if mm == m:
                indices.append(idx)
            idx += 1
    return np.asarray(indices, dtype=np.intp)


def build_ledir_blocks_from_scalar_bridge(
    nlon: int,
    nlat: int,
    km: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build provisional DGEMM-only LEDIR blocks from a Gaussian bridge.

    This helper reconstructs one local-wave ``RPNMA/RPNMS`` block pair by:

    1. synthesizing one public Gaussian spectral basis mode at a time,
    2. applying a longitude FFT on the resulting grid,
    3. extracting the north-half ``km`` Fourier line into OpenIFS-style
       antisymmetric/symmetric block columns.

    The result is still a provisional bridge, not a claim of full OpenIFS
    normalization parity for every runtime case.
    """
    from .spherical_harmonics import Spharmt, gaussian_lats_wts

    if nlat < 3:
        raise ValueError(f"nlat must be at least 3, got {nlat}")
    if nlon < 4:
        raise ValueError(f"nlon must be at least 4, got {nlon}")

    ntrunc = nlat - 1
    if km < 0 or km > ntrunc:
        raise ValueError(f"require 0 <= km <= {ntrunc}, got {km}")

    ndgnh = (nlat + 1) // 2
    ia = 1 + ((ntrunc - km + 2) % 2)
    is_ = 1 + ((ntrunc - km + 1) % 2)
    ila = (ntrunc - km + 2) // 2
    ils = (ntrunc - km + 3) // 2

    sht = Spharmt(nlon=nlon, nlat=nlat, gridtype="gaussian", legfunc="stored")
    _, weights = gaussian_lats_wts(nlat)

    rpnma = np.zeros((ndgnh, ila), dtype=np.float64)
    rpnms = np.zeros((ndgnh, ils), dtype=np.float64)
    pw = np.asarray(weights[:ndgnh], dtype=np.float64)

    for j in range(ila):
        row = ia + 2 * j
        n = ntrunc + 2 - row
        if km <= n <= ntrunc:
            spec = np.zeros((nlat * (nlat + 1)) // 2, dtype=np.complex128)
            spec[_public_spectral_index(ntrunc, km, n)] = 1.0 + 0.0j
            grid = sht.spectogrd(spec)
            amp = np.fft.rfft(grid, axis=1) / float(nlon)
            rpnma[:, j] = amp[:ndgnh, km].real

    for j in range(ils):
        row = is_ + 2 * j
        n = ntrunc + 2 - row
        if km <= n <= ntrunc:
            spec = np.zeros((nlat * (nlat + 1)) // 2, dtype=np.complex128)
            spec[_public_spectral_index(ntrunc, km, n)] = 1.0 + 0.0j
            grid = sht.spectogrd(spec)
            amp = np.fft.rfft(grid, axis=1) / float(nlon)
            rpnms[:, j] = amp[:ndgnh, km].real

    return rpnma, rpnms, pw


def build_ledir_blocks_from_local_legendre(
    nlat: int,
    km: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build one provisional OpenIFS-style ``RPNMA/RPNMS`` pair locally.

    This helper is a narrower, bridge-free first step toward the OpenIFS
    ``SULEG -> RPNMA/RPNMS`` architecture. It does **not** implement the full
    OpenIFS setup path; instead it reconstructs the same block layout directly
    from Gaussian latitudes, weights, and public stored Legendre values.

    Current scope
    -------------
    - one truncation defined by ``nlat - 1``
    - one zonal wavenumber ``km``
    - north-half Gaussian rows only
    - no butterfly / MPI / legpol file IO
    """
    from ._spherepack import getlegfunc
    from .spherical_harmonics import gaussian_lats_wts

    if nlat < 3:
        raise ValueError(f"nlat must be at least 3, got {nlat}")

    ntrunc = nlat - 1
    if km < 0 or km > ntrunc:
        raise ValueError(f"require 0 <= km <= {ntrunc}, got {km}")

    ndgnh = (nlat + 1) // 2
    ia = 1 + ((ntrunc - km + 2) % 2)
    is_ = 1 + ((ntrunc - km + 1) % 2)
    ila = (ntrunc - km + 2) // 2
    ils = (ntrunc - km + 3) // 2

    lats_deg, weights = gaussian_lats_wts(nlat)
    north_lats = np.asarray(lats_deg[:ndgnh], dtype=np.float64)
    pw = np.asarray(weights[:ndgnh], dtype=np.float64)

    ncoeff = ((ntrunc + 1) * (ntrunc + 2)) // 2
    legfunc_all = np.empty((ndgnh, ncoeff), dtype=np.float64)
    for ilat, lat in enumerate(north_lats):
        legfunc_all[ilat, :] = np.asarray(
            getlegfunc(float(lat), int(ntrunc)),
            dtype=np.float64,
        )

    rpnma = np.zeros((ndgnh, ila), dtype=np.float64)
    rpnms = np.zeros((ndgnh, ils), dtype=np.float64)

    start_idx = 0
    for m in range(ntrunc + 1):
        block = ntrunc - m + 1
        if m == km:
            for local_n in range(block):
                coeff_index = start_idx + local_n
                n = m + local_n
                row = ntrunc + 2 - n
                if (row % 2) == (ia % 2):
                    col = (row - ia) // 2
                    if 0 <= col < ila:
                        rpnma[:, col] = legfunc_all[:, coeff_index]
                else:
                    col = (row - is_) // 2
                    if 0 <= col < ils:
                        rpnms[:, col] = legfunc_all[:, coeff_index]
            break
        start_idx += block

    return rpnma, rpnms, pw


def _supolf_tables(
    ntrunc: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build the small OpenIFS `TPM_POL` scalar tables needed by `SUPOLF`."""
    jn = np.arange(ntrunc + 1, dtype=np.float64)
    dfa = np.zeros((ntrunc + 1,), dtype=np.float64)
    dfb = np.zeros((ntrunc + 1,), dtype=np.float64)
    dff = np.zeros((ntrunc + 1,), dtype=np.float64)
    dfg = np.zeros((ntrunc + 1,), dtype=np.float64)
    valid = jn >= 1.0
    jj = jn[valid]
    dfa[valid] = 1.0 / np.sqrt(jj * (jj + 1.0))
    dfb[valid] = np.sqrt((2.0 * jj + 1.0) / (jj * (jj + 1.0)))
    dff[valid] = (2.0 * jj - 1.0) / jj
    dfg[valid] = (jj - 1.0) / jj
    return dfa, dfb, dff, jn


def _supolf_scalar(mu: float, ntrunc: int, km: int, cheap: int) -> np.ndarray:
    """Python diagnostic clone of the narrow `SUPOLF` branch used by `SULEG`."""
    if cheap not in (1, 2, 3):
        raise ValueError(f"cheap must be 1, 2, or 3, got {cheap}")
    if km < 0 or km > ntrunc:
        raise ValueError(f"require 0 <= km <= {ntrunc}, got {km}")

    dfa, dfb, dff, dfg = _supolf_tables(ntrunc)
    ddpol = np.zeros((ntrunc + 1,), dtype=np.float64)

    x = float(mu)
    sin2 = max(0.0, 1.0 - x * x)
    sin1 = np.sqrt(sin2)
    eps = np.finfo(np.float64).eps
    inv_sin1 = 0.0 if sin1 <= eps else 1.0 / sin1

    dlkm2 = 1.0
    dlkm1 = x

    if km == 0:
        ddpol[0] = dlkm2
        if ntrunc >= 1:
            ddpol[1] = dlkm1 * dfb[1] / dfa[1]
        for jn in range(2, ntrunc + 1):
            dlk = dff[jn] * x * dlkm1 - dfg[jn] * dlkm2
            ddpol[jn] = dlk * dfb[jn] / dfa[jn]
            dlkm2 = dlkm1
            dlkm1 = dlk
        return ddpol

    if km == 1:
        ddpol[0] = 0.0
        if ntrunc >= 1:
            ddpol[1] = sin1 * dfb[1]
        for jn in range(2, ntrunc + 1):
            dlk = dff[jn] * x * dlkm1 - dfg[jn] * dlkm2
            dl1 = float(jn) * (dlkm1 - x * dlk) * inv_sin1
            ddpol[jn] = dl1 * dfb[jn]
            dlkm2 = dlkm1
            dlkm1 = dlk
        return ddpol

    zscale = 1.0e100
    ziscale = 1.0e-100
    iend = km // 2
    zlsita = 1.0
    icorr3 = 0
    for _ in range(1, iend + 1):
        zlsita *= sin2
        if abs(zlsita) < ziscale:
            zlsita *= zscale
            icorr3 += 1
    if km % 2 == 1:
        zlsita *= sin1

    zfac = 1.0
    for jn in range(1, km):
        zfac *= np.sqrt(float(2 * jn - 1))
        zfac /= np.sqrt(float(2 * jn))
    zfac *= np.sqrt(float(2 * km - 1))

    zfac0 = 1.0
    zfac1 = 1.0
    max_ic = min(ntrunc - km, 3)
    icorr = np.full((ntrunc + 1,), icorr3, dtype=np.int32)

    for ic in range(max_ic + 1):
        zfac0 *= float(2 * km + ic)
        if ic == 0:
            zmult = zfac
        elif ic == 1:
            zfac *= float(2 * km + ic)
            zmult = zfac * x
        elif ic == 2:
            zmult = 0.5 * zfac * ((2.0 * km + 3.0) * x * x - 1.0)
        else:
            zfac *= float(2 * km + ic)
            zmult = (x * zfac * ((2.0 * km + 5.0) * x * x - 3.0)) / 6.0
        idx = km + ic
        ddpol[idx] = zlsita * zmult * np.sqrt(2.0 * (float(idx) + 0.5) * zfac1 / zfac0)
        zfac1 *= float(ic + 1)

    def _dcl(k: int) -> float:
        return np.sqrt(
            (
                float(k - km + 1)
                * float(k - km + 2)
                * float(k + km + 1)
                * float(k + km + 2)
            )
            / (
                float(2 * k + 1)
                * float(2 * k + 3)
                * float(2 * k + 3)
                * float(2 * k + 5)
            )
        )

    def _ddl(k: int) -> float:
        return (2.0 * float(k) * float(k + 1) - 2.0 * float(km * km) - 1.0) / (
            float(2 * k - 1) * float(2 * k + 3)
        )

    if cheap == 2:
        start = km + 2
        seed_parity = 0
    elif cheap == 3:
        start = km + 3
        seed_parity = 1
    else:
        start = km + 2
        seed_parity = None

    if cheap in (2, 3):
        for jn in range(start, ntrunc - 1, 2):
            if abs(ddpol[jn - 2]) > zscale:
                ddpol[jn - 2] /= zscale
                ddpol[jn] /= zscale
                icorr[jn - 2 :] -= 1
            ddpol[jn + 2] = (
                (x * x - _ddl(jn)) * ddpol[jn] - _dcl(jn - 2) * ddpol[jn - 2]
            ) / _dcl(jn)
        parity_start = km if seed_parity == 0 else km + 1
        for jn in range(parity_start, ntrunc + 1, 2):
            for _ in range(int(icorr[jn])):
                ddpol[jn] /= zscale
                if ddpol[jn] < eps:
                    ddpol[jn] = eps
        return ddpol

    for jn in range(km + 2, ntrunc - 1):
        if abs(ddpol[jn - 2]) > zscale:
            ddpol[jn - 2] /= zscale
            ddpol[jn - 1] /= zscale
            ddpol[jn] /= zscale
            ddpol[jn + 1] /= zscale
            icorr[jn - 2 :] -= 1
        ddpol[jn + 2] = (
            (x * x - _ddl(jn)) * ddpol[jn] - _dcl(jn - 2) * ddpol[jn - 2]
        ) / _dcl(jn)
    for jn in range(km, ntrunc + 1):
        for _ in range(int(icorr[jn])):
            ddpol[jn] /= zscale
            if ddpol[jn] < eps:
                ddpol[jn] = eps
    return ddpol


def build_ledir_blocks_from_supolf(
    nlat: int,
    km: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build one local `RPNMA/RPNMS` pair using a Python `SUPOLF` diagnostic clone.

    This is intentionally still a diagnostic helper:

    - the OpenIFS row/column structure is reproduced closely enough to recover
      the expected `LEDIR` peak-row placement on the current small reference
      case
    - the remaining mismatch against the trusted Gaussian-bridge builder is
      still a near-uniform normalization factor, so this helper is not yet
      promoted into the direct pilot chain
    """
    from .spherical_harmonics import gaussian_lats_wts

    if nlat < 3:
        raise ValueError(f"nlat must be at least 3, got {nlat}")

    ntrunc = nlat - 1
    if km < 0 or km > ntrunc:
        raise ValueError(f"require 0 <= km <= {ntrunc}, got {km}")

    ndgnh = (nlat + 1) // 2
    ila = (ntrunc - km + 2) // 2
    ils = (ntrunc - km + 3) // 2

    lats_deg, weights = gaussian_lats_wts(nlat)
    mu_north = np.sin(np.deg2rad(np.asarray(lats_deg[:ndgnh], dtype=np.float64)))
    pw = np.asarray(weights[:ndgnh], dtype=np.float64)

    if ((ntrunc - km) % 2) == 0:
        inmax_a = ntrunc + 1
        inmax_s = ntrunc
    else:
        inmax_a = ntrunc
        inmax_s = ntrunc + 1

    rpnma = np.zeros((ndgnh, ila), dtype=np.float64)
    rpnms = np.zeros((ndgnh, ils), dtype=np.float64)

    for jgl, mu in enumerate(mu_north):
        anti = _supolf_scalar(float(mu), int(inmax_a), int(km), cheap=3)
        sym = _supolf_scalar(float(mu), int(inmax_s), int(km), cheap=2)
        for ji in range(ila):
            jn = km + 2 * ji + 1
            if jn <= ntrunc:
                rpnma[jgl, ila - ji - 1] = anti[jn]
        for ji in range(ils):
            jn = km + 2 * ji
            if jn <= ntrunc:
                rpnms[jgl, ils - ji - 1] = sym[jn]

    return rpnma, rpnms, pw


def _build_ledir_blocks_for_direct_pilot(
    nlat: int,
    km: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the currently preferred direct-pilot `RPNMA/RPNMS` blocks.

    This keeps the runtime on the more OpenIFS-like local `SUPOLF` path while
    applying the currently validated small-case normalization needed to close
    the existing `ldfou2 -> ledir -> poa1_to_vordiv` pilot chain against the
    trusted bridge reference.
    """
    rpnma, rpnms, pw = build_ledir_blocks_from_supolf(nlat, km)
    scale = np.sqrt(0.5)
    rpnma *= scale
    rpnms *= scale
    return rpnma, rpnms, pw


def validate_nloen(nloen: np.ndarray) -> tuple[int, int]:
    """Validate reduced-grid longitude counts through the native stub."""
    return _backend.validate_nloen(
        np.ascontiguousarray(nloen, dtype=np.int32),
    )


def create_setup(
    nloen: np.ndarray,
    mu: np.ndarray,
    weights: np.ndarray,
    rsphere: float,
):
    """Create a native reduced-grid setup handle."""
    return _backend.create_setup(
        np.ascontiguousarray(nloen, dtype=np.int32),
        np.ascontiguousarray(mu, dtype=np.float64),
        np.ascontiguousarray(weights, dtype=np.float64),
        float(rsphere),
    )


def describe_setup(handle) -> dict:
    """Describe a native reduced-grid setup handle."""
    return _backend.describe_setup(handle)


def destroy_setup(handle) -> None:
    """Destroy a native reduced-grid setup handle."""
    _backend.destroy_setup(handle)


def scalar_analysis_stub(
    handle,
    datagrid: np.ndarray,
    ntrunc: int,
) -> tuple[np.ndarray, int]:
    """Call the native scalar-analysis stub and return ectrans-packed spectrum plus status."""
    return _backend.scalar_analysis_stub(
        handle,
        np.ascontiguousarray(datagrid, dtype=np.float64),
        int(ntrunc),
    )


def scalar_synthesis_stub(
    handle,
    dataspec: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Call the native scalar-synthesis stub on ectrans-packed spectrum and return grid plus status."""
    return _backend.scalar_synthesis_stub(
        handle,
        np.ascontiguousarray(dataspec, dtype=np.float64),
    )


def scalar_fourier_stub(
    handle,
    datagrid: np.ndarray,
    mmax: int,
) -> tuple[np.ndarray, int]:
    """Call the native reduced-grid Fourier helper and return complex Fourier blocks plus status."""
    return _backend.scalar_fourier_stub(
        handle,
        np.ascontiguousarray(datagrid, dtype=np.float64),
        int(mmax),
    )


def scalar_block_solve_stub(
    handle,
    basis: np.ndarray,
    observed: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Call the native reduced-grid weighted complex block solver and return the solution plus status."""
    return _backend.scalar_block_solve_stub(
        handle,
        np.ascontiguousarray(basis, dtype=np.complex128),
        np.ascontiguousarray(observed, dtype=np.complex128),
    )


def weighted_block_solve_stub(
    weights: np.ndarray,
    basis: np.ndarray,
    observed: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Call the native general weighted complex block solver and return the solution plus status."""
    return _backend.weighted_block_solve_stub(
        np.ascontiguousarray(weights, dtype=np.float64),
        np.ascontiguousarray(basis, dtype=np.complex128),
        np.ascontiguousarray(observed, dtype=np.complex128),
    )


def vrtdiv_analysis_stub(
    handle,
    ugrid: np.ndarray,
    vgrid: np.ndarray,
    ntrunc: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Call the native vector-analysis stub and return vrt/div spectra plus status."""
    return _backend.vrtdiv_analysis_stub(
        handle,
        np.ascontiguousarray(ugrid, dtype=np.float64),
        np.ascontiguousarray(vgrid, dtype=np.float64),
        int(ntrunc),
    )


def uv_synthesis_stub(
    handle,
    vrtspec: np.ndarray,
    divspec: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Call the native wind-synthesis stub and return grid winds plus status."""
    return _backend.uv_synthesis_stub(
        handle,
        np.ascontiguousarray(vrtspec, dtype=np.complex128),
        np.ascontiguousarray(divspec, dtype=np.complex128),
    )


def gradient_synthesis_stub(
    handle,
    chispec: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Call the native gradient-synthesis stub and return grid gradients plus status."""
    return _backend.gradient_synthesis_stub(
        handle,
        np.ascontiguousarray(chispec, dtype=np.complex128),
    )


def vordiv_to_uv(
    vrtspec: np.ndarray,
    divspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Call the experimental ectrans vordiv->uv spectral helper."""
    return _backend.vordiv_to_uv(
        np.ascontiguousarray(vrtspec, dtype=np.complex128),
        np.ascontiguousarray(divspec, dtype=np.complex128),
        float(rsphere),
    )


def prfi1b_uv_block(
    km: int,
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, int]:
    """Build a native PRFI1B-style UV `POA1` block from public spectral `u/v`."""
    return _backend.prfi1b_uv_block(
        int(km),
        np.ascontiguousarray(uspec, dtype=np.complex128),
        np.ascontiguousarray(vspec, dtype=np.complex128),
        float(rsphere),
    )


def vd2uv_uv_block(
    km: int,
    vrtspec: np.ndarray,
    divspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, int]:
    """Build a native VD2UV-internal UV `POA1` block from public `vrt/div` spectra."""
    return _backend.vd2uv_uv_block(
        int(km),
        np.ascontiguousarray(vrtspec, dtype=np.complex128),
        np.ascontiguousarray(divspec, dtype=np.complex128),
        float(rsphere),
    )


def uv_to_vordiv(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Call the experimental ectrans uv-to-vordiv spectral helper."""
    return _backend.uv_to_vordiv(
        np.ascontiguousarray(uspec, dtype=np.complex128),
        np.ascontiguousarray(vspec, dtype=np.complex128),
        float(rsphere),
    )


def create_uv_to_vordiv_block_setup(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
):
    """Create a cached native blockwise ``uv -> vrt/div`` setup handle.

    The cached basis depends on spectral geometry and ``rsphere`` only.
    Trailing batch dimensions on the input sample are accepted and ignored for
    the setup itself.
    """
    return _backend.create_uv_to_vordiv_block_setup(
        np.ascontiguousarray(uspec, dtype=np.complex128),
        np.ascontiguousarray(vspec, dtype=np.complex128),
        float(rsphere),
    )


def destroy_uv_to_vordiv_block_setup(handle) -> None:
    """Destroy a cached native blockwise ``uv -> vrt/div`` setup handle."""
    _backend.destroy_uv_to_vordiv_block_setup(handle)


def uv_to_vordiv_block_native_with_setup(
    handle,
    uspec: np.ndarray,
    vspec: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``vrt/div`` with a prebuilt native blockwise setup handle."""
    vrt, div, ierror = _backend.uv_to_vordiv_block_native_with_setup(
        handle,
        np.ascontiguousarray(uspec, dtype=np.complex128),
        np.ascontiguousarray(vspec, dtype=np.complex128),
    )
    if ierror != 0:
        raise RuntimeError(
            f"uv_to_vordiv_block_native_with_setup backend failed with ierror={ierror}"
        )
    return vrt, div


def ldfou2_uv_scaling(
    ntrunc: int,
    km: int,
    paia: np.ndarray,
    rsphere: float,
    psia: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Apply the experimental ectrans LDFOU2 uv scaling helper."""
    return _backend.ldfou2_uv_scaling(
        int(ntrunc),
        int(km),
        np.ascontiguousarray(paia, dtype=np.float64),
        float(rsphere),
        np.ascontiguousarray(psia, dtype=np.float64),
    )


def ledir_dgemm(
    ntrunc: int,
    km: int,
    paia: np.ndarray,
    psia: np.ndarray,
    rpnma: np.ndarray,
    rpnms: np.ndarray,
    pw: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Apply the experimental DGEMM-only OpenIFS LEDIR kernel."""
    return _backend.ledir_dgemm(
        int(ntrunc),
        int(km),
        np.ascontiguousarray(paia, dtype=np.float64),
        np.ascontiguousarray(psia, dtype=np.float64),
        np.ascontiguousarray(rpnma, dtype=np.float64),
        np.ascontiguousarray(rpnms, dtype=np.float64),
        np.ascontiguousarray(pw, dtype=np.float64),
    )


def poa1_to_vordiv(
    ntrunc: int,
    km: int,
    poa1: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Apply the experimental OpenIFS UVTVD+UPDSPB stage to one POA1 block."""
    return _backend.poa1_to_vordiv(
        int(ntrunc),
        int(km),
        np.ascontiguousarray(poa1, dtype=np.float64),
        float(rsphere),
    )


def direct_vordiv_pilot_from_uv_blocks(
    nlon: int,
    nlat: int,
    km: int,
    paia: np.ndarray,
    rsphere: float,
    psia: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Run the current provisional single-``km`` direct vordiv pilot chain.

    This helper is intentionally narrow. The current block builder is local and
    more OpenIFS-like than the older Gaussian-bridge path, but it still uses a
    small normalization calibration that is only validated on the current
    focused reference cases.

    1. build calibrated local ``SUPOLF``-based ``RPNMA/RPNMS`` blocks,
    2. apply ``LDFOU2`` scaling,
    3. apply DGEMM-only ``LEDIR``,
    4. apply ``UVTVD + UPDSPB`` through ``poa1_to_vordiv``.

    Returns
    -------
    poa1, vrtspec, divspec, pw, ierror
        ``pw`` is returned so callers can inspect or reuse the Gaussian weights
        used by the provisional block constructor.
    """
    ntrunc = int(nlat) - 1
    _ = int(nlon)
    rpnma, rpnms, pw = _build_ledir_blocks_for_direct_pilot(
        int(nlat),
        int(km),
    )
    paia_scaled, psia_scaled, ierror = ldfou2_uv_scaling(
        ntrunc,
        int(km),
        paia,
        float(rsphere),
        psia,
    )
    if ierror != 0:
        return (
            np.zeros((0, 0), dtype=np.float64),
            np.zeros((0,), dtype=np.complex128),
            np.zeros((0,), dtype=np.complex128),
            pw,
            ierror,
        )

    poa1, ierror = ledir_dgemm(
        ntrunc,
        int(km),
        paia_scaled,
        psia_scaled,
        rpnma,
        rpnms,
        pw,
    )
    if ierror != 0:
        return (
            poa1,
            np.zeros((0,), dtype=np.complex128),
            np.zeros((0,), dtype=np.complex128),
            pw,
            ierror,
        )

    vrtspec, divspec, ierror = poa1_to_vordiv(
        ntrunc,
        int(km),
        poa1,
        float(rsphere),
    )
    return poa1, vrtspec, divspec, pw, ierror


def direct_vordiv_pilot_sweep_from_uv_blocks(
    nlon: int,
    nlat: int,
    uv_blocks_by_km: list[tuple[np.ndarray, np.ndarray]],
    rsphere: float,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray, np.ndarray, int]:
    """Run the current provisional direct vordiv pilot for ``km = 0..ntrunc``.

    Notes
    -----
    This helper currently targets the single-field pilot shape only, i.e.
    each ``(paia, psia)`` block pair must have leading dimension ``4``.
    That keeps the sweep aligned with the currently validated
    ``poa1_to_vordiv(...)`` wrapper behavior.
    """
    ntrunc = int(nlat) - 1
    ncoeff = ((ntrunc + 1) * (ntrunc + 2)) // 2

    if len(uv_blocks_by_km) != ntrunc + 1:
        raise ValueError(
            f"uv_blocks_by_km must have length {ntrunc + 1}, got {len(uv_blocks_by_km)}"
        )

    poa1_by_km: list[np.ndarray] = []
    vrtspec_total = np.zeros((ncoeff,), dtype=np.complex128)
    divspec_total = np.zeros((ncoeff,), dtype=np.complex128)
    pw_ref = np.zeros(((nlat + 1) // 2,), dtype=np.float64)

    for km, (paia, psia) in enumerate(uv_blocks_by_km):
        paia_arr = np.ascontiguousarray(paia, dtype=np.float64)
        psia_arr = np.ascontiguousarray(psia, dtype=np.float64)

        if paia_arr.ndim != 2 or psia_arr.ndim != 2:
            raise ValueError("each paia/psia block pair must contain rank-2 arrays")
        if paia_arr.shape != psia_arr.shape:
            raise ValueError(
                f"paia/psia block shapes must match for km={km}, got "
                f"{paia_arr.shape} and {psia_arr.shape}"
            )
        if paia_arr.shape[0] != 4:
            raise ValueError(
                "current sweep pilot expects single-field blocks with first "
                f"dimension 4, got {paia_arr.shape[0]} for km={km}"
            )

        poa1, vrtspec, divspec, pw, ierror = direct_vordiv_pilot_from_uv_blocks(
            nlon,
            nlat,
            km,
            paia_arr,
            rsphere,
            psia_arr,
        )
        if ierror != 0:
            return poa1_by_km, vrtspec_total, divspec_total, pw_ref, ierror

        if km == 0:
            pw_ref = pw
        elif not np.array_equal(pw, pw_ref):
            raise ValueError(
                "provisional Gaussian bridge returned inconsistent pw across km"
            )

        poa1_by_km.append(poa1)
        vrtspec_total += vrtspec
        divspec_total += divspec

    return poa1_by_km, vrtspec_total, divspec_total, pw_ref, 0


def uv_to_vordiv_linear_pilot(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``vrt/div`` spectra by explicitly inverting ``vordiv_to_uv``.

    This is a trusted but expensive small-problem reference helper. It builds
    the dense real operator mapping ``[vrt, div]`` to ``[u, v]`` using repeated
    calls to :func:`vordiv_to_uv`, then solves the corresponding linear system.
    """
    uspec_arr = np.ascontiguousarray(uspec, dtype=np.complex128)
    vspec_arr = np.ascontiguousarray(vspec, dtype=np.complex128)

    if uspec_arr.shape != vspec_arr.shape:
        raise ValueError(
            f"uspec and vspec must have the same shape, got {uspec_arr.shape} and {vspec_arr.shape}"
        )
    if uspec_arr.ndim != 1:
        raise ValueError("uv_to_vordiv_linear_pilot currently expects rank-1 spectra")

    ncoeff = int(uspec_arr.shape[0])
    _infer_ntrunc_from_ncoeff(ncoeff)
    operator = np.zeros((4 * ncoeff, 4 * ncoeff), dtype=np.float64)

    for idx in range(ncoeff):
        basis_vrt = np.zeros((ncoeff,), dtype=np.complex128)
        basis_div = np.zeros((ncoeff,), dtype=np.complex128)

        basis_vrt[idx] = 1.0 + 0.0j
        u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
        if ierr != 0:
            raise RuntimeError(
                f"vordiv_to_uv failed while building vrt-real basis at index {idx}"
            )
        operator[:ncoeff, idx] = u.real
        operator[ncoeff : 2 * ncoeff, idx] = u.imag
        operator[2 * ncoeff : 3 * ncoeff, idx] = v.real
        operator[3 * ncoeff :, idx] = v.imag

        basis_vrt[idx] = 0.0 + 1.0j
        u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
        if ierr != 0:
            raise RuntimeError(
                f"vordiv_to_uv failed while building vrt-imag basis at index {idx}"
            )
        operator[:ncoeff, ncoeff + idx] = u.real
        operator[ncoeff : 2 * ncoeff, ncoeff + idx] = u.imag
        operator[2 * ncoeff : 3 * ncoeff, ncoeff + idx] = v.real
        operator[3 * ncoeff :, ncoeff + idx] = v.imag

        basis_vrt[idx] = 0.0 + 0.0j
        basis_div[idx] = 1.0 + 0.0j
        u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
        if ierr != 0:
            raise RuntimeError(
                f"vordiv_to_uv failed while building div-real basis at index {idx}"
            )
        operator[:ncoeff, 2 * ncoeff + idx] = u.real
        operator[ncoeff : 2 * ncoeff, 2 * ncoeff + idx] = u.imag
        operator[2 * ncoeff : 3 * ncoeff, 2 * ncoeff + idx] = v.real
        operator[3 * ncoeff :, 2 * ncoeff + idx] = v.imag

        basis_div[idx] = 0.0 + 1.0j
        u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
        if ierr != 0:
            raise RuntimeError(
                f"vordiv_to_uv failed while building div-imag basis at index {idx}"
            )
        operator[:ncoeff, 3 * ncoeff + idx] = u.real
        operator[ncoeff : 2 * ncoeff, 3 * ncoeff + idx] = u.imag
        operator[2 * ncoeff : 3 * ncoeff, 3 * ncoeff + idx] = v.real
        operator[3 * ncoeff :, 3 * ncoeff + idx] = v.imag

    rhs = np.concatenate(
        [
            uspec_arr.real,
            uspec_arr.imag,
            vspec_arr.real,
            vspec_arr.imag,
        ]
    )
    col_norms = np.linalg.norm(operator, axis=0)
    active = col_norms > 1e-18
    reduced_operator = operator[:, active]
    reduced_solution, _, _, _ = np.linalg.lstsq(reduced_operator, rhs, rcond=None)
    solution = np.zeros((4 * ncoeff,), dtype=np.float64)
    solution[active] = reduced_solution
    vrt = solution[:ncoeff] + 1j * solution[ncoeff : 2 * ncoeff]
    div = solution[2 * ncoeff : 3 * ncoeff] + 1j * solution[3 * ncoeff :]
    return vrt, div


def uv_to_vordiv_block_linear_pilot(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``vrt/div`` by solving one independent zonal-wavenumber block at a time.

    This helper keeps the same trusted forward operator as
    :func:`uv_to_vordiv_linear_pilot`, but respects the native ectrans structure
    more closely by solving one public ``m`` block at a time instead of building
    one full dense operator for the whole triangular spectrum.
    """
    uspec_arr = np.ascontiguousarray(uspec, dtype=np.complex128)
    vspec_arr = np.ascontiguousarray(vspec, dtype=np.complex128)

    if uspec_arr.shape != vspec_arr.shape:
        raise ValueError(
            f"uspec and vspec must have the same shape, got {uspec_arr.shape} and {vspec_arr.shape}"
        )
    if uspec_arr.ndim != 1:
        raise ValueError(
            "uv_to_vordiv_block_linear_pilot currently expects rank-1 spectra"
        )

    ncoeff = int(uspec_arr.shape[0])
    ntrunc = _infer_ntrunc_from_ncoeff(ncoeff)
    vrt = np.zeros((ncoeff,), dtype=np.complex128)
    div = np.zeros((ncoeff,), dtype=np.complex128)

    for m in range(ntrunc + 1):
        block_indices = _public_indices_for_zonal_wavenumber(ntrunc, m)
        block_size = int(block_indices.size)
        operator = np.zeros((4 * block_size, 4 * block_size), dtype=np.float64)

        for local_idx, global_idx in enumerate(block_indices):
            basis_vrt = np.zeros((ncoeff,), dtype=np.complex128)
            basis_div = np.zeros((ncoeff,), dtype=np.complex128)

            basis_vrt[global_idx] = 1.0 + 0.0j
            u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
            if ierr != 0:
                raise RuntimeError(
                    f"vordiv_to_uv failed while building m={m} vrt-real basis "
                    f"at global index {global_idx}"
                )
            operator[:block_size, local_idx] = u[block_indices].real
            operator[block_size : 2 * block_size, local_idx] = u[block_indices].imag
            operator[2 * block_size : 3 * block_size, local_idx] = v[block_indices].real
            operator[3 * block_size :, local_idx] = v[block_indices].imag

            basis_vrt[global_idx] = 0.0 + 1.0j
            u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
            if ierr != 0:
                raise RuntimeError(
                    f"vordiv_to_uv failed while building m={m} vrt-imag basis "
                    f"at global index {global_idx}"
                )
            operator[:block_size, block_size + local_idx] = u[block_indices].real
            operator[block_size : 2 * block_size, block_size + local_idx] = u[
                block_indices
            ].imag
            operator[2 * block_size : 3 * block_size, block_size + local_idx] = v[
                block_indices
            ].real
            operator[3 * block_size :, block_size + local_idx] = v[block_indices].imag

            basis_vrt[global_idx] = 0.0 + 0.0j
            basis_div[global_idx] = 1.0 + 0.0j
            u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
            if ierr != 0:
                raise RuntimeError(
                    f"vordiv_to_uv failed while building m={m} div-real basis "
                    f"at global index {global_idx}"
                )
            operator[:block_size, 2 * block_size + local_idx] = u[block_indices].real
            operator[block_size : 2 * block_size, 2 * block_size + local_idx] = u[
                block_indices
            ].imag
            operator[2 * block_size : 3 * block_size, 2 * block_size + local_idx] = v[
                block_indices
            ].real
            operator[3 * block_size :, 2 * block_size + local_idx] = v[
                block_indices
            ].imag

            basis_div[global_idx] = 0.0 + 1.0j
            u, v, ierr = vordiv_to_uv(basis_vrt, basis_div, rsphere)
            if ierr != 0:
                raise RuntimeError(
                    f"vordiv_to_uv failed while building m={m} div-imag basis "
                    f"at global index {global_idx}"
                )
            operator[:block_size, 3 * block_size + local_idx] = u[block_indices].real
            operator[block_size : 2 * block_size, 3 * block_size + local_idx] = u[
                block_indices
            ].imag
            operator[2 * block_size : 3 * block_size, 3 * block_size + local_idx] = v[
                block_indices
            ].real
            operator[3 * block_size :, 3 * block_size + local_idx] = v[
                block_indices
            ].imag

        rhs = np.concatenate(
            [
                uspec_arr[block_indices].real,
                uspec_arr[block_indices].imag,
                vspec_arr[block_indices].real,
                vspec_arr[block_indices].imag,
            ]
        )
        col_norms = np.linalg.norm(operator, axis=0)
        active = col_norms > 1e-18
        reduced_operator = operator[:, active]
        reduced_solution, _, _, _ = np.linalg.lstsq(reduced_operator, rhs, rcond=None)
        solution = np.zeros((4 * block_size,), dtype=np.float64)
        solution[active] = reduced_solution

        vrt[block_indices] = (
            solution[:block_size] + 1j * solution[block_size : 2 * block_size]
        )
        div[block_indices] = solution[2 * block_size : 3 * block_size] + 1j * (
            solution[3 * block_size :]
        )

    return vrt, div


def uv_to_vordiv_poa1_block_linear_pilot(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``vrt/div`` by inverting the current PRFI/POA1 boundary per ``m`` block.

    This pilot is narrower than :func:`uv_to_vordiv_block_linear_pilot`: for
    each zonal wavenumber it builds the linear map from a forward direct-chain
    ``POA1`` uv block to the currently observed native PRFI-style uv block,
    solves that small least-squares problem, and then reuses
    :func:`poa1_to_vordiv` to return spectra.

    It is therefore closer to the active ectrans direct-chain semantics than
    the pure spectral block-linear reference, while still staying in Python.
    """
    with UvToVordivPoa1BlockSolver(uspec, vspec, rsphere) as solver:
        return solver.solve(uspec, vspec)


def uv_to_vordiv_poa1_block_cached_pilot(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the narrowed PRFI/POA1 pilot through the module-level cache."""
    cache = _get_uv_to_vordiv_poa1_block_pilot_cache()
    return cache.solve(uspec, vspec, float(rsphere))


def uv_to_vordiv_block_native_pilot(
    uspec: np.ndarray,
    vspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``vrt/div`` with the native weighted block solver per zonal block.

    This stays on the same block structure as
    :func:`uv_to_vordiv_block_linear_pilot`, but replaces the NumPy least-squares
    solve with the existing native complex block solver already used elsewhere in
    the experimental reduced-grid backend.
    """
    uspec_arr = np.ascontiguousarray(uspec, dtype=np.complex128)
    vspec_arr = np.ascontiguousarray(vspec, dtype=np.complex128)
    if uspec_arr.shape != vspec_arr.shape:
        raise ValueError(
            f"uspec and vspec must have the same shape, got {uspec_arr.shape} and {vspec_arr.shape}"
        )
    cache = _get_uv_to_vordiv_block_native_pilot_cache()
    return cache.solve(uspec_arr, vspec_arr, float(rsphere))
