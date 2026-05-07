"""Thin Python wrapper around the experimental native ectrans backend."""

from __future__ import annotations

import importlib

import numpy as np

_backend = importlib.import_module("skyborn.spharm._ectrans_backend")


def _public_spectral_index(ntrunc: int, m: int, n: int) -> int:
    """Return the public triangular coefficient index for one ``(m, n)``."""
    idx = 0
    for mm in range(ntrunc + 1):
        for nn in range(mm, ntrunc + 1):
            if mm == m and nn == n:
                return idx
            idx += 1
    raise ValueError(f"invalid (m, n)=({m}, {n}) for ntrunc={ntrunc}")


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

    This helper is intentionally narrow and stays on the existing provisional
    bridge:

    1. build provisional ``RPNMA/RPNMS`` blocks from Gaussian synthesis+FFT,
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
    rpnma, rpnms, pw = build_ledir_blocks_from_scalar_bridge(
        int(nlon),
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
