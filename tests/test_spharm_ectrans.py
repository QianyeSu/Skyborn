from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

import skyborn.spharm as spharm_pkg  # noqa: E402
from skyborn.spharm import _spherepack  # noqa: E402
from skyborn.spharm.ectrans_backend_api import (  # noqa: E402
    UvToVordivBlockSolver,
    build_ledir_blocks_from_scalar_bridge,
    create_uv_to_vordiv_block_setup,
    destroy_uv_to_vordiv_block_setup,
    direct_vordiv_pilot_from_uv_blocks,
    direct_vordiv_pilot_sweep_from_uv_blocks,
    ldfou2_uv_scaling,
    ledir_dgemm,
    poa1_to_vordiv,
    scalar_analysis_stub,
    scalar_block_solve_stub,
    scalar_fourier_stub,
    scalar_synthesis_stub,
    uv_to_vordiv,
    uv_to_vordiv_block_linear_pilot,
    uv_to_vordiv_block_native_pilot,
    uv_to_vordiv_block_native_with_setup,
    uv_to_vordiv_linear_pilot,
    vordiv_to_uv,
    weighted_block_solve_stub,
)
from skyborn.spharm.reduced_gaussian import (
    _ECTRANS_BACKEND_AVAILABLE,
    ReducedGaussianGrid,
    ReducedGaussianSpharmt,
    ReducedGaussianValidationError,
)
from skyborn.spharm.reduced_gaussian import regrid as reduced_regrid  # noqa: E402
from skyborn.spharm.reduced_gaussian import regriduv as reduced_regriduv


def _single_set_ectrans_offsets(ntrunc: int) -> tuple[np.ndarray, np.ndarray]:
    """Return OpenIFS-style single-set NASM0/NPMT offsets for one local wave set."""
    nasm0 = np.empty(ntrunc + 1, dtype=np.int64)
    npm_t = np.empty(ntrunc + 1, dtype=np.int64)

    spec_offset = 0
    leg_offset = 0
    for m in range(ntrunc + 1):
        nasm0[m] = spec_offset
        npm_t[m] = leg_offset
        spec_offset += 2 * (ntrunc - m + 1)
        leg_offset += ntrunc + 2 - m
    return nasm0, npm_t


def _python_vordiv_to_uv_reference(
    vrtspec: np.ndarray,
    divspec: np.ndarray,
    rsphere: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Pure-Python transcription of the narrowed OpenIFS vordiv->uv pilot chain."""
    vrtspec = np.asarray(vrtspec, dtype=np.complex128)
    divspec = np.asarray(divspec, dtype=np.complex128)
    if vrtspec.shape != divspec.shape:
        raise ValueError("vrtspec and divspec must have the same shape")

    ncoeff = vrtspec.shape[0]
    ntrunc = int(-1.5 + 0.5 * np.sqrt(1.0 + 8.0 * ncoeff))
    if ((ntrunc + 1) * (ntrunc + 2)) // 2 != ncoeff:
        raise ValueError("invalid triangular coefficient count")

    extra_shape = vrtspec.shape[1:]
    nt = int(np.prod(extra_shape, dtype=int)) if extra_shape else 1
    vrt_flat = vrtspec.reshape(ncoeff, nt)
    div_flat = divspec.reshape(ncoeff, nt)

    nasm0, npm_t = _single_set_ectrans_offsets(ntrunc)
    nspec2 = int(2 * ncoeff)
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2)

    pspvor = np.zeros((nt, nspec2), dtype=np.float64)
    pspdiv = np.zeros((nt, nspec2), dtype=np.float64)
    pspu = np.zeros((nt, nspec2), dtype=np.float64)
    pspv = np.zeros((nt, nspec2), dtype=np.float64)

    idx = 0
    for m in range(ntrunc + 1):
        for n in range(m, ntrunc + 1):
            inm = int(nasm0[m] + 2 * (n - m))
            pspvor[:, inm] = vrt_flat[idx].real
            pspvor[:, inm + 1] = vrt_flat[idx].imag
            pspdiv[:, inm] = div_flat[idx].real
            pspdiv[:, inm + 1] = div_flat[idx].imag
            idx += 1

    rn = np.arange(-1, ntrunc + 4, dtype=np.float64)
    rlapin = np.zeros(ntrunc + 4, dtype=np.float64)
    for n in range(1, ntrunc + 3):
        rlapin[n + 1] = -(rsphere * rsphere) / float(n * (n + 1))

    repsnm = []
    for m in range(ntrunc + 1):
        for n in range(m, ntrunc + 3):
            repsnm.append(np.sqrt((n * n - m * m) / float(4 * n * n - 1)))
    repsnm = np.asarray(repsnm, dtype=np.float64)

    for m in range(ntrunc + 1):
        zei = np.zeros(ntrunc + 3, dtype=np.float64)
        if m > 0:
            zei[:m] = 0.0
        km_loc = m + 1
        for n in range(m, ntrunc + 3):
            zei[n] = repsnm[int(npm_t[m] + n)]

        zia = np.zeros((nlei1, 8 * nt), dtype=np.float64)
        ilcm = ntrunc + 1 - m
        ioff = int(nasm0[m])

        for j in range(1, ilcm + 1):
            inm = ioff + (ilcm - j) * 2
            row = j + 1
            for fld in range(nt):
                ir = 2 * fld
                ii = ir + 1
                zia[row, ir] = pspvor[fld, inm]
                zia[row, ii] = pspvor[fld, inm + 1]
                zia[row, 2 * nt + ir] = pspdiv[fld, inm]
                zia[row, 2 * nt + ii] = pspdiv[fld, inm + 1]

        zn = np.zeros(ntrunc + 6, dtype=np.float64)
        zlap = np.zeros(ntrunc + 6, dtype=np.float64)
        zeps = np.zeros(ntrunc + 6, dtype=np.float64)
        ismax = ntrunc
        for n in range(m - 1, ismax + 3):
            ij = ismax + 3 - n
            zn[ij + 1] = rn[n + 1]
            zlap[ij + 1] = rlapin[n + 1]
            if n >= 0:
                zeps[ij + 1] = zei[n]
        zn[1] = rn[ismax + 4]

        if m == 0:
            for fld in range(nt):
                ir = 2 * fld
                for ji in range(2, ismax + 4 - m):
                    zia[ji - 1, 4 * nt + ir] = (
                        zn[ji + 2] * zeps[ji + 1] * zlap[ji + 2] * zia[ji, ir]
                        - zn[ji - 1] * zeps[ji] * zlap[ji] * zia[ji - 2, ir]
                    )
                    zia[ji - 1, 6 * nt + ir] = (
                        -zn[ji + 2] * zeps[ji + 1] * zlap[ji + 2] * zia[ji, 2 * nt + ir]
                        + zn[ji - 1] * zeps[ji] * zlap[ji] * zia[ji - 2, 2 * nt + ir]
                    )
        else:
            zkm = float(m)
            for fld in range(nt):
                ir = 2 * fld
                ii = ir + 1
                for ji in range(2, ismax + 4 - m):
                    zia[ji - 1, 4 * nt + ir] = (
                        -zkm * zlap[ji + 1] * zia[ji - 1, 2 * nt + ii]
                        + zn[ji + 2] * zeps[ji + 1] * zlap[ji + 2] * zia[ji, ir]
                        - zn[ji - 1] * zeps[ji] * zlap[ji] * zia[ji - 2, ir]
                    )
                    zia[ji - 1, 4 * nt + ii] = (
                        zkm * zlap[ji + 1] * zia[ji - 1, 2 * nt + ir]
                        + zn[ji + 2] * zeps[ji + 1] * zlap[ji + 2] * zia[ji, ii]
                        - zn[ji - 1] * zeps[ji] * zlap[ji] * zia[ji - 2, ii]
                    )
                    zia[ji - 1, 6 * nt + ir] = (
                        -zkm * zlap[ji + 1] * zia[ji - 1, ii]
                        - zn[ji + 2]
                        * zeps[ji + 1]
                        * zlap[ji + 2]
                        * zia[ji, 2 * nt + ir]
                        + zn[ji - 1] * zeps[ji] * zlap[ji] * zia[ji - 2, 2 * nt + ir]
                    )
                    zia[ji - 1, 6 * nt + ii] = (
                        zkm * zlap[ji + 1] * zia[ji - 1, ir]
                        - zn[ji + 2]
                        * zeps[ji + 1]
                        * zlap[ji + 2]
                        * zia[ji, 2 * nt + ii]
                        + zn[ji - 1] * zeps[ji] * zlap[ji] * zia[ji - 2, 2 * nt + ii]
                    )

        inv_ra = 1.0 / rsphere
        for j in range(1, ilcm + 1):
            inm = ioff + (ilcm - j) * 2
            row = j + 1
            for fld in range(nt):
                ir = 2 * fld
                ii = ir + 1
                pspu[fld, inm] = zia[row, 4 * nt + ir] * inv_ra
                pspu[fld, inm + 1] = zia[row, 4 * nt + ii] * inv_ra
                pspv[fld, inm] = zia[row, 6 * nt + ir] * inv_ra
                pspv[fld, inm + 1] = zia[row, 6 * nt + ii] * inv_ra

    uspec = np.empty((ncoeff, nt), dtype=np.complex128)
    vspec = np.empty((ncoeff, nt), dtype=np.complex128)
    idx = 0
    for m in range(ntrunc + 1):
        for n in range(m, ntrunc + 1):
            inm = int(nasm0[m] + 2 * (n - m))
            uspec[idx, :] = pspu[:, inm] + 1j * pspu[:, inm + 1]
            vspec[idx, :] = pspv[:, inm] + 1j * pspv[:, inm + 1]
            idx += 1

    if extra_shape:
        return uspec.reshape((ncoeff,) + extra_shape), vspec.reshape(
            (ncoeff,) + extra_shape
        )
    return uspec[:, 0], vspec[:, 0]


def _spectral_index(ntrunc: int, m: int, n: int) -> int:
    """Return the public triangular coefficient index for one (m, n)."""
    idx = 0
    for mm in range(ntrunc + 1):
        for nn in range(mm, ntrunc + 1):
            if mm == m and nn == n:
                return idx
            idx += 1
    raise ValueError(f"invalid (m, n)=({m}, {n}) for ntrunc={ntrunc}")


def test_reduced_gaussian_grid_basic_metadata():
    grid = ReducedGaussianGrid(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    assert grid.ndgl == 5
    assert grid.ngptot == 116
    assert grid.max_nlon == 28


def test_reduced_gaussian_public_exports_are_available_from_spharm_package():
    assert spharm_pkg.ReducedGaussianGrid is ReducedGaussianGrid
    assert spharm_pkg.ReducedGaussianSpharmt is ReducedGaussianSpharmt
    assert spharm_pkg.regrid is not None
    assert spharm_pkg.regriduv is not None


def test_reduced_gaussian_spharmt_repr_contains_key_sizes():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    text = repr(backend)
    assert "ndgl=5" in text
    assert "ngptot=116" in text


def test_reduced_gaussian_spharmt_scalar_prototype_returns_result():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    field = np.zeros((116,), dtype=np.float64)
    spec = backend.grdtospec(field)
    assert spec.shape == (15,)
    np.testing.assert_allclose(spec, 0.0, rtol=0.0, atol=1e-12)
    backend.close()


def test_reduced_gaussian_vector_prototype_matches_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend._native_uv_synthesis_available = False
    backend._native_vrtdiv_analysis_available = False
    backend._native_gradient_synthesis_available = False
    bridge = backend._get_bridge_spharmt()

    spec = np.zeros((6,), dtype=np.complex128)
    spec[1] = 0.25 - 0.1j
    spec[2] = -0.2 + 0.05j
    spec[4] = 0.15 + 0.12j

    full_scalar = bridge.spectogrd(spec)
    packed_scalar = backend._bridge_to_packed_grid(full_scalar, "test packed scalar")
    rebuilt_full_scalar = backend._packed_to_bridge_grid(
        backend._validate_packed_grid_data(packed_scalar, "test packed scalar")[1],
        (),
    )
    np.testing.assert_allclose(rebuilt_full_scalar, full_scalar, rtol=0.0, atol=5e-8)

    u_expected, v_expected = bridge.getgrad(spec)
    packed_u = backend._bridge_to_packed_grid(u_expected, "test packed u")
    packed_v = backend._bridge_to_packed_grid(v_expected, "test packed v")

    vrt_expected, div_expected = bridge.getvrtdivspec(u_expected, v_expected, ntrunc=2)
    vrt_actual, div_actual = backend.getvrtdivspec(packed_u, packed_v, ntrunc=2)
    np.testing.assert_allclose(vrt_actual, vrt_expected, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(div_actual, div_expected, rtol=0.0, atol=1e-10)

    u_actual, v_actual = backend.getuv(vrt_expected, div_expected)
    np.testing.assert_allclose(u_actual, packed_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(v_actual, packed_v, rtol=0.0, atol=1e-10)

    grad_u_actual, grad_v_actual = backend.getgrad(spec)
    np.testing.assert_allclose(grad_u_actual, packed_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(grad_v_actual, packed_v, rtol=0.0, atol=1e-10)

    backend.close()


def test_reduced_gaussian_default_native_vector_capabilities_follow_backend():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))

    assert backend._supports_native_uv_synthesis() is _ECTRANS_BACKEND_AVAILABLE
    assert backend._supports_native_vrtdiv_analysis() is _ECTRANS_BACKEND_AVAILABLE
    assert backend._supports_native_gradient_synthesis() is _ECTRANS_BACKEND_AVAILABLE
    backend.close()


def test_reduced_gaussian_native_vector_pipeline_recovers_full_grid_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()

    chi_spec = np.zeros((6,), dtype=np.complex128)
    chi_spec[1] = 0.25 - 0.1j
    chi_spec[2] = -0.2 + 0.05j
    chi_spec[4] = 0.15 + 0.12j

    psi_spec = np.zeros((6,), dtype=np.complex128)
    psi_spec[1] = 0.21 - 0.04j
    psi_spec[3] = -0.13 + 0.07j
    psi_spec[4] = 0.09 - 0.02j

    full_u_div, full_v_div = bridge.getgrad(chi_spec)
    packed_u_div = backend._bridge_to_packed_grid(full_u_div, "native div u")
    packed_v_div = backend._bridge_to_packed_grid(full_v_div, "native div v")
    psi_ref_div = bridge.getpsispec(full_u_div, full_v_div, ntrunc=2)
    chi_ref_div = bridge.getchispec(full_u_div, full_v_div, ntrunc=2)

    vrt_spec = _spherepack.lap(psi_spec.reshape(6, 1), backend.rsphere).reshape(
        6,
    )
    full_u_vrt, full_v_vrt = bridge.getuv(vrt_spec, np.zeros_like(vrt_spec))
    packed_u_vrt = backend._bridge_to_packed_grid(full_u_vrt, "native vrt u")
    packed_v_vrt = backend._bridge_to_packed_grid(full_v_vrt, "native vrt v")
    psi_ref_vrt = bridge.getpsispec(full_u_vrt, full_v_vrt, ntrunc=2)
    chi_ref_vrt = bridge.getchispec(full_u_vrt, full_v_vrt, ntrunc=2)

    backend._native_uv_synthesis_available = True
    backend._native_vrtdiv_analysis_available = True
    backend._native_gradient_synthesis_available = True

    backend._vector_analysis_basis_cache.clear()
    psi_div = backend.getpsispec(packed_u_div, packed_v_div, ntrunc=2)
    chi_div = backend.getchispec(packed_u_div, packed_v_div, ntrunc=2)
    np.testing.assert_allclose(psi_div, psi_ref_div, rtol=0.0, atol=2e-8)
    np.testing.assert_allclose(chi_div, chi_ref_div, rtol=0.0, atol=3e-8)

    backend._vector_analysis_basis_cache.clear()
    psi_vrt = backend.getpsispec(packed_u_vrt, packed_v_vrt, ntrunc=2)
    chi_vrt = backend.getchispec(packed_u_vrt, packed_v_vrt, ntrunc=2)
    np.testing.assert_allclose(psi_vrt, psi_ref_vrt, rtol=0.0, atol=3e-8)
    np.testing.assert_allclose(chi_vrt, chi_ref_vrt, rtol=0.0, atol=1e-8)

    backend.close()


def test_vordiv_to_uv_matches_python_reference():
    sht = spharm_pkg.Spharmt(nlon=32, nlat=17, gridtype="gaussian", legfunc="stored")
    rng = np.random.default_rng(20260507)
    ncoeff = (sht.nlat * (sht.nlat + 1)) // 2
    vrtspec = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    divspec = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    uspec_ref, vspec_ref = _python_vordiv_to_uv_reference(vrtspec, divspec, sht.rsphere)

    uspec_pilot, vspec_pilot, ierror = vordiv_to_uv(vrtspec, divspec, sht.rsphere)

    assert ierror == 0
    np.testing.assert_allclose(uspec_pilot, uspec_ref, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(vspec_pilot, vspec_ref, rtol=0.0, atol=1e-10)


def test_vordiv_to_uv_supports_trailing_dimensions():
    sht = spharm_pkg.Spharmt(nlon=24, nlat=13, gridtype="gaussian", legfunc="stored")
    rng = np.random.default_rng(20260508)
    ncoeff = (sht.nlat * (sht.nlat + 1)) // 2
    vrtspec = rng.standard_normal((ncoeff, 2, 2)).astype(
        np.float64
    ) + 1j * rng.standard_normal((ncoeff, 2, 2)).astype(np.float64)
    divspec = rng.standard_normal((ncoeff, 2, 2)).astype(
        np.float64
    ) + 1j * rng.standard_normal((ncoeff, 2, 2)).astype(np.float64)

    uspec_pilot, vspec_pilot, ierror = vordiv_to_uv(vrtspec, divspec, sht.rsphere)

    assert ierror == 0
    assert uspec_pilot.shape == vrtspec.shape
    assert vspec_pilot.shape == divspec.shape


def test_vordiv_to_uv_rejects_mismatched_shapes():
    vrtspec = np.zeros((10,), dtype=np.complex128)
    divspec = np.zeros((10, 2), dtype=np.complex128)

    try:
        vordiv_to_uv(vrtspec, divspec, 6.3712e6)
    except Exception as exc:
        assert "same rank" in str(exc) or "same shape" in str(exc)
    else:
        raise AssertionError("Expected shape validation error for vordiv_to_uv")


def test_ldfou2_uv_scaling_matches_analytic_gaussian_factor():
    ntrunc = 4
    km = 0
    kf_uv = 1
    rsphere = 6.3712e6
    nlat = ntrunc + 1
    ndgnh = (nlat + 1) // 2

    paia = np.ones((4 * kf_uv, ndgnh), dtype=np.float64)
    psia = np.ones((4 * kf_uv, ndgnh), dtype=np.float64)
    paia_out, psia_out, ierror = ldfou2_uv_scaling(ntrunc, km, paia, rsphere, psia)

    assert ierror == 0
    lats_deg, _ = spharm_pkg.gaussian_lats_wts(nlat)
    north_lats = lats_deg[:ndgnh]
    expected = 1.0 / (rsphere * np.cos(np.deg2rad(north_lats)))

    np.testing.assert_allclose(paia_out[0], expected, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(psia_out[0], expected, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(
        paia_out, np.broadcast_to(expected, paia_out.shape), rtol=0.0, atol=2e-12
    )
    np.testing.assert_allclose(
        psia_out, np.broadcast_to(expected, psia_out.shape), rtol=0.0, atol=2e-12
    )


def test_ledir_dgemm_matches_manual_small_matrix_case():
    ntrunc = 4
    km = 1
    kfc = 2
    kdglu = 3
    ila = (ntrunc - km + 2) // 2
    ils = (ntrunc - km + 3) // 2
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2)

    paia = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
        ],
        dtype=np.float64,
    )
    psia = np.array(
        [
            [0.5, 1.0, 1.5],
            [2.0, 2.5, 3.0],
        ],
        dtype=np.float64,
    )
    rpnma = np.array(
        [
            [1.0, 0.1],
            [0.2, 1.5],
            [0.3, 0.4],
        ],
        dtype=np.float64,
    )
    rpnms = np.array(
        [
            [0.6, 0.0, 0.2],
            [0.1, 0.7, 0.4],
            [0.2, 0.3, 0.5],
        ],
        dtype=np.float64,
    )
    pw = np.array([1.0, 2.0, 0.5], dtype=np.float64)

    actual, ierror = ledir_dgemm(ntrunc, km, paia, psia, rpnma, rpnms, pw)

    assert ierror == 0
    assert actual.shape == (nlei1, kfc)

    expected = np.zeros((nlei1, kfc), dtype=np.float64)
    ia = 1 + ((ntrunc - km + 2) % 2)
    is_ = 1 + ((ntrunc - km + 1) % 2)
    for jk in range(kfc):
        weighted_a = paia[jk] * pw
        weighted_s = psia[jk] * pw
        zca = rpnma.T @ weighted_a
        zcs = rpnms.T @ weighted_s
        for j in range(ila):
            expected[ia - 1 + 2 * j, jk] = zca[j]
        for j in range(ils):
            expected[is_ - 1 + 2 * j, jk] = zcs[j]

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-12)


def test_ledir_dgemm_reproduces_gaussian_bridge_scalar_mode_rows():
    nlat = 17
    nlon = 32
    ntrunc = nlat - 1
    km = 2
    ndgnh = (nlat + 1) // 2
    ila = (ntrunc - km + 2) // 2
    ils = (ntrunc - km + 3) // 2
    ia = 1 + ((ntrunc - km + 2) % 2)
    is_ = 1 + ((ntrunc - km + 1) % 2)
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2)

    rpnma, rpnms, pw = build_ledir_blocks_from_scalar_bridge(nlon, nlat, km)

    for n in range(km, ntrunc + 1):
        row = ntrunc + 2 - n
        if row % 2 == ia % 2:
            paia = np.zeros((1, ndgnh), dtype=np.float64)
            psia = np.zeros((1, ndgnh), dtype=np.float64)
            col = (row - ia) // 2
            paia[0, :] = rpnma[:, col] / pw
            actual, ierror = ledir_dgemm(ntrunc, km, paia, psia, rpnma, rpnms, pw)
        else:
            paia = np.zeros((1, ndgnh), dtype=np.float64)
            psia = np.zeros((1, ndgnh), dtype=np.float64)
            col = (row - is_) // 2
            psia[0, :] = rpnms[:, col] / pw
            actual, ierror = ledir_dgemm(ntrunc, km, paia, psia, rpnma, rpnms, pw)

        assert ierror == 0
        assert actual.shape == (nlei1, 1)
        peak_row = int(np.argmax(np.abs(actual[:, 0])))
        assert peak_row == row - 1


def test_build_ledir_blocks_from_scalar_bridge_returns_expected_shapes():
    nlat = 17
    nlon = 32
    km = 2

    rpnma, rpnms, pw = build_ledir_blocks_from_scalar_bridge(nlon, nlat, km)

    assert rpnma.shape == (9, 8)
    assert rpnms.shape == (9, 8)
    assert pw.shape == (9,)
    assert np.all(np.isfinite(rpnma))
    assert np.all(np.isfinite(rpnms))
    assert np.all(pw > 0.0)


def test_poa1_to_vordiv_places_single_mode_on_expected_spectral_index():
    ntrunc = 4
    km = 2
    rsphere = 6.3712e6
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2)
    poa1 = np.zeros((nlei1, 4), dtype=np.float64)

    ia = 1 + ((ntrunc - km + 2) % 2)
    target_row = ia
    poa1[target_row - 1, 0] = 1.0

    vrtspec, divspec, ierror = poa1_to_vordiv(ntrunc, km, poa1, rsphere)

    assert ierror == 0
    assert vrtspec.shape == (((ntrunc + 1) * (ntrunc + 2)) // 2,)
    assert divspec.shape == (((ntrunc + 1) * (ntrunc + 2)) // 2,)
    assert np.count_nonzero(np.abs(vrtspec) > 0.0) > 0


def test_ldfou2_ledir_poa1_chain_preserves_target_poa1_peak_row():
    nlat = 17
    nlon = 32
    ntrunc = nlat - 1
    km = 2
    rsphere = 6.3712e6
    ndgnh = (nlat + 1) // 2
    ia = 1 + ((ntrunc - km + 2) % 2)
    is_ = 1 + ((ntrunc - km + 1) % 2)

    rpnma, rpnms, pw = build_ledir_blocks_from_scalar_bridge(nlon, nlat, km)
    lats_deg, _ = spharm_pkg.gaussian_lats_wts(nlat)
    scaling = 1.0 / (rsphere * np.cos(np.deg2rad(lats_deg[:ndgnh])))

    target_n = 4
    target_row = ntrunc + 2 - target_n
    paia = np.zeros((4, ndgnh), dtype=np.float64)
    psia = np.zeros((4, ndgnh), dtype=np.float64)

    if target_row % 2 == ia % 2:
        col = (target_row - ia) // 2
        desired_scaled = rpnma[:, col] / pw
        paia[0, :] = desired_scaled / scaling
    else:
        col = (target_row - is_) // 2
        desired_scaled = rpnms[:, col] / pw
        psia[0, :] = desired_scaled / scaling

    paia_scaled, psia_scaled, ierr = ldfou2_uv_scaling(ntrunc, km, paia, rsphere, psia)
    assert ierr == 0

    poa1, ierr = ledir_dgemm(ntrunc, km, paia_scaled, psia_scaled, rpnma, rpnms, pw)
    assert ierr == 0
    assert int(np.argmax(np.abs(poa1[:, 0]))) == target_row - 1

    vrtspec, divspec, ierr = poa1_to_vordiv(ntrunc, km, poa1, rsphere)
    assert ierr == 0
    assert np.count_nonzero(np.abs(vrtspec) > 1e-10) > 0
    assert np.count_nonzero(np.abs(divspec) > 1e-10) > 0

    expected_vrt_peak = _spectral_index(ntrunc, km, target_n + 1)
    expected_div_peak = _spectral_index(ntrunc, km, target_n)
    assert int(np.argmax(np.abs(vrtspec))) == expected_vrt_peak
    assert int(np.argmax(np.abs(divspec))) == expected_div_peak


def test_direct_vordiv_pilot_from_uv_blocks_matches_manual_chain():
    nlat = 17
    nlon = 32
    ntrunc = nlat - 1
    km = 2
    rsphere = 6.3712e6
    ndgnh = (nlat + 1) // 2
    ia = 1 + ((ntrunc - km + 2) % 2)
    is_ = 1 + ((ntrunc - km + 1) % 2)

    rpnma, rpnms, pw_manual = build_ledir_blocks_from_scalar_bridge(nlon, nlat, km)
    lats_deg, _ = spharm_pkg.gaussian_lats_wts(nlat)
    scaling = 1.0 / (rsphere * np.cos(np.deg2rad(lats_deg[:ndgnh])))

    target_n = 4
    target_row = ntrunc + 2 - target_n
    paia = np.zeros((4, ndgnh), dtype=np.float64)
    psia = np.zeros((4, ndgnh), dtype=np.float64)

    if target_row % 2 == ia % 2:
        col = (target_row - ia) // 2
        desired_scaled = rpnma[:, col] / pw_manual
        paia[0, :] = desired_scaled / scaling
    else:
        col = (target_row - is_) // 2
        desired_scaled = rpnms[:, col] / pw_manual
        psia[0, :] = desired_scaled / scaling

    paia_scaled, psia_scaled, ierr = ldfou2_uv_scaling(ntrunc, km, paia, rsphere, psia)
    assert ierr == 0
    poa1_manual, ierr = ledir_dgemm(
        ntrunc, km, paia_scaled, psia_scaled, rpnma, rpnms, pw_manual
    )
    assert ierr == 0
    vrtspec_manual, divspec_manual, ierr = poa1_to_vordiv(
        ntrunc, km, poa1_manual, rsphere
    )
    assert ierr == 0

    poa1_actual, vrtspec_actual, divspec_actual, pw_actual, ierr = (
        direct_vordiv_pilot_from_uv_blocks(nlon, nlat, km, paia, rsphere, psia)
    )
    assert ierr == 0
    np.testing.assert_allclose(pw_actual, pw_manual, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(poa1_actual, poa1_manual, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(vrtspec_actual, vrtspec_manual, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(divspec_actual, divspec_manual, rtol=0.0, atol=1e-12)


def test_direct_vordiv_pilot_sweep_matches_sum_of_single_km_pilots():
    nlat = 9
    nlon = 16
    ntrunc = nlat - 1
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260507)

    uv_blocks_by_km = []
    for km in range(ntrunc + 1):
        ndgnh = (nlat + 1) // 2
        paia = rng.standard_normal((4, ndgnh)).astype(np.float64)
        psia = rng.standard_normal((4, ndgnh)).astype(np.float64)
        uv_blocks_by_km.append((paia, psia))

    poa1_by_km, vrt_total, div_total, pw_total, ierr = (
        direct_vordiv_pilot_sweep_from_uv_blocks(
            nlon,
            nlat,
            uv_blocks_by_km,
            rsphere,
        )
    )
    assert ierr == 0
    assert len(poa1_by_km) == ntrunc + 1
    assert pw_total.shape == ((nlat + 1) // 2,)

    vrt_manual = np.zeros_like(vrt_total)
    div_manual = np.zeros_like(div_total)
    pw_manual = None
    for km, (paia, psia) in enumerate(uv_blocks_by_km):
        poa1, vrt, div, pw, ierr = direct_vordiv_pilot_from_uv_blocks(
            nlon, nlat, km, paia, rsphere, psia
        )
        assert ierr == 0
        assert poa1.shape == poa1_by_km[km].shape
        np.testing.assert_allclose(poa1, poa1_by_km[km], rtol=0.0, atol=1e-12)
        vrt_manual += vrt
        div_manual += div
        if pw_manual is None:
            pw_manual = pw
        else:
            np.testing.assert_allclose(pw, pw_manual, rtol=0.0, atol=0.0)

    np.testing.assert_allclose(vrt_total, vrt_manual, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(div_total, div_manual, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(pw_total, pw_manual, rtol=0.0, atol=0.0)


def test_direct_vordiv_pilot_sweep_validates_length_and_block_shape():
    nlat = 9
    nlon = 16
    rsphere = 6.3712e6
    ndgnh = (nlat + 1) // 2

    with np.testing.assert_raises_regex(ValueError, "must have length 9"):
        direct_vordiv_pilot_sweep_from_uv_blocks(
            nlon,
            nlat,
            [(np.zeros((4, ndgnh)), np.zeros((4, ndgnh)))],
            rsphere,
        )

    bad_blocks = [
        (np.zeros((4, ndgnh), dtype=np.float64), np.zeros((4, ndgnh), dtype=np.float64))
        for _ in range(nlat)
    ]
    bad_blocks[3] = (
        np.zeros((2, ndgnh), dtype=np.float64),
        np.zeros((2, ndgnh), dtype=np.float64),
    )

    with np.testing.assert_raises_regex(ValueError, "first dimension 4"):
        direct_vordiv_pilot_sweep_from_uv_blocks(
            nlon,
            nlat,
            bad_blocks,
            rsphere,
        )


def test_direct_vordiv_pilot_sweep_currently_differs_from_uv_spectral_roundtrip():
    nlat = 9
    nlon = 16
    ntrunc = nlat - 1
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260508)

    uv_blocks_by_km = []
    for km in range(ntrunc + 1):
        ndgnh = (nlat + 1) // 2
        paia = rng.standard_normal((4, ndgnh)).astype(np.float64)
        psia = rng.standard_normal((4, ndgnh)).astype(np.float64)
        uv_blocks_by_km.append((paia, psia))

    _, vrt_direct, div_direct, _, ierr = direct_vordiv_pilot_sweep_from_uv_blocks(
        nlon,
        nlat,
        uv_blocks_by_km,
        rsphere,
    )
    assert ierr == 0

    uspec, vspec, ierr = vordiv_to_uv(vrt_direct, div_direct, rsphere)
    assert ierr == 0
    vrt_back, div_back, ierr = uv_to_vordiv(uspec, vspec, rsphere)
    assert ierr == 0

    vrt_abs = np.abs(vrt_back - vrt_direct)
    div_abs = np.abs(div_back - div_direct)

    assert vrt_abs.max() > 1e-3
    assert div_abs.max() > 1e-3
    assert np.abs(vrt_direct).max() < 1e-4
    assert np.abs(div_direct).max() < 1e-4
    assert np.abs(vrt_back).max() > 1e-1
    assert np.abs(div_back).max() > 1e-1


def test_uv_to_vordiv_linear_pilot_reconstructs_random_uv_spectra():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260510)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0
    vrt_back, div_back = uv_to_vordiv_linear_pilot(uspec, vspec, rsphere)
    uspec_back, vspec_back, ierr = vordiv_to_uv(vrt_back, div_back, rsphere)
    assert ierr == 0

    np.testing.assert_allclose(uspec_back, uspec, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(vspec_back, vspec, rtol=1e-13, atol=2e-8)
    assert np.abs(vrt_back[0]) < 1e-12
    assert np.abs(div_back[0]) < 1e-12


def test_uv_to_vordiv_linear_pilot_reconstructs_direct_sweep_uv_spectra():
    nlat = 9
    nlon = 16
    ntrunc = nlat - 1
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260511)

    uv_blocks_by_km = []
    for km in range(ntrunc + 1):
        ndgnh = (nlat + 1) // 2
        paia = rng.standard_normal((4, ndgnh)).astype(np.float64)
        psia = rng.standard_normal((4, ndgnh)).astype(np.float64)
        uv_blocks_by_km.append((paia, psia))

    _, vrt_direct, div_direct, _, ierr = direct_vordiv_pilot_sweep_from_uv_blocks(
        nlon,
        nlat,
        uv_blocks_by_km,
        rsphere,
    )
    assert ierr == 0
    uspec, vspec, ierr = vordiv_to_uv(vrt_direct, div_direct, rsphere)
    assert ierr == 0

    vrt_back, div_back = uv_to_vordiv_linear_pilot(uspec, vspec, rsphere)
    uspec_back, vspec_back, ierr = vordiv_to_uv(vrt_back, div_back, rsphere)
    assert ierr == 0

    np.testing.assert_allclose(uspec_back, uspec, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(vspec_back, vspec, rtol=1e-13, atol=2e-8)
    assert np.abs(vrt_back[0]) < 1e-12
    assert np.abs(div_back[0]) < 1e-12


def test_uv_to_vordiv_block_linear_pilot_matches_dense_linear_reference():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260512)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    vrt_dense, div_dense = uv_to_vordiv_linear_pilot(uspec, vspec, rsphere)
    vrt_block, div_block = uv_to_vordiv_block_linear_pilot(uspec, vspec, rsphere)

    np.testing.assert_allclose(vrt_block, vrt_dense, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(div_block, div_dense, rtol=1e-13, atol=2e-8)


def test_uv_to_vordiv_block_linear_pilot_reconstructs_random_uv_spectra():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260513)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    vrt_back, div_back = uv_to_vordiv_block_linear_pilot(uspec, vspec, rsphere)
    uspec_back, vspec_back, ierr = vordiv_to_uv(vrt_back, div_back, rsphere)
    assert ierr == 0

    np.testing.assert_allclose(uspec_back, uspec, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(vspec_back, vspec, rtol=1e-13, atol=2e-8)
    assert np.abs(vrt_back[0]) < 1e-12
    assert np.abs(div_back[0]) < 1e-12


def test_uv_to_vordiv_native_still_differs_from_block_linear_reference():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260514)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    vrt_ref, div_ref = uv_to_vordiv_block_linear_pilot(uspec, vspec, rsphere)
    vrt_native, div_native, ierr = uv_to_vordiv(uspec, vspec, rsphere)
    assert ierr == 0

    assert np.max(np.abs(vrt_native - vrt_ref)) > 1e-3
    assert np.max(np.abs(div_native - div_ref)) > 1e-3


def test_uv_to_vordiv_block_native_pilot_matches_block_linear_reference():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260515)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    vrt_ref, div_ref = uv_to_vordiv_block_linear_pilot(uspec, vspec, rsphere)
    vrt_native, div_native = uv_to_vordiv_block_native_pilot(uspec, vspec, rsphere)

    np.testing.assert_allclose(vrt_native, vrt_ref, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(div_native, div_ref, rtol=1e-13, atol=2e-8)


def test_uv_to_vordiv_block_native_pilot_reconstructs_random_uv_spectra():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260516)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    vrt_back, div_back = uv_to_vordiv_block_native_pilot(uspec, vspec, rsphere)
    uspec_back, vspec_back, ierr = vordiv_to_uv(vrt_back, div_back, rsphere)
    assert ierr == 0

    np.testing.assert_allclose(uspec_back, uspec, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(vspec_back, vspec, rtol=1e-13, atol=2e-8)
    assert np.abs(vrt_back[0]) < 1e-12
    assert np.abs(div_back[0]) < 1e-12


def test_uv_to_vordiv_block_native_with_setup_matches_block_linear_reference():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260517)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    handle = create_uv_to_vordiv_block_setup(uspec, vspec, rsphere)
    try:
        vrt_cached, div_cached = uv_to_vordiv_block_native_with_setup(
            handle, uspec, vspec
        )
    finally:
        destroy_uv_to_vordiv_block_setup(handle)

    vrt_ref, div_ref = uv_to_vordiv_block_linear_pilot(uspec, vspec, rsphere)
    np.testing.assert_allclose(vrt_cached, vrt_ref, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(div_cached, div_ref, rtol=1e-13, atol=2e-8)


def test_uv_to_vordiv_block_native_with_setup_reuses_handle_and_destroy_closes_it():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260518)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    handle = create_uv_to_vordiv_block_setup(uspec, vspec, rsphere)
    vrt_first, div_first = uv_to_vordiv_block_native_with_setup(handle, uspec, vspec)
    vrt_second, div_second = uv_to_vordiv_block_native_with_setup(handle, uspec, vspec)

    np.testing.assert_allclose(vrt_second, vrt_first, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(div_second, div_first, rtol=0.0, atol=0.0)

    destroy_uv_to_vordiv_block_setup(handle)

    with np.testing.assert_raises_regex(RuntimeError, "setup handle is closed"):
        uv_to_vordiv_block_native_with_setup(handle, uspec, vspec)


def test_uv_to_vordiv_block_solver_reuses_cached_handle():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260521)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0

    solver = UvToVordivBlockSolver(uspec, vspec, rsphere)
    try:
        vrt_first, div_first = solver.solve(uspec, vspec)
        vrt_second, div_second = solver.solve(uspec, vspec)
    finally:
        solver.close()

    np.testing.assert_allclose(vrt_second, vrt_first, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(div_second, div_first, rtol=0.0, atol=0.0)
    assert solver.closed


def test_uv_to_vordiv_block_solver_context_manager_closes_and_matches_reference():
    nlat = 9
    rsphere = 6.3712e6
    rng = np.random.default_rng(20260522)
    ncoeff = (nlat * (nlat + 1)) // 2

    vrt = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)
    div = rng.standard_normal(ncoeff).astype(np.float64) + 1j * rng.standard_normal(
        ncoeff
    ).astype(np.float64)

    uspec, vspec, ierr = vordiv_to_uv(vrt, div, rsphere)
    assert ierr == 0
    vrt_ref, div_ref = uv_to_vordiv_block_linear_pilot(uspec, vspec, rsphere)

    with UvToVordivBlockSolver(uspec, vspec, rsphere) as solver:
        assert not solver.closed
        assert solver.rsphere == rsphere
        vrt_actual, div_actual = solver.solve(uspec, vspec)

    assert solver.closed
    np.testing.assert_allclose(vrt_actual, vrt_ref, rtol=1e-13, atol=2e-8)
    np.testing.assert_allclose(div_actual, div_ref, rtol=1e-13, atol=2e-8)

    with np.testing.assert_raises_regex(RuntimeError, "solver is closed"):
        solver.solve(uspec, vspec)


def test_reduced_gaussian_spharmt_vector_shape_validation():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    u = np.zeros((116,), dtype=np.float64)
    v = np.zeros((116, 2), dtype=np.float64)
    spec = np.zeros((15,), dtype=np.complex128)
    spec_2d = np.zeros((15, 2), dtype=np.complex128)

    try:
        backend.getvrtdivspec(u, v)
    except Exception as exc:
        assert "same shape" in str(exc)
    else:
        raise AssertionError("Expected same-shape validation error for vector analysis")

    try:
        backend.getuv(spec, spec_2d)
    except Exception as exc:
        assert "same shape" in str(exc)
    else:
        raise AssertionError("Expected same-shape validation error for wind synthesis")


def test_reduced_gaussian_scalar_ectrans_pack_roundtrip():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.array(
        [
            1.0 + 2.0j,
            3.0 + 4.0j,
            5.0 + 6.0j,
            7.0 + 8.0j,
            9.0 + 10.0j,
            11.0 + 12.0j,
        ],
        dtype=np.complex128,
    )

    packed, ntrunc, extra_shape = backend._pack_scalar_spectrum_for_ectrans(
        spec, "test_pack"
    )
    restored = backend._unpack_scalar_spectrum_from_ectrans(packed, "test_unpack")

    assert ntrunc == 2
    assert extra_shape == ()
    assert packed.shape == (12,)
    np.testing.assert_allclose(packed[0::2], spec.real)
    np.testing.assert_allclose(packed[1::2], spec.imag)
    np.testing.assert_allclose(restored, spec)


def test_reduced_gaussian_validate_ectrans_scalar_packed_layout():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    packed = np.zeros((12, 3), dtype=np.float64)

    nt, ntrunc, normalized, extra_shape = backend._validate_ectrans_scalar_packed(
        packed, "packed_layout"
    )

    assert nt == 3
    assert ntrunc == 2
    assert normalized.shape == (12, 3)
    assert extra_shape == (3,)


def test_reduced_gaussian_native_setup_handle_lifecycle():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))

    handle = backend._require_setup_handle()
    assert handle is backend._require_setup_handle()

    backend.close()
    assert backend._setup_handle is None

    handle = backend._require_setup_handle()
    assert handle is not None
    backend.close()


def test_reduced_gaussian_native_setup_layout_metadata():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))

    description = backend._describe_setup()

    assert description["ndgl"] == 5
    assert description["ngptot"] == 116
    assert description["max_nlon"] == 28
    np.testing.assert_array_equal(
        description["nloen"],
        np.array([20, 24, 28, 24, 20], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        description["lat_offsets"],
        np.array([0, 20, 44, 72, 96], dtype=np.int32),
    )
    assert description["mu"].shape == (5,)
    assert description["weights"].shape == (5,)
    assert np.all(np.diff(description["mu"]) < 0.0)
    assert np.all(description["weights"] > 0.0)
    np.testing.assert_allclose(description["weights"].sum(), 2.0, rtol=0.0, atol=1e-14)

    backend.close()


def test_reduced_gaussian_close_is_idempotent():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend.close()
    backend.close()


def test_reduced_gaussian_spectogrd_constant_field_prototype():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.zeros((6,), dtype=np.complex128)
    spec[0] = np.sqrt(2.0)

    grid = backend.spectogrd(spec)

    assert grid.shape == (116,)
    np.testing.assert_allclose(grid, 1.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_native_scalar_synthesis_stub_returns_success():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.zeros((6,), dtype=np.complex128)
    spec[0] = np.sqrt(2.0)
    packed_spec, _, _ = backend._pack_scalar_spectrum_for_ectrans(
        spec, "native stub test"
    )

    datagrid, ierror = scalar_synthesis_stub(
        backend._require_setup_handle(),
        np.asarray(packed_spec, dtype=np.float64),
    )

    assert ierror == 0
    assert datagrid.shape == (116,)
    np.testing.assert_allclose(datagrid, 1.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_native_scalar_fourier_stub_matches_expected_modes():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    lat_offsets = backend._get_lat_offsets()
    grid = np.zeros((backend.ngptot, 2), dtype=np.float64)

    for ilat, nlon in enumerate(backend.nloen):
        offset = int(lat_offsets[ilat])
        lon = 2.0 * np.pi * np.arange(int(nlon), dtype=np.float64) / float(nlon)
        grid[offset : offset + int(nlon), 0] = 1.0
        grid[offset : offset + int(nlon), 1] = np.cos(lon)

    fourier, ierror = scalar_fourier_stub(
        backend._require_setup_handle(),
        grid,
        2,
    )

    assert ierror == 0
    assert fourier.shape == (backend.ndgl, 3, 2)
    np.testing.assert_allclose(fourier[:, 0, 0], 1.0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(fourier[:, 1:, 0], 0.0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(fourier[:, 0, 1], 0.0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(fourier[:, 1, 1], 0.5, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(fourier[:, 2, 1], 0.0, rtol=0.0, atol=1e-12)
    backend.close()


def test_reduced_gaussian_native_scalar_block_solver_matches_numpy_lstsq():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    _, weights, _ = backend._get_gaussian_metadata()
    sqrt_weights = np.sqrt(weights)[:, None]

    basis_block = np.array(
        [
            [1.0 + 0.0j, 0.2 - 0.1j],
            [0.8 + 0.1j, -0.3 + 0.2j],
            [0.5 - 0.2j, 0.4 + 0.3j],
            [0.2 + 0.3j, 0.6 - 0.2j],
            [0.1 - 0.1j, 0.7 + 0.1j],
        ],
        dtype=np.complex128,
    )
    observed = np.array(
        [
            [0.5 + 0.2j, -0.2 + 0.1j],
            [0.3 - 0.1j, 0.1 + 0.2j],
            [-0.2 + 0.4j, 0.3 - 0.3j],
            [0.1 + 0.3j, 0.2 + 0.0j],
            [0.4 - 0.2j, -0.1 + 0.5j],
        ],
        dtype=np.complex128,
    )

    native_solution, ierror = scalar_block_solve_stub(
        backend._require_setup_handle(),
        basis_block,
        observed,
    )
    weighted_basis = sqrt_weights * basis_block
    weighted_observed = sqrt_weights * observed
    numpy_solution, _, _, _ = np.linalg.lstsq(
        weighted_basis,
        weighted_observed,
        rcond=None,
    )

    assert ierror == 0
    np.testing.assert_allclose(native_solution, numpy_solution, rtol=0.0, atol=1e-12)
    backend.close()


def test_reduced_gaussian_native_weighted_block_solver_matches_numpy_lstsq():
    weights = np.array([0.2, 0.4, 0.3, 0.7, 0.5, 0.6], dtype=np.float64)
    sqrt_weights = np.sqrt(weights)[:, None]

    basis_block = np.array(
        [
            [1.0 + 0.0j, 0.2 - 0.1j, -0.1 + 0.3j],
            [0.8 + 0.1j, -0.3 + 0.2j, 0.5 - 0.2j],
            [0.5 - 0.2j, 0.4 + 0.3j, -0.2 + 0.1j],
            [0.2 + 0.3j, 0.6 - 0.2j, 0.1 + 0.4j],
            [0.1 - 0.1j, 0.7 + 0.1j, -0.3 + 0.2j],
            [0.4 + 0.2j, -0.5 + 0.2j, 0.6 - 0.1j],
        ],
        dtype=np.complex128,
    )
    observed = np.array(
        [
            [0.5 + 0.2j, -0.2 + 0.1j],
            [0.3 - 0.1j, 0.1 + 0.2j],
            [-0.2 + 0.4j, 0.3 - 0.3j],
            [0.1 + 0.3j, 0.2 + 0.0j],
            [0.4 - 0.2j, -0.1 + 0.5j],
            [-0.3 + 0.1j, 0.6 - 0.4j],
        ],
        dtype=np.complex128,
    )

    native_solution, ierror = weighted_block_solve_stub(
        weights,
        basis_block,
        observed,
    )
    weighted_basis = sqrt_weights * basis_block
    weighted_observed = sqrt_weights * observed
    numpy_solution, _, _, _ = np.linalg.lstsq(
        weighted_basis,
        weighted_observed,
        rcond=None,
    )

    assert ierror == 0
    np.testing.assert_allclose(native_solution, numpy_solution, rtol=0.0, atol=1e-12)


def test_reduced_gaussian_native_scalar_analysis_stub_returns_success():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grid = np.ones((backend.ngptot,), dtype=np.float64)

    dataspec_packed, ierror = scalar_analysis_stub(
        backend._require_setup_handle(),
        grid,
        2,
    )
    spec = backend._unpack_scalar_spectrum_from_ectrans(
        dataspec_packed,
        "native scalar analysis test",
    )

    assert ierror == 0
    assert dataspec_packed.shape == (12,)
    np.testing.assert_allclose(spec[0], np.sqrt(2.0), rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(spec[1:], 0.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_native_scalar_analysis_matches_prototype_on_nontrivial_field():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    lat_offsets = backend._get_lat_offsets()
    field = np.zeros((backend.ngptot,), dtype=np.float64)

    for ilat, nlon in enumerate(backend.nloen):
        offset = int(lat_offsets[ilat])
        lon = 2.0 * np.pi * np.arange(int(nlon), dtype=np.float64) / float(nlon)
        field[offset : offset + int(nlon)] = (
            (ilat + 1) * 0.1 + 0.7 * np.cos(lon) - 0.2 * np.sin(2.0 * lon)
        )

    _, normalized, extra_shape = backend._validate_packed_grid_data(
        field,
        "native prototype compare",
    )
    prototype = backend._grdtospec_prototype(normalized, 2, extra_shape)
    native = backend.grdtospec(field, ntrunc=2)

    np.testing.assert_allclose(native, prototype, rtol=0.0, atol=1e-8)
    backend.close()


def test_reduced_gaussian_grdtospec_constant_field_prototype():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grid = np.ones((116,), dtype=np.float64)

    spec = backend.grdtospec(grid, ntrunc=2)

    assert spec.shape == (6,)
    np.testing.assert_allclose(spec[0], np.sqrt(2.0), rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(spec[1:], 0.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_grdtospec_recovers_native_basis_vector():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec_in = np.zeros((6,), dtype=np.complex128)
    spec_in[4] = 1.0 + 0.0j

    grid = backend.spectogrd(spec_in)
    spec_out = backend.grdtospec(grid, ntrunc=2)

    np.testing.assert_allclose(spec_out, spec_in, rtol=0.0, atol=1e-5)
    backend.close()


def test_reduced_gaussian_scalar_constant_roundtrip_subset():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec_in = np.zeros((6,), dtype=np.complex128)
    spec_in[0] = np.sqrt(2.0)

    grid = backend.spectogrd(spec_in)
    spec_out = backend.grdtospec(grid, ntrunc=2)

    np.testing.assert_allclose(spec_out[0], spec_in[0], rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(spec_out[1:], 0.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_spectogrd_multi_field_shape_and_values():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    spec = np.zeros((6, 2), dtype=np.complex128)
    spec[0, 0] = np.sqrt(2.0)
    spec[0, 1] = 2.0 * np.sqrt(2.0)

    grid = backend.spectogrd(spec)

    assert grid.shape == (116, 2)
    np.testing.assert_allclose(grid[:, 0], 1.0, rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(grid[:, 1], 2.0, rtol=0.0, atol=1e-6)
    backend.close()


def test_reduced_gaussian_psichi_and_smoothing_match_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()

    base_spec = np.zeros((6,), dtype=np.complex128)
    base_spec[1] = 0.25 - 0.1j
    base_spec[2] = -0.2 + 0.05j
    base_spec[4] = 0.15 + 0.12j

    full_u, full_v = bridge.getgrad(base_spec)
    packed_u = backend._bridge_to_packed_grid(full_u, "test packed u")
    packed_v = backend._bridge_to_packed_grid(full_v, "test packed v")

    psi_spec_expected = bridge.getpsispec(full_u, full_v, ntrunc=2)
    chi_spec_expected = bridge.getchispec(full_u, full_v, ntrunc=2)
    psi_spec_actual = backend.getpsispec(packed_u, packed_v, ntrunc=2)
    chi_spec_actual = backend.getchispec(packed_u, packed_v, ntrunc=2)
    np.testing.assert_allclose(psi_spec_actual, psi_spec_expected, rtol=0.0, atol=5e-8)
    np.testing.assert_allclose(chi_spec_actual, chi_spec_expected, rtol=0.0, atol=5e-8)

    psi_expected = backend._bridge_to_packed_grid(
        bridge.spectogrd(psi_spec_expected),
        "test packed psi",
    )
    chi_expected = backend._bridge_to_packed_grid(
        bridge.spectogrd(chi_spec_expected),
        "test packed chi",
    )
    psi_actual = backend.getpsi(packed_u, packed_v, ntrunc=2)
    chi_actual = backend.getchi(packed_u, packed_v, ntrunc=2)
    both_psi, both_chi = backend.getpsichi(packed_u, packed_v, ntrunc=2)
    np.testing.assert_allclose(psi_actual, psi_expected, rtol=0.0, atol=5e-8)
    np.testing.assert_allclose(chi_actual, chi_expected, rtol=0.0, atol=2e-7)
    np.testing.assert_allclose(both_psi, psi_expected, rtol=0.0, atol=5e-8)
    np.testing.assert_allclose(both_chi, chi_expected, rtol=0.0, atol=2e-7)

    smooth = np.array([1.0, 0.8, 0.6, 0.4, 0.2], dtype=np.float64)
    full_scalar = bridge.spectogrd(base_spec)
    packed_scalar = backend._bridge_to_packed_grid(
        full_scalar, "test packed scalar smooth"
    )
    smoothed_expected = backend._bridge_to_packed_grid(
        bridge.specsmooth(full_scalar, smooth),
        "test packed scalar smooth expected",
    )
    smoothed_actual = backend.specsmooth(packed_scalar, smooth)
    np.testing.assert_allclose(
        smoothed_actual, smoothed_expected, rtol=0.0, atol=1.2e-7
    )

    backend.close()


def test_reduced_gaussian_local_lap_and_invlap_match_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()

    spec = np.zeros((6, 2), dtype=np.complex128)
    spec[1, 0] = 0.25 - 0.1j
    spec[2, 0] = -0.2 + 0.05j
    spec[4, 0] = 0.15 + 0.12j
    spec[0, 1] = 0.5 + 0.0j
    spec[3, 1] = -0.08 + 0.03j
    spec[5, 1] = 0.04 - 0.02j

    lap_actual = backend._lapspec(spec, "test_lap")
    invlap_actual = backend._invlapspec(spec, "test_invlap")
    lap_expected = bridge._restore_spectral_shape(
        bridge._call_spherepack_safely(
            __import__(
                "skyborn.spharm.spherical_harmonics", fromlist=["_spherepack"]
            )._spherepack.lap,
            bridge._validate_spectral_data(spec, "test_lap")[2],
            bridge.rsphere,
            operation_name="bridge lap",
        ),
        (2,),
    )
    invlap_expected = bridge._invlapspec(spec, "test_invlap")

    np.testing.assert_allclose(lap_actual, lap_expected, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(invlap_actual, invlap_expected, rtol=0.0, atol=1e-12)
    backend.close()


def test_reduced_gaussian_psichi_pair_paths_only_analyze_once(monkeypatch):
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    calls = {"getvrtdivspec": 0}

    vrtspec = np.zeros((6,), dtype=np.complex128)
    divspec = np.zeros((6,), dtype=np.complex128)
    vrtspec[1] = 0.3 - 0.1j
    divspec[2] = -0.2 + 0.05j

    def fake_getvrtdivspec(ugrid, vgrid, ntrunc=None):
        calls["getvrtdivspec"] += 1
        return vrtspec, divspec

    monkeypatch.setattr(backend, "getvrtdivspec", fake_getvrtdivspec)

    psi_spec, chi_spec = backend.getpsichispec(
        np.zeros((backend.ngptot,), dtype=np.float64),
        np.zeros((backend.ngptot,), dtype=np.float64),
        ntrunc=2,
    )
    assert calls["getvrtdivspec"] == 1

    psi_grid, chi_grid = backend.getpsichi(
        np.zeros((backend.ngptot,), dtype=np.float64),
        np.zeros((backend.ngptot,), dtype=np.float64),
        ntrunc=2,
    )
    assert calls["getvrtdivspec"] == 2
    np.testing.assert_allclose(psi_spec, backend._invlapspec(vrtspec, "test psi spec"))
    np.testing.assert_allclose(chi_spec, backend._invlapspec(divspec, "test chi spec"))
    assert psi_grid.shape == (backend.ngptot,)
    assert chi_grid.shape == (backend.ngptot,)
    backend.close()


def test_reduced_gaussian_pair_bridge_helpers_match_separate_paths():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()
    lat_offsets = backend._get_lat_offsets()

    packed_a = np.zeros((backend.ngptot,), dtype=np.float64)
    packed_b = np.zeros((backend.ngptot,), dtype=np.float64)
    for ilat, nlon in enumerate(backend.nloen):
        offset = int(lat_offsets[ilat])
        lon = 2.0 * np.pi * np.arange(int(nlon), dtype=np.float64) / float(nlon)
        packed_a[offset : offset + int(nlon)] = 0.5 * (ilat + 1) + np.cos(lon)
        packed_b[offset : offset + int(nlon)] = -0.25 * (ilat + 1) + np.sin(2.0 * lon)

    _, normalized_a, extra_shape = backend._validate_packed_grid_data(
        packed_a,
        "pair bridge a",
    )
    _, normalized_b, _ = backend._validate_packed_grid_data(packed_b, "pair bridge b")

    full_a_expected = backend._packed_to_bridge_grid(normalized_a, extra_shape)
    full_b_expected = backend._packed_to_bridge_grid(normalized_b, extra_shape)
    full_a_actual, full_b_actual = backend._packed_pair_to_bridge_grids(
        normalized_a,
        normalized_b,
        extra_shape,
    )
    np.testing.assert_allclose(full_a_actual, full_a_expected, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(full_b_actual, full_b_expected, rtol=0.0, atol=1e-12)

    packed_a_expected = backend._bridge_to_packed_grid(
        full_a_expected, "pair bridge out"
    )
    packed_b_expected = backend._bridge_to_packed_grid(
        full_b_expected, "pair bridge out"
    )
    packed_a_actual, packed_b_actual = backend._bridge_pair_to_packed_grids(
        full_a_expected,
        full_b_expected,
        "pair bridge out",
    )
    np.testing.assert_allclose(packed_a_actual, packed_a_expected, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(packed_b_actual, packed_b_expected, rtol=0.0, atol=1e-12)
    assert (
        full_a_actual.shape
        == bridge.spectogrd(np.zeros((6,), dtype=np.complex128)).shape
    )
    backend.close()


def test_reduced_gaussian_vector_paths_skip_unimplemented_native_stubs(monkeypatch):
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend._native_uv_synthesis_available = False
    backend._native_vrtdiv_analysis_available = False
    backend._native_gradient_synthesis_available = False

    def fail_if_called(*args, **kwargs):
        raise AssertionError("native vector stub should not be called when disabled")

    monkeypatch.setattr(
        "skyborn.spharm.reduced_gaussian.uv_synthesis_stub",
        fail_if_called,
    )
    monkeypatch.setattr(
        "skyborn.spharm.reduced_gaussian.vrtdiv_analysis_stub",
        fail_if_called,
    )
    monkeypatch.setattr(
        "skyborn.spharm.reduced_gaussian.gradient_synthesis_stub",
        fail_if_called,
    )

    bridge = backend._get_bridge_spharmt()
    spec = np.zeros((6,), dtype=np.complex128)
    spec[1] = 0.25 - 0.1j
    spec[2] = -0.2 + 0.05j
    spec[4] = 0.15 + 0.12j
    full_u, full_v = bridge.getgrad(spec)
    packed_u = backend._bridge_to_packed_grid(full_u, "skip native u")
    packed_v = backend._bridge_to_packed_grid(full_v, "skip native v")
    vrt, div = bridge.getvrtdivspec(full_u, full_v, ntrunc=2)

    u_actual, v_actual = backend.getuv(vrt, div)
    gu_actual, gv_actual = backend.getgrad(spec)
    vrt_actual, div_actual = backend.getvrtdivspec(packed_u, packed_v, ntrunc=2)

    np.testing.assert_allclose(u_actual, packed_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(v_actual, packed_v, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(gu_actual, packed_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(gv_actual, packed_v, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(vrt_actual, vrt, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(div_actual, div, rtol=0.0, atol=1e-10)
    backend.close()


def test_reduced_gaussian_native_uv_synthesis_matches_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend._native_uv_synthesis_available = True
    backend._native_vrtdiv_analysis_available = False
    backend._native_gradient_synthesis_available = False
    bridge = backend._get_bridge_spharmt()

    spec = np.zeros((6,), dtype=np.complex128)
    spec[1] = 0.25 - 0.1j
    spec[2] = -0.2 + 0.05j
    spec[4] = 0.15 + 0.12j

    full_u, full_v = bridge.getgrad(spec)
    vrtspec, divspec = bridge.getvrtdivspec(full_u, full_v, ntrunc=2)
    packed_u_expected = backend._bridge_to_packed_grid(full_u, "native uv expected u")
    packed_v_expected = backend._bridge_to_packed_grid(full_v, "native uv expected v")

    packed_u_actual, packed_v_actual = backend.getuv(vrtspec, divspec)

    np.testing.assert_allclose(
        packed_u_actual,
        packed_u_expected,
        rtol=0.0,
        atol=2.5e-7,
    )
    np.testing.assert_allclose(
        packed_v_actual,
        packed_v_expected,
        rtol=0.0,
        atol=2.5e-7,
    )
    backend.close()


def test_reduced_gaussian_native_gradient_synthesis_matches_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend._native_uv_synthesis_available = True
    backend._native_vrtdiv_analysis_available = False
    backend._native_gradient_synthesis_available = True
    bridge = backend._get_bridge_spharmt()

    spec = np.zeros((6,), dtype=np.complex128)
    spec[1] = 0.25 - 0.1j
    spec[2] = -0.2 + 0.05j
    spec[4] = 0.15 + 0.12j

    full_u_expected, full_v_expected = bridge.getgrad(spec)
    packed_u_expected = backend._bridge_to_packed_grid(
        full_u_expected,
        "native grad expected u",
    )
    packed_v_expected = backend._bridge_to_packed_grid(
        full_v_expected,
        "native grad expected v",
    )

    packed_u_actual, packed_v_actual = backend.getgrad(spec)

    np.testing.assert_allclose(
        packed_u_actual,
        packed_u_expected,
        rtol=0.0,
        atol=2.5e-7,
    )
    np.testing.assert_allclose(
        packed_v_actual,
        packed_v_expected,
        rtol=0.0,
        atol=2.5e-7,
    )
    backend.close()


def test_reduced_gaussian_native_vrtdiv_analysis_matches_bridge_reference():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    backend._native_uv_synthesis_available = True
    backend._native_vrtdiv_analysis_available = True
    backend._native_gradient_synthesis_available = True
    bridge = backend._get_bridge_spharmt()

    spec = np.zeros((6,), dtype=np.complex128)
    spec[1] = 0.25 - 0.1j
    spec[2] = -0.2 + 0.05j
    spec[4] = 0.15 + 0.12j

    full_u, full_v = bridge.getgrad(spec)
    packed_u = backend._bridge_to_packed_grid(full_u, "native vrtdiv input u")
    packed_v = backend._bridge_to_packed_grid(full_v, "native vrtdiv input v")
    vrt_expected, div_expected = bridge.getvrtdivspec(full_u, full_v, ntrunc=2)

    vrt_actual, div_actual = backend.getvrtdivspec(packed_u, packed_v, ntrunc=2)

    np.testing.assert_allclose(
        vrt_actual,
        vrt_expected,
        rtol=0.0,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        div_actual,
        div_expected,
        rtol=0.0,
        atol=1e-10,
    )
    backend.close()


def test_reduced_gaussian_regrid_rejects_bad_input_and_smooth_shape():
    grdin = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grdout = ReducedGaussianSpharmt(np.array([16, 20, 24, 20, 16], dtype=np.int32))

    with np.testing.assert_raises_regex(
        ReducedGaussianValidationError,
        "leading dimension",
    ):
        reduced_regrid(grdin, grdout, np.zeros((115,), dtype=np.float64))

    with np.testing.assert_raises_regex(
        ReducedGaussianValidationError,
        "smooth must be rank 1 with size 5",
    ):
        reduced_regrid(
            grdin,
            grdout,
            np.zeros((grdin.ngptot,), dtype=np.float64),
            smooth=np.ones((4,), dtype=np.float32),
        )

    with np.testing.assert_raises_regex(
        ReducedGaussianValidationError,
        "ntrunc must be between 0 and 4",
    ):
        reduced_regrid(
            grdin,
            grdout,
            np.zeros((grdin.ngptot,), dtype=np.float64),
            ntrunc=5,
        )

    grdin.close()
    grdout.close()


def test_reduced_gaussian_regrid_identity_matches_direct_roundtrip():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    lat_offsets = backend._get_lat_offsets()
    field = np.zeros((backend.ngptot,), dtype=np.float64)

    for ilat, nlon in enumerate(backend.nloen):
        offset = int(lat_offsets[ilat])
        lon = 2.0 * np.pi * np.arange(int(nlon), dtype=np.float64) / float(nlon)
        field[offset : offset + int(nlon)] = (
            0.5 * (ilat + 1) + 0.75 * np.cos(lon) - 0.2 * np.sin(2.0 * lon)
        )

    direct = backend.spectogrd(backend.grdtospec(field))
    regridded = reduced_regrid(backend, backend, field)

    np.testing.assert_allclose(regridded, direct, rtol=0.0, atol=1e-10)
    backend.close()


def test_reduced_gaussian_regrid_smooth_matches_manual_spectral_chain():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    lat_offsets = backend._get_lat_offsets()
    field = np.zeros((backend.ngptot,), dtype=np.float64)

    for ilat, nlon in enumerate(backend.nloen):
        offset = int(lat_offsets[ilat])
        lon = 2.0 * np.pi * np.arange(int(nlon), dtype=np.float64) / float(nlon)
        field[offset : offset + int(nlon)] = (
            0.2 * (ilat + 1) + np.cos(lon) - 0.35 * np.sin(2.0 * lon)
        )

    smooth = np.array([1.0, 0.75, 0.4, 0.15, 0.0], dtype=np.float32)
    manual = backend.spectogrd(
        backend._apply_spectral_smoothing(
            backend.grdtospec(field, ntrunc=4),
            smooth,
            operation_name="manual regrid smooth",
            expected_size=backend.ndgl,
        )
    )
    actual = reduced_regrid(backend, backend, field, ntrunc=4, smooth=smooth)

    np.testing.assert_allclose(actual, manual, rtol=0.0, atol=1e-10)
    backend.close()


def test_reduced_gaussian_regrid_between_packed_layouts_matches_spectral_reference():
    grdin = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grdout = ReducedGaussianSpharmt(
        np.array([12, 16, 20, 24, 20, 16, 12], dtype=np.int32)
    )
    bridge_in = grdin._get_bridge_spharmt()
    bridge_out = grdout._get_bridge_spharmt()
    lat_offsets = grdin._get_lat_offsets()

    field = np.zeros((grdin.ngptot,), dtype=np.float64)
    for ilat, nlon in enumerate(grdin.nloen):
        offset = int(lat_offsets[ilat])
        lon = 2.0 * np.pi * np.arange(int(nlon), dtype=np.float64) / float(nlon)
        field[offset : offset + int(nlon)] = (
            1.0
            + 0.3 * np.cos(lon)
            - 0.15 * np.sin(2.0 * lon)
            + 0.05 * np.cos(3.0 * lon)
            + 0.1 * ilat
        )

    actual = reduced_regrid(grdin, grdout, field)

    full_in = grdin._bridge_to_packed_grid(
        grdin._packed_to_bridge_grid(
            grdin._validate_packed_grid_data(field, "regrid field")[1],
            (),
        ),
        "regrid bridge roundtrip",
    )
    np.testing.assert_allclose(full_in, field, rtol=0.0, atol=5e-8)

    bridge_field = grdin._packed_to_bridge_grid(
        grdin._validate_packed_grid_data(field, "regrid field")[1],
        (),
    )
    ref_spec = bridge_in.grdtospec(bridge_field, ntrunc=4)
    ref_full = bridge_out.spectogrd(ref_spec)
    expected = grdout._bridge_to_packed_grid(ref_full, "regrid expected")

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=4e-7)
    assert actual.shape == (grdout.ngptot,)
    grdin.close()
    grdout.close()


def test_reduced_gaussian_regriduv_rejects_bad_shapes_and_ntrunc():
    grdin = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grdout = ReducedGaussianSpharmt(np.array([16, 20, 24, 20, 16], dtype=np.int32))
    u = np.zeros((grdin.ngptot,), dtype=np.float64)
    v = np.zeros((grdin.ngptot, 2), dtype=np.float64)

    with np.testing.assert_raises_regex(
        ReducedGaussianValidationError,
        "same shape",
    ):
        reduced_regriduv(grdin, grdout, u, v)

    with np.testing.assert_raises_regex(
        ReducedGaussianValidationError,
        "between 0 and 4",
    ):
        reduced_regriduv(
            grdin,
            grdout,
            np.zeros((grdin.ngptot,), dtype=np.float64),
            np.zeros((grdin.ngptot,), dtype=np.float64),
            ntrunc=5,
        )

    grdin.close()
    grdout.close()


def test_reduced_gaussian_regriduv_identity_matches_direct_vrtdiv_roundtrip():
    backend = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    bridge = backend._get_bridge_spharmt()

    chi_spec = np.zeros((6,), dtype=np.complex128)
    chi_spec[1] = 0.25 - 0.1j
    chi_spec[2] = -0.2 + 0.05j
    chi_spec[4] = 0.15 + 0.12j

    psi_spec = np.zeros((6,), dtype=np.complex128)
    psi_spec[1] = 0.21 - 0.04j
    psi_spec[3] = -0.13 + 0.07j
    psi_spec[4] = 0.09 - 0.02j

    divspec = _spherepack.lap(chi_spec.reshape(6, 1), backend.rsphere).reshape(
        6,
    )
    vrtspec = _spherepack.lap(psi_spec.reshape(6, 1), backend.rsphere).reshape(
        6,
    )
    full_u, full_v = bridge.getuv(vrtspec, divspec)
    packed_u = backend._bridge_to_packed_grid(full_u, "regriduv identity u")
    packed_v = backend._bridge_to_packed_grid(full_v, "regriduv identity v")

    expected_u, expected_v = backend.getuv(
        *backend.getvrtdivspec(packed_u, packed_v, ntrunc=2)
    )
    actual_u, actual_v = reduced_regriduv(
        backend, backend, packed_u, packed_v, ntrunc=2
    )

    np.testing.assert_allclose(actual_u, expected_u, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(actual_v, expected_v, rtol=0.0, atol=1e-10)
    backend.close()


def test_reduced_gaussian_regriduv_between_layouts_matches_bridge_reference():
    grdin = ReducedGaussianSpharmt(np.array([20, 24, 28, 24, 20], dtype=np.int32))
    grdout = ReducedGaussianSpharmt(
        np.array([12, 16, 20, 24, 20, 16, 12], dtype=np.int32)
    )
    bridge_in = grdin._get_bridge_spharmt()
    bridge_out = grdout._get_bridge_spharmt()

    chi_spec = np.zeros((6,), dtype=np.complex128)
    chi_spec[1] = 0.25 - 0.1j
    chi_spec[2] = -0.2 + 0.05j
    chi_spec[4] = 0.15 + 0.12j

    psi_spec = np.zeros((6,), dtype=np.complex128)
    psi_spec[1] = 0.21 - 0.04j
    psi_spec[3] = -0.13 + 0.07j
    psi_spec[4] = 0.09 - 0.02j

    divspec = _spherepack.lap(chi_spec.reshape(6, 1), grdin.rsphere).reshape(
        6,
    )
    vrtspec = _spherepack.lap(psi_spec.reshape(6, 1), grdin.rsphere).reshape(
        6,
    )
    full_u, full_v = bridge_in.getuv(vrtspec, divspec)
    packed_u = grdin._bridge_to_packed_grid(full_u, "regriduv input u")
    packed_v = grdin._bridge_to_packed_grid(full_v, "regriduv input v")

    actual_u, actual_v = reduced_regriduv(grdin, grdout, packed_u, packed_v, ntrunc=2)

    ref_vrt, ref_div = bridge_in.getvrtdivspec(full_u, full_v, ntrunc=2)
    ref_u_full, ref_v_full = bridge_out.getuv(ref_vrt, ref_div)
    expected_u = grdout._bridge_to_packed_grid(ref_u_full, "regriduv expected u")
    expected_v = grdout._bridge_to_packed_grid(ref_v_full, "regriduv expected v")

    np.testing.assert_allclose(actual_u, expected_u, rtol=0.0, atol=3e-7)
    np.testing.assert_allclose(actual_v, expected_v, rtol=0.0, atol=3e-7)
    assert actual_u.shape == (grdout.ngptot,)
    assert actual_v.shape == (grdout.ngptot,)
    grdin.close()
    grdout.close()
