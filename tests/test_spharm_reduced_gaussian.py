"""Tests for packed reduced-Gaussian scalar spherical harmonic transforms."""

import sys
from pathlib import Path

import numpy as np
import pytest

src_path = str(Path(__file__).parent.parent / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

try:
    from skyborn.spharm import ReducedGaussianSpharmt, Spharmt

    SPHARM_AVAILABLE = True
except Exception:
    ReducedGaussianSpharmt = None
    Spharmt = None
    SPHARM_AVAILABLE = False


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_matches_full_gaussian_scalar_analysis_and_synthesis():
    nlat, nlon = 9, 24
    pl = np.full(nlat, nlon, dtype=np.int32)
    rng = np.random.default_rng(20260509)
    grid = rng.standard_normal((nlat, nlon, 2)).astype(np.float32)

    regular = Spharmt(nlon, nlat, gridtype="gaussian", legfunc="stored")
    reduced = ReducedGaussianSpharmt(pl)

    packed = grid.reshape(nlat * nlon, 2)
    spec_regular = regular.grdtospec(grid, ntrunc=6)
    spec_reduced = reduced.grdtospec(packed, ntrunc=6)

    np.testing.assert_allclose(spec_reduced, spec_regular, rtol=2e-5, atol=2e-6)

    grid_regular = regular.spectogrd(spec_regular).reshape(nlat * nlon, 2)
    grid_reduced = reduced.spectogrd(spec_regular)

    np.testing.assert_allclose(grid_reduced, grid_regular, rtol=2e-5, atol=2e-6)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_low_truncation_roundtrip_on_varying_pl():
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    ntrunc = 5
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(20260509)
    spec = (
        rng.standard_normal((nspec, 3)) + 1j * rng.standard_normal((nspec, 3))
    ).astype(np.complex64)
    spec[: ntrunc + 1] = spec[: ntrunc + 1].real.astype(np.complex64)

    reduced = ReducedGaussianSpharmt(pl)
    packed = reduced.spectogrd(spec)
    spec_back = reduced.grdtospec(packed, ntrunc=ntrunc)

    np.testing.assert_allclose(spec_back, spec, rtol=3e-5, atol=3e-5)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_matches_full_gaussian_vector_analysis_and_synthesis():
    nlat, nlon = 9, 24
    pl = np.full(nlat, nlon, dtype=np.int32)
    rng = np.random.default_rng(20260509)
    u = rng.standard_normal((nlat, nlon, 2)).astype(np.float32)
    v = rng.standard_normal((nlat, nlon, 2)).astype(np.float32)

    regular = Spharmt(nlon, nlat, gridtype="gaussian", legfunc="stored")
    reduced = ReducedGaussianSpharmt(pl)

    vrt_regular, div_regular = regular.getvrtdivspec(u, v, ntrunc=6)
    vrt_reduced, div_reduced = reduced.getvrtdivspec(
        u.reshape(nlat * nlon, 2), v.reshape(nlat * nlon, 2), ntrunc=6
    )
    vrt_only = reduced.getvrtspec(
        u.reshape(nlat * nlon, 2), v.reshape(nlat * nlon, 2), ntrunc=6
    )
    div_only = reduced.getdivspec(
        u.reshape(nlat * nlon, 2), v.reshape(nlat * nlon, 2), ntrunc=6
    )

    np.testing.assert_allclose(vrt_reduced, vrt_regular, rtol=2e-6, atol=1e-12)
    np.testing.assert_allclose(div_reduced, div_regular, rtol=2e-6, atol=1e-12)
    np.testing.assert_allclose(vrt_only, vrt_regular, rtol=2e-6, atol=1e-12)
    np.testing.assert_allclose(div_only, div_regular, rtol=2e-6, atol=1e-12)

    u_regular, v_regular = regular.getuv(vrt_regular, div_regular)
    u_reduced, v_reduced = reduced.getuv(vrt_regular, div_regular)

    np.testing.assert_allclose(
        u_reduced, u_regular.reshape(nlat * nlon, 2), rtol=2e-6, atol=1e-6
    )
    np.testing.assert_allclose(
        v_reduced, v_regular.reshape(nlat * nlon, 2), rtol=2e-6, atol=1e-6
    )


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_vector_low_truncation_roundtrip_on_varying_pl():
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    ntrunc = 5
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(20260510)
    vrt = (
        rng.standard_normal((nspec, 2)) + 1j * rng.standard_normal((nspec, 2))
    ).astype(np.complex64)
    div = (
        rng.standard_normal((nspec, 2)) + 1j * rng.standard_normal((nspec, 2))
    ).astype(np.complex64)
    vrt[: ntrunc + 1] = vrt[: ntrunc + 1].real.astype(np.complex64)
    div[: ntrunc + 1] = div[: ntrunc + 1].real.astype(np.complex64)
    vrt[0] = 0.0
    div[0] = 0.0

    reduced = ReducedGaussianSpharmt(pl)
    u, v = reduced.getuv(vrt, div)
    vrt_back, div_back = reduced.getvrtdivspec(u, v, ntrunc=ntrunc)
    vrt_only = reduced.getvrtspec(u, v, ntrunc=ntrunc)
    div_only = reduced.getdivspec(u, v, ntrunc=ntrunc)

    np.testing.assert_allclose(vrt_back, vrt, rtol=3e-5, atol=3e-5)
    np.testing.assert_allclose(div_back, div, rtol=3e-5, atol=3e-5)
    np.testing.assert_allclose(vrt_only, vrt, rtol=3e-5, atol=3e-5)
    np.testing.assert_allclose(div_only, div, rtol=3e-5, atol=3e-5)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_single_component_analysis_skips_paired_helper(
    monkeypatch,
):
    nlat, nlon = 9, 24
    pl = np.full(nlat, nlon, dtype=np.int32)
    rng = np.random.default_rng(20260516)
    u = rng.standard_normal((nlat * nlon, 2)).astype(np.float32)
    v = rng.standard_normal((nlat * nlon, 2)).astype(np.float32)
    reduced = ReducedGaussianSpharmt(pl)

    expected_vrt, expected_div = reduced.getvrtdivspec(u, v, ntrunc=6)
    monkeypatch.setattr(
        reduced,
        "getvrtdivspec",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("single-component analysis should not use getvrtdivspec")
        ),
    )

    vrt_only = reduced.getvrtspec(u, v, ntrunc=6)
    div_only = reduced.getdivspec(u, v, ntrunc=6)

    np.testing.assert_allclose(vrt_only, expected_vrt, rtol=0, atol=1e-14)
    np.testing.assert_allclose(div_only, expected_div, rtol=0, atol=1e-14)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_getgrad_pair_matches_individual_calls():
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    ntrunc = 5
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(20260511)
    spec_a = (
        rng.standard_normal((nspec, 2)) + 1j * rng.standard_normal((nspec, 2))
    ).astype(np.complex64)
    spec_b = (
        rng.standard_normal((nspec, 2)) + 1j * rng.standard_normal((nspec, 2))
    ).astype(np.complex64)
    spec_a[: ntrunc + 1] = spec_a[: ntrunc + 1].real.astype(np.complex64)
    spec_b[: ntrunc + 1] = spec_b[: ntrunc + 1].real.astype(np.complex64)

    reduced = ReducedGaussianSpharmt(pl)
    a_u, a_v = reduced.getgrad(spec_a)
    b_u, b_v = reduced.getgrad(spec_b)
    pair_a_u, pair_a_v, pair_b_u, pair_b_v = reduced.getgrad_pair(spec_a, spec_b)

    np.testing.assert_allclose(pair_a_u, a_u, rtol=2e-6, atol=1e-12)
    np.testing.assert_allclose(pair_a_v, a_v, rtol=2e-6, atol=1e-12)
    np.testing.assert_allclose(pair_b_u, b_u, rtol=2e-6, atol=1e-12)
    np.testing.assert_allclose(pair_b_v, b_v, rtol=2e-6, atol=1e-12)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_getgrad_prefers_single_component_kernel(monkeypatch):
    import skyborn.spharm.reduced_gaussian as reduced_module

    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    ntrunc = 5
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(20260518)
    spec = (
        rng.standard_normal((nspec, 1)) + 1j * rng.standard_normal((nspec, 1))
    ).astype(np.complex64)
    spec[: ntrunc + 1] = spec[: ntrunc + 1].real.astype(np.complex64)

    reduced = ReducedGaussianSpharmt(pl)
    direct_calls = []

    def fake_getgrad(
        spec_a,
        pl_arg,
        basis,
        dbasis,
        sin_theta,
        ntrunc_arg,
        ngptot,
        rsphere,
    ):
        direct_calls.append(spec_a.shape)
        nt = spec_a.shape[1]
        zeros = np.zeros((ngptot, nt), dtype=np.float32, order="F")
        return zeros, zeros.copy(order="F"), 0

    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_getgrad",
        fake_getgrad,
    )
    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_getgrad_pair",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("getgrad should prefer reduced_gaussian_getgrad")
        ),
    )

    reduced.getgrad(spec)
    reduced.getgrad(spec)

    assert direct_calls == [(nspec, 1), (nspec, 1)]


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_spectogrd_pair_matches_individual_calls():
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    ntrunc = 5
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(20260512)
    spec_a = (
        rng.standard_normal((nspec, 2)) + 1j * rng.standard_normal((nspec, 2))
    ).astype(np.complex64)
    spec_b = (
        rng.standard_normal((nspec, 2)) + 1j * rng.standard_normal((nspec, 2))
    ).astype(np.complex64)
    spec_a[: ntrunc + 1] = spec_a[: ntrunc + 1].real.astype(np.complex64)
    spec_b[: ntrunc + 1] = spec_b[: ntrunc + 1].real.astype(np.complex64)

    reduced = ReducedGaussianSpharmt(pl)
    grid_a = reduced.spectogrd(spec_a)
    grid_b = reduced.spectogrd(spec_b)
    pair_a, pair_b = reduced._spectogrd_pair(spec_a, spec_b)

    np.testing.assert_allclose(pair_a, grid_a, rtol=0, atol=0)
    np.testing.assert_allclose(pair_b, grid_b, rtol=0, atol=0)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_shape_validation():
    reduced = ReducedGaussianSpharmt(np.array([8, 12, 8], dtype=np.int32))

    with pytest.raises(ValueError, match="packed first dimension"):
        reduced.grdtospec(np.zeros(10, dtype=np.float32))

    with pytest.raises(ValueError, match="invalid triangular spectral size"):
        reduced.spectogrd(np.zeros(5, dtype=np.complex64))
