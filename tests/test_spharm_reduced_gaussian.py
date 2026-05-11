"""Tests for packed reduced-Gaussian scalar spherical harmonic transforms."""

import sys
from pathlib import Path

import numpy as np
import pytest

src_path = str(Path(__file__).parent.parent / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

try:
    from skyborn.spharm import ReducedGaussianSpharmt, Spharmt, SpheremackError

    SPHARM_AVAILABLE = True
except Exception:
    ReducedGaussianSpharmt = None
    Spharmt = None
    SpheremackError = None
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
def test_reduced_gaussian_precision_aliases_and_float64_inputs():
    pl = np.array([8, 10, 12, 10, 8], dtype=np.int32)
    reduced = ReducedGaussianSpharmt(pl, precision="float64")
    field = np.random.default_rng(0).standard_normal(int(pl.sum())).astype(np.float64)

    spec = reduced.grdtospec(field, ntrunc=2)
    back = reduced.spectogrd(spec)

    assert reduced.precision == "double"
    assert spec.dtype == np.complex128
    assert back.dtype == np.float64


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


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_validation_cache_and_restore_branches(monkeypatch):
    import skyborn.spharm.reduced_gaussian as reduced_module

    cache = {}
    first = np.array([1.0], dtype=np.float32)
    reduced_module._cache_basis_array(cache, (1, 1), first, max_entries=1)
    reduced_module._cache_basis_array(
        cache, (1, 1), np.array([2.0], dtype=np.float32), max_entries=1
    )
    assert cache[(1, 1)] is first
    reduced_module._cache_basis_array(
        cache, (2, 2), np.array([3.0], dtype=np.float32), max_entries=1
    )
    assert (1, 1) not in cache

    reduced_module._GLOBAL_GAUSSIAN_GEOMETRY_CACHE.pop(7, None)
    geometry = reduced_module._cached_gaussian_geometry(7)
    cached_geometry = reduced_module._cached_gaussian_geometry(7)
    assert cached_geometry[0] is geometry[0]
    assert cached_geometry[1] is geometry[1]

    with pytest.raises(ValueError, match="rank 1"):
        ReducedGaussianSpharmt(np.array([[8, 8, 8]], dtype=np.int32))
    with pytest.raises(ValueError, match="at least 3 rows"):
        ReducedGaussianSpharmt(np.array([8, 8], dtype=np.int32))
    with pytest.raises(ValueError, match=">= 4 points"):
        ReducedGaussianSpharmt(np.array([8, 3, 8], dtype=np.int32))
    with pytest.raises(ValueError, match='legfunc must be "stored" or "computed"'):
        ReducedGaussianSpharmt(np.array([8, 10, 8], dtype=np.int32), legfunc="bad")
    with pytest.raises(ValueError, match='precision must be "auto"'):
        ReducedGaussianSpharmt(np.array([8, 10, 8], dtype=np.int32), precision="bad")
    with pytest.raises(ValueError, match="rsphere must be positive"):
        ReducedGaussianSpharmt(np.array([8, 10, 8], dtype=np.int32), rsphere=0.0)

    reduced = ReducedGaussianSpharmt(np.array([8, 10, 8], dtype=np.int32))
    assert "ReducedGaussianSpharmt" in repr(reduced)
    reduced.close()
    reduced.temporary = 1
    del reduced.temporary
    with pytest.raises(AttributeError, match="rebind"):
        reduced.nlat = 99
    with pytest.raises(AttributeError, match="unbind"):
        del reduced.nlat
    with pytest.raises(ValueError, match="ntrunc must be between"):
        reduced._validate_ntrunc(-1)
    with pytest.raises(ValueError, match="at least a rank 1 packed array"):
        reduced._validate_grid_data(np.array(1.0), "op")
    with pytest.raises(ValueError, match="at least a rank 1 spectrum"):
        reduced._validate_spectral_data(np.array(1.0 + 0.0j), "op")

    too_large = np.zeros(10, dtype=np.complex64)
    with pytest.raises(ValueError, match="ntrunc too large"):
        reduced.spectogrd(too_large)

    spec = np.ones((3, 1), dtype=np.complex64, order="F")
    grid = np.ones((reduced.npoints, 1), dtype=np.float32, order="F")
    assert reduced._restore_spectral_shape(spec, (), None).dtype == np.complex64
    assert reduced._restore_grid_shape(grid, (), None).dtype == np.float32

    reduced_single = ReducedGaussianSpharmt(
        np.array([8, 10, 8], dtype=np.int32), precision="single"
    )
    assert reduced_single._public_real_dtype(np.ones(1, dtype=np.float64)) == np.dtype(
        np.float32
    )
    assert reduced_single._restore_spectral_shape(spec, (), None).dtype == np.complex64
    assert reduced_single._restore_grid_shape(grid, (), None).dtype == np.float32

    reduced_double = ReducedGaussianSpharmt(
        np.array([8, 10, 8], dtype=np.int32), precision="double"
    )
    assert reduced_double._restore_spectral_shape(spec, (), None).dtype == np.complex128
    assert reduced_double._restore_grid_shape(grid, (), None).dtype == np.float64
    assert reduced_double._public_real_dtype(np.ones(1, dtype=np.float32)) == np.dtype(
        np.float64
    )

    reduced_auto = ReducedGaussianSpharmt(np.array([8, 10, 8], dtype=np.int32))
    assert reduced_auto._public_real_dtype(np.ones(1, dtype=np.complex128)) == np.dtype(
        np.float64
    )
    assert reduced_auto._public_complex_dtype(
        np.ones(1, dtype=np.complex128)
    ) == np.dtype(np.complex128)

    original_validate = reduced._validate_grid_data

    grid_validate_calls = {"count": 0}

    def mismatched_grid_extra(data, operation_name):
        grid_validate_calls["count"] += 1
        nt, normalized, extra_shape = original_validate(data, operation_name)
        if grid_validate_calls["count"] == 2:
            return nt, normalized, extra_shape + (1,)
        return nt, normalized, extra_shape

    with pytest.raises(ValueError, match="same shape"):
        reduced._validate_paired_grid_data(grid, grid[:, :0], "op")

    monkeypatch.setattr(reduced, "_validate_grid_data", mismatched_grid_extra)
    with pytest.raises(ValueError, match="consistent dimensions"):
        reduced._validate_paired_grid_data(grid, grid.copy(order="F"), "op")

    with pytest.raises(ValueError, match="same shape"):
        reduced._validate_paired_spectral_data(spec, spec[:2], "op")

    original_validate_spec = reduced._validate_spectral_data

    spec_validate_calls = {"count": 0}

    def mismatched_spec_extra(data, operation_name):
        spec_validate_calls["count"] += 1
        nt, ntrunc, normalized, extra_shape = original_validate_spec(
            data, operation_name
        )
        if spec_validate_calls["count"] == 2:
            extra_shape = extra_shape + (1,)
        return nt, ntrunc, normalized, extra_shape

    monkeypatch.setattr(reduced, "_validate_spectral_data", mismatched_spec_extra)
    with pytest.raises(ValueError, match="consistent dimensions"):
        reduced._validate_paired_spectral_data(spec, spec.copy(order="F"), "op")


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_basis_error_and_fallback_branches(monkeypatch):
    import skyborn.spharm.reduced_gaussian as reduced_module

    pl = np.array([8, 10, 8], dtype=np.int32)
    reduced = ReducedGaussianSpharmt(pl)
    key = (reduced.nlat, 1)

    reduced_module._GLOBAL_BASIS_CACHE.pop(key, None)
    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_legendre_basis",
        lambda *args: (_ for _ in ()).throw(RuntimeError("basis boom")),
    )
    with pytest.raises(SpheremackError, match="Legendre setup failed"):
        reduced._basis(1)

    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_legendre_basis",
        lambda *args: (np.zeros((1, 1), dtype=np.float32), 7),
    )
    with pytest.raises(SpheremackError, match="error code 7"):
        reduced._basis(1)

    cached = np.ones((1, 1), dtype=np.float32, order="F")
    reduced._basis_cache[1] = cached
    assert reduced._basis(1) is cached

    dbasis = np.ones((1, 1), dtype=np.float32, order="F")
    reduced._dbasis_cache[1] = dbasis
    assert reduced._dbasis(1) is dbasis

    reduced._dbasis_cache.clear()
    reduced_module._GLOBAL_DBASIS_CACHE.pop(key, None)
    monkeypatch.setattr(reduced, "_basis", lambda ntrunc: cached)
    monkeypatch.setattr(
        reduced_module,
        "_HAS_REDUCED_GAUSSIAN_LEGENDRE_DERIVATIVE_FROM_BASIS",
        False,
    )
    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_legendre_derivative_basis",
        lambda *args: (dbasis, 0),
    )
    assert reduced._dbasis(1).shape == dbasis.shape

    reduced._dbasis_cache.clear()
    reduced_module._GLOBAL_DBASIS_CACHE.pop(key, None)
    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_legendre_derivative_basis",
        lambda *args: (_ for _ in ()).throw(RuntimeError("dbasis boom")),
    )
    with pytest.raises(SpheremackError, match="derivative setup failed"):
        reduced._dbasis(1)

    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_legendre_derivative_basis",
        lambda *args: (dbasis, 8),
    )
    with pytest.raises(SpheremackError, match="error code 8"):
        reduced._dbasis(1)


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_native_wrapper_error_branches(monkeypatch):
    import skyborn.spharm.reduced_gaussian as reduced_module

    pl = np.array([8, 10, 8], dtype=np.int32)
    reduced = ReducedGaussianSpharmt(pl)
    ntrunc = 1
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    grid = np.zeros((reduced.npoints, 2), dtype=np.float32, order="F")
    spec = np.zeros((nspec, 2), dtype=np.complex64, order="F")
    out_grid = np.zeros_like(grid)
    out_spec = np.zeros_like(spec)
    pair_out = (out_grid, out_grid.copy(order="F"))
    grad_pair_out = (
        out_grid,
        out_grid.copy(order="F"),
        out_grid.copy(order="F"),
        out_grid.copy(order="F"),
    )

    monkeypatch.setattr(reduced, "_basis", lambda n: np.ones((1, 1), dtype=np.float32))
    monkeypatch.setattr(
        reduced, "_basis_transposed", lambda n: np.ones((1, 1), dtype=np.float32)
    )
    monkeypatch.setattr(
        reduced, "_dbasis_transposed", lambda n: np.ones((1, 1), dtype=np.float32)
    )

    cases = [
        (
            "reduced_gaussian_grdtospec",
            lambda: reduced._analyze_scalar(grid, ntrunc),
            "grdtospec",
            (out_spec,),
        ),
        (
            "reduced_gaussian_spectogrd",
            lambda: reduced._synthesize_scalar(spec, ntrunc),
            "spectogrd",
            (out_grid,),
        ),
        (
            "reduced_gaussian_spectogrd_pair",
            lambda: reduced._synthesize_scalar_pair(spec, spec, ntrunc),
            "paired spectogrd",
            pair_out,
        ),
        (
            "reduced_gaussian_getvrtdivspec",
            lambda: reduced._analyze_wind(grid, grid, ntrunc),
            "vector analysis",
            (out_spec, out_spec.copy(order="F")),
        ),
        (
            "reduced_gaussian_getvrtspec",
            lambda: reduced._analyze_vorticity(grid, grid, ntrunc),
            "vorticity analysis",
            (out_spec,),
        ),
        (
            "reduced_gaussian_getdivspec",
            lambda: reduced._analyze_divergence(grid, grid, ntrunc),
            "divergence analysis",
            (out_spec,),
        ),
        (
            "reduced_gaussian_getgrad",
            lambda: reduced._synthesize_gradient(spec, ntrunc),
            "gradient synthesis",
            pair_out,
        ),
        (
            "reduced_gaussian_getgrad_pair",
            lambda: reduced._synthesize_gradient_pair(spec, spec, ntrunc),
            "paired gradient synthesis",
            grad_pair_out,
        ),
    ]

    for func_name, call, message, outputs in cases:
        monkeypatch.setattr(
            reduced_module._spherepack,
            func_name,
            lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        with pytest.raises(SpheremackError, match=message):
            call()

        monkeypatch.setattr(
            reduced_module._spherepack,
            func_name,
            lambda *args, _outputs=outputs, **kwargs: (*_outputs, 9),
        )
        with pytest.raises(SpheremackError, match="error code 9"):
            call()

    monkeypatch.setattr(reduced_module, "_HAS_REDUCED_GAUSSIAN_SPECTOGRD_PAIR", False)
    scalar_calls = {"count": 0}

    def fake_scalar(normalized, ntrunc_arg):
        scalar_calls["count"] += 1
        return np.zeros((reduced.npoints, normalized.shape[1]), dtype=np.float32)

    monkeypatch.setattr(reduced, "_synthesize_scalar", fake_scalar)
    reduced._synthesize_scalar_pair(spec, spec, ntrunc)
    assert scalar_calls["count"] == 2

    monkeypatch.setattr(reduced_module, "_HAS_REDUCED_GAUSSIAN_GETGRAD", False)
    monkeypatch.setattr(
        reduced_module._spherepack,
        "reduced_gaussian_getgrad_pair",
        lambda *args: (*grad_pair_out, 0),
    )
    ugrad, vgrad = reduced._synthesize_gradient(spec, ntrunc)
    assert ugrad.shape == out_grid.shape
    assert vgrad.shape == out_grid.shape


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
def test_reduced_gaussian_invlapspec_and_specsmooth_branches():
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    reduced = ReducedGaussianSpharmt(pl)
    rng = np.random.default_rng(20260520)
    field = rng.standard_normal(reduced.npoints).astype(np.float32)
    spec = reduced.grdtospec(field, ntrunc=5)

    invlap = reduced._invlapspec(spec)
    assert invlap.shape == spec.shape

    with pytest.raises(ValueError, match="smooth must be rank 1"):
        reduced.specsmooth(field, np.ones((reduced.nlat, 1), dtype=np.float32))

    smoothed = reduced.specsmooth(field, np.ones(reduced.nlat, dtype=np.float32))
    assert smoothed.shape == field.shape
