"""Tests for packed reduced-Gaussian windspharm interface."""

import sys
from pathlib import Path

import numpy as np
import pytest

src_path = str(Path(__file__).parent.parent / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

try:
    from skyborn.windspharm import ReducedVectorWind, VectorWind

    WINDSPHARM_AVAILABLE = True
except Exception:
    ReducedVectorWind = None
    VectorWind = None
    WINDSPHARM_AVAILABLE = False


@pytest.mark.skipif(not WINDSPHARM_AVAILABLE, reason="windspharm module not available")
def test_reduced_vectorwind_matches_full_gaussian_vectorwind():
    nlat, nlon = 9, 24
    pl = np.full(nlat, nlon, dtype=np.int32)
    rng = np.random.default_rng(20260510)
    u = rng.standard_normal((nlat, nlon, 2)).astype(np.float32)
    v = rng.standard_normal((nlat, nlon, 2)).astype(np.float32)
    scalar = rng.standard_normal((nlat, nlon, 2)).astype(np.float32)

    regular = VectorWind(u, v, gridtype="gaussian")
    reduced = ReducedVectorWind(
        u.reshape(nlat * nlon, 2),
        v.reshape(nlat * nlon, 2),
        pl,
    )

    pairs = [
        (regular.vrtdiv(truncation=6), reduced.vrtdiv(truncation=6)),
        (regular.sfvp(truncation=6), reduced.sfvp(truncation=6)),
        (regular.helmholtz(truncation=6), reduced.helmholtz(truncation=6)),
        (
            regular.gradient(scalar, truncation=6),
            reduced.gradient(scalar.reshape(nlat * nlon, 2), truncation=6),
        ),
    ]

    for regular_outputs, reduced_outputs in pairs:
        for regular_output, reduced_output in zip(regular_outputs, reduced_outputs):
            np.testing.assert_allclose(
                reduced_output,
                regular_output.reshape(nlat * nlon, 2),
                rtol=2e-4,
                atol=5e-1,
            )

    np.testing.assert_allclose(
        reduced.truncate(scalar.reshape(nlat * nlon, 2), truncation=6),
        regular.truncate(scalar, truncation=6).reshape(nlat * nlon, 2),
        rtol=2e-4,
        atol=5e-1,
    )


@pytest.mark.skipif(not WINDSPHARM_AVAILABLE, reason="windspharm module not available")
def test_reduced_vectorwind_precision_aliases():
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    rng = np.random.default_rng(20260512)
    u = rng.standard_normal((int(pl.sum()), 2)).astype(np.float64)
    v = rng.standard_normal((int(pl.sum()), 2)).astype(np.float64)

    reduced = ReducedVectorWind(u, v, pl, precision="float64")
    assert reduced.s.precision == "double"
    assert reduced.magnitude().dtype == np.float64


@pytest.mark.skipif(not WINDSPHARM_AVAILABLE, reason="windspharm module not available")
def test_reduced_vectorwind_varying_pl_smoke():
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    rng = np.random.default_rng(20260510)
    u = rng.standard_normal((int(pl.sum()), 2)).astype(np.float32)
    v = rng.standard_normal((int(pl.sum()), 2)).astype(np.float32)

    reduced = ReducedVectorWind(u, v, pl)

    outputs = []
    outputs.append(reduced.magnitude())
    outputs.extend(reduced.vrtdiv(truncation=5))
    outputs.append(reduced.vorticity(truncation=5))
    outputs.append(reduced.divergence(truncation=5))
    outputs.append(reduced.planetaryvorticity())
    outputs.append(reduced.absolutevorticity(truncation=5))
    outputs.extend(reduced.sfvp(truncation=5))
    outputs.extend(reduced.helmholtz(truncation=5))
    outputs.extend(reduced.gradient(outputs[2], truncation=5))
    outputs.append(reduced.truncate(outputs[2], truncation=5))
    outputs.append(reduced.rossbywavesource(truncation=5))

    for output in outputs:
        assert output.shape == u.shape
        assert np.isfinite(output).all()
        assert output.dtype == np.float32


@pytest.mark.skipif(not WINDSPHARM_AVAILABLE, reason="windspharm module not available")
def test_reduced_vrtdiv_batches_scalar_synthesis(monkeypatch):
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    rng = np.random.default_rng(20260513)
    u = rng.standard_normal((int(pl.sum()), 2)).astype(np.float32)
    v = rng.standard_normal((int(pl.sum()), 2)).astype(np.float32)
    reduced = ReducedVectorWind(u, v, pl)

    monkeypatch.setattr(
        reduced.s,
        "spectogrd",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("vrtdiv should use paired scalar synthesis")
        ),
    )

    vrt, div = reduced.vrtdiv(truncation=5)

    assert vrt.shape == u.shape
    assert div.shape == u.shape
    assert np.isfinite(vrt).all()
    assert np.isfinite(div).all()


@pytest.mark.skipif(not WINDSPHARM_AVAILABLE, reason="windspharm module not available")
def test_reduced_vectorwind_reuses_cached_spectra(monkeypatch):
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    rng = np.random.default_rng(20260514)
    u = rng.standard_normal(int(pl.sum())).astype(np.float32)
    v = rng.standard_normal(int(pl.sum())).astype(np.float32)
    reduced = ReducedVectorWind(u, v, pl)

    calls = {"vrtdiv": 0}
    original_getvrtdivspec = reduced.s.getvrtdivspec

    def counting_getvrtdivspec(*args, **kwargs):
        calls["vrtdiv"] += 1
        return original_getvrtdivspec(*args, **kwargs)

    monkeypatch.setattr(reduced.s, "getvrtdivspec", counting_getvrtdivspec)
    monkeypatch.setattr(
        reduced.s,
        "getvrtspec",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("cached vorticity spectrum should be reused")
        ),
    )
    monkeypatch.setattr(
        reduced.s,
        "getdivspec",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("cached divergence spectrum should be reused")
        ),
    )

    vrt, div = reduced.vrtdiv(truncation=5)
    psi = reduced.streamfunction(truncation=5)
    chi = reduced.velocitypotential(truncation=5)

    assert vrt.shape == u.shape
    assert div.shape == u.shape
    assert psi.shape == u.shape
    assert chi.shape == u.shape
    assert calls["vrtdiv"] == 1


@pytest.mark.skipif(not WINDSPHARM_AVAILABLE, reason="windspharm module not available")
def test_reduced_vectorwind_single_component_uses_single_spectrum(monkeypatch):
    nlat = 9
    pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    rng = np.random.default_rng(20260515)
    u = rng.standard_normal(int(pl.sum())).astype(np.float32)
    v = rng.standard_normal(int(pl.sum())).astype(np.float32)

    reduced = ReducedVectorWind(u, v, pl)
    vrt_calls = {"count": 0}
    original_getvrtspec = reduced.s.getvrtspec

    def counting_getvrtspec(*args, **kwargs):
        vrt_calls["count"] += 1
        return original_getvrtspec(*args, **kwargs)

    monkeypatch.setattr(reduced.s, "getvrtspec", counting_getvrtspec)
    monkeypatch.setattr(
        reduced.s,
        "getvrtdivspec",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("single-component diagnostics should not compute vrt/div")
        ),
    )
    vrt = reduced.vorticity(truncation=5)
    assert vrt_calls["count"] == 1

    reduced = ReducedVectorWind(u, v, pl)
    div_calls = {"count": 0}
    original_getdivspec = reduced.s.getdivspec

    def counting_getdivspec(*args, **kwargs):
        div_calls["count"] += 1
        return original_getdivspec(*args, **kwargs)

    monkeypatch.setattr(reduced.s, "getdivspec", counting_getdivspec)
    monkeypatch.setattr(
        reduced.s,
        "getvrtdivspec",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("single-component diagnostics should not compute vrt/div")
        ),
    )
    div = reduced.divergence(truncation=5)
    assert div_calls["count"] == 1

    assert vrt.shape == u.shape
    assert div.shape == u.shape
    assert np.isfinite(vrt).all()
    assert np.isfinite(div).all()
