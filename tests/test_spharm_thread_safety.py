"""Regression tests for thread-safety / SIMD-correctness fixes in spharm."""

import sys
from concurrent.futures import ThreadPoolExecutor
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

pytestmark = pytest.mark.skipif(
    not SPHARM_AVAILABLE, reason="spharm module not available"
)


def _grid_and_pl(variant):
    if variant == 0:
        pl = np.array([16, 20, 24, 28, 32, 28, 24, 20, 16], dtype=np.int32)
    else:
        pl = np.array([12, 18, 24, 18, 12], dtype=np.int32)
    rng = np.random.default_rng(100 + variant)
    return pl, rng.standard_normal((int(pl.sum()), 2)).astype(np.float32)


def test_reduced_gaussian_concurrent_calls_with_different_grids_match_serial():
    """Callers must never observe each other's FFT tables or work arrays."""
    cases = [_grid_and_pl(0), _grid_and_pl(1)]
    expected = []
    for pl, field in cases:
        sht = ReducedGaussianSpharmt(pl)
        expected.append(sht.grdtospec(field, ntrunc=3))

    def run(idx):
        pl, field = cases[idx % 2]
        sht = ReducedGaussianSpharmt(pl)
        return idx % 2, sht.grdtospec(field, ntrunc=3)

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(run, range(64)))

    for which, spec in results:
        np.testing.assert_array_equal(spec, expected[which])


def test_reduced_gaussian_vector_analysis_is_repeatable_across_calls():
    """Per-call vector work arrays must not leak state between calls."""
    pl = np.array([16, 20, 24, 20, 16], dtype=np.int32)
    rng = np.random.default_rng(7)
    n = int(pl.sum())
    u = rng.standard_normal((n, 2)).astype(np.float32)
    v = rng.standard_normal((n, 2)).astype(np.float32)
    sht = ReducedGaussianSpharmt(pl)

    vrt1, div1 = sht.getvrtdivspec(u, v, ntrunc=3)
    # Interleave a call with a larger truncation to dirty any shared storage.
    sht.getvrtdivspec(u * 2.0, v * 3.0, ntrunc=4)
    vrt2, div2 = sht.getvrtdivspec(u, v, ntrunc=3)

    np.testing.assert_array_equal(vrt1, vrt2)
    np.testing.assert_array_equal(div1, div2)


@pytest.mark.parametrize("nlat,nlon", [(33, 64), (48, 96)])
def test_gaussian_computed_legendre_matches_stored(nlat, nlon):
    """The ``legfunc='computed'`` vector path (vhagc) must match ``stored``."""
    rng = np.random.default_rng(nlat)
    u = rng.standard_normal((nlat, nlon, 3)).astype(np.float32)
    v = rng.standard_normal((nlat, nlon, 3)).astype(np.float32)

    computed = Spharmt(nlon, nlat, gridtype="gaussian", legfunc="computed")
    stored = Spharmt(nlon, nlat, gridtype="gaussian", legfunc="stored")

    vrt_c, div_c = computed.getvrtdivspec(u, v, ntrunc=15)
    vrt_s, div_s = stored.getvrtdivspec(u, v, ntrunc=15)

    scale = max(float(np.max(np.abs(vrt_s))), float(np.max(np.abs(div_s))))
    np.testing.assert_allclose(vrt_c, vrt_s, rtol=0, atol=2e-4 * scale)
    np.testing.assert_allclose(div_c, div_s, rtol=0, atol=2e-4 * scale)
