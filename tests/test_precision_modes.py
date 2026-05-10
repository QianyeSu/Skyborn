import sys
from pathlib import Path

import numpy as np
import pytest

src_path = str(Path(__file__).parent.parent / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

try:
    from skyborn.spharm import ReducedGaussianSpharmt, Spharmt
    from skyborn.windspharm import ReducedVectorWind, VectorWind

    AVAILABLE = True
except Exception:
    AVAILABLE = False


@pytest.mark.skipif(not AVAILABLE, reason="skyborn spharm/windspharm unavailable")
def test_spharm_precision_controls_public_dtype():
    sht_auto = Spharmt(24, 9, gridtype="gaussian", precision="auto")
    sht_single = Spharmt(24, 9, gridtype="gaussian", precision="single")
    sht_double = Spharmt(24, 9, gridtype="gaussian", precision="double")

    grid32 = np.random.default_rng(0).standard_normal((9, 24)).astype(np.float32)
    grid64 = grid32.astype(np.float64)

    assert sht_auto.grdtospec(grid32).dtype == np.complex64
    assert sht_auto.grdtospec(grid64).dtype == np.complex128
    assert sht_single.grdtospec(grid64).dtype == np.complex64
    assert sht_double.grdtospec(grid32).dtype == np.complex128

    spec32 = sht_auto.grdtospec(grid32)
    assert sht_auto.spectogrd(spec32).dtype == np.float32
    assert sht_single.spectogrd(spec32).dtype == np.float32
    assert sht_double.spectogrd(spec32).dtype == np.float64


@pytest.mark.skipif(not AVAILABLE, reason="skyborn spharm/windspharm unavailable")
def test_reduced_precision_controls_public_dtype():
    pl = np.full(9, 24, dtype=np.int32)
    rg_auto = ReducedGaussianSpharmt(pl, precision="auto")
    rg_single = ReducedGaussianSpharmt(pl, precision="single")
    rg_double = ReducedGaussianSpharmt(pl, precision="double")

    grid32 = (
        np.random.default_rng(1).standard_normal((int(pl.sum()),)).astype(np.float32)
    )
    grid64 = grid32.astype(np.float64)
    spec32 = np.random.default_rng(2).standard_normal((45,)).astype(np.complex64)
    spec64 = spec32.astype(np.complex128)

    assert rg_auto.grdtospec(grid32).dtype == np.complex64
    assert rg_auto.grdtospec(grid64).dtype == np.complex128
    assert rg_single.grdtospec(grid64).dtype == np.complex64
    assert rg_double.grdtospec(grid32).dtype == np.complex128

    assert rg_auto.spectogrd(spec32).dtype == np.float32
    assert rg_single.spectogrd(spec64).dtype == np.float32
    assert rg_double.spectogrd(spec32).dtype == np.float64


@pytest.mark.skipif(not AVAILABLE, reason="skyborn spharm/windspharm unavailable")
def test_windspharm_precision_is_forwarded():
    u = np.random.default_rng(3).standard_normal((9, 24)).astype(np.float32)
    v = np.random.default_rng(4).standard_normal((9, 24)).astype(np.float32)

    vw_double = VectorWind(u, v, gridtype="gaussian", precision="double")
    vw_single = VectorWind(u, v, gridtype="gaussian", precision="single")

    assert vw_double.vorticity().dtype == np.float64
    assert vw_double.divergence().dtype == np.float64
    assert vw_single.vorticity().dtype == np.float32
    assert vw_single.divergence().dtype == np.float32


@pytest.mark.skipif(not AVAILABLE, reason="skyborn spharm/windspharm unavailable")
def test_windspharm_double_overrides_explicit_restore_dtype():
    u = np.random.default_rng(5).standard_normal((9, 24)).astype(np.float32)
    v = np.random.default_rng(6).standard_normal((9, 24)).astype(np.float32)
    vw = VectorWind(u, v, gridtype="gaussian", precision="double")
    chi = np.random.default_rng(7).standard_normal((9, 24)).astype(np.float32)
    ugrad, vgrad = vw.gradient(chi)
    assert ugrad.dtype == np.float64
    assert vgrad.dtype == np.float64
