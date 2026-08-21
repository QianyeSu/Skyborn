import importlib

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose

import skyborn.calc.cape as cape_package
import skyborn.calc.cape.core as cape_core

# ---------------------------------------------------------------------------
# Reference sounding (metpy 1.7.1 docstring case). Dewpoint derived from RH
# with the same Magnus formula used by the Fortran backend.
# ---------------------------------------------------------------------------
_SOUNDING_P = np.array(
    [
        1008.0,
        1000.0,
        950.0,
        900.0,
        850.0,
        800.0,
        750.0,
        700.0,
        650.0,
        600.0,
        550.0,
        500.0,
        450.0,
        400.0,
        350.0,
        300.0,
        250.0,
        200.0,
        175.0,
        150.0,
        125.0,
        100.0,
        80.0,
        70.0,
        60.0,
        50.0,
        40.0,
        30.0,
        25.0,
        20.0,
    ]
)
_SOUNDING_T = np.array(
    [
        29.3,
        28.1,
        23.5,
        20.9,
        18.4,
        15.9,
        13.1,
        10.1,
        6.7,
        3.1,
        -0.5,
        -4.5,
        -9.0,
        -14.8,
        -21.5,
        -29.7,
        -40.0,
        -52.4,
        -59.2,
        -66.5,
        -74.1,
        -78.5,
        -76.0,
        -71.6,
        -66.7,
        -61.3,
        -56.3,
        -51.7,
        -50.7,
        -47.5,
    ]
)
_SOUNDING_RH = np.array(
    [
        0.85,
        0.65,
        0.36,
        0.39,
        0.82,
        0.72,
        0.75,
        0.86,
        0.65,
        0.22,
        0.52,
        0.66,
        0.64,
        0.20,
        0.05,
        0.75,
        0.76,
        0.45,
        0.25,
        0.48,
        0.76,
        0.88,
        0.56,
        0.88,
        0.39,
        0.67,
        0.15,
        0.04,
        0.94,
        0.35,
    ]
)

# metpy 1.7.1 results (Romps LCL + LSODA integrator). The Skyborn backend
# uses the Bolton LCL + fixed-step RK4, so a 5% relative tolerance applies.
_METPY_CAPE = 4830.746079707783
_METPY_CIN = 0.0
_METPY_LFC = 913.55
_METPY_EL = 112.25


def _sounding_dewpoint():
    es = 6.112 * np.exp(17.67 * _SOUNDING_T / (_SOUNDING_T + 243.5))
    e = _SOUNDING_RH * es
    return 243.5 * np.log(e / 6.112) / (17.67 - np.log(e / 6.112))


def _sample_profile():
    """Compact unstable profile with non-zero CIN (stable low levels)."""
    pressure = np.array(
        [
            1000.0,
            950.0,
            900.0,
            850.0,
            800.0,
            750.0,
            700.0,
            650.0,
            600.0,
            550.0,
            500.0,
            450.0,
            400.0,
            350.0,
            300.0,
            250.0,
            200.0,
            150.0,
            100.0,
        ]
    )
    temperature = np.array(
        [
            30.0,
            29.0,
            28.0,
            27.0,
            18.0,
            15.0,
            12.0,
            9.0,
            6.0,
            3.0,
            0.0,
            -4.0,
            -9.0,
            -15.0,
            -23.0,
            -33.0,
            -46.0,
            -64.0,
            -80.0,
        ]
    )
    dewpoint = np.array(
        [
            24.0,
            23.0,
            22.0,
            21.0,
            14.0,
            11.0,
            8.0,
            5.0,
            2.0,
            -1.0,
            -4.0,
            -8.0,
            -13.0,
            -19.0,
            -27.0,
            -37.0,
            -50.0,
            -68.0,
            -85.0,
        ]
    )
    return pressure, temperature, dewpoint


def _sample_grid():
    pressure, temperature, dewpoint = _sample_profile()
    pressure_3d = np.broadcast_to(pressure[:, None, None], (19, 2, 3)).copy()
    temperature_3d = np.broadcast_to(temperature[:, None, None], (19, 2, 3)).copy()
    dewpoint_3d = np.broadcast_to(dewpoint[:, None, None], (19, 2, 3)).copy()
    return pressure_3d, temperature_3d, dewpoint_3d


# ---------------------------------------------------------------------------
# Package exports
# ---------------------------------------------------------------------------
def test_calc_package_exports_cape_submodule_and_function():
    calc_module = importlib.reload(importlib.import_module("skyborn.calc"))

    assert calc_module.cape is cape_package
    assert calc_module.calculate_cape_cin is cape_package.calculate_cape_cin
    assert calc_module.calculate_parcel_profile is cape_package.calculate_parcel_profile
    assert (
        calc_module.calculate_most_unstable_parcel
        is cape_package.calculate_most_unstable_parcel
    )
    assert (
        calc_module.calculate_most_unstable_cape_cin
        is cape_package.calculate_most_unstable_cape_cin
    )


# ---------------------------------------------------------------------------
# Most unstable parcel / CAPE
# ---------------------------------------------------------------------------
def test_most_unstable_parcel_matches_metpy():
    """metpy docstring: depth=50 hPa -> (1008, 29.3, 26.4767, 0)."""
    pressure = _SOUNDING_P
    temperature = _SOUNDING_T
    dewpoint = _sounding_dewpoint()

    mup_p, mup_t, mup_td, mup_idx = cape_package.calculate_most_unstable_parcel(
        pressure, temperature, dewpoint, depth=50.0
    )
    assert_allclose(mup_p, 1008.0, rtol=1e-6)
    assert_allclose(mup_t, 29.3, rtol=1e-6)
    # dewpoint differs slightly (Magnus vs metpy formula), ~0.05 K tolerance
    assert_allclose(mup_td, 26.4767, atol=0.1)
    assert mup_idx == 0


def test_most_unstable_cape_cin_matches_metpy():
    """metpy: most_unstable_cape_cin(p, T, Td) = (4830.75, 0) for this sounding."""
    pressure = _SOUNDING_P
    temperature = _SOUNDING_T
    dewpoint = _sounding_dewpoint()

    mucape, mucin, _, _ = cape_package.calculate_most_unstable_cape_cin(
        pressure, temperature, dewpoint, depth=300.0
    )
    assert_allclose(mucape, _METPY_CAPE, rtol=0.05, atol=50.0)
    assert_allclose(mucin, 0.0, atol=1.0)


def test_most_unstable_grid_matches_profile():
    """Grid results must equal the per-column 1D results."""
    pressure, temperature, dewpoint = _sample_profile()
    p3 = np.broadcast_to(pressure[:, None, None], (19, 2, 3)).copy()
    t3 = np.broadcast_to(temperature[:, None, None], (19, 2, 3)).copy()
    td3 = np.broadcast_to(dewpoint[:, None, None], (19, 2, 3)).copy()

    prof = cape_package.calculate_most_unstable_parcel(
        pressure, temperature, dewpoint, depth=300.0
    )
    grid = cape_package.calculate_most_unstable_parcel(p3, t3, td3, depth=300.0)
    # grid returns (p3, t3, td3, idx3): the index field stores the selected
    # level index at that level's position; other levels stay -1
    assert grid[3][int(prof[3]), 0, 0] == prof[3]
    assert_allclose(grid[0][int(prof[3]), 0, 0], prof[0], rtol=1e-12)

    mucape_prof = cape_package.calculate_most_unstable_cape_cin(
        pressure, temperature, dewpoint, depth=300.0
    )
    mucape_grid = cape_package.calculate_most_unstable_cape_cin(
        p3, t3, td3, depth=300.0
    )
    assert_allclose(mucape_grid[0][0, 0], mucape_prof[0], rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Public API vs backend consistency
# ---------------------------------------------------------------------------
def test_cape_public_api_matches_backend_profile():
    pressure, temperature, dewpoint = _sample_profile()

    public_result = cape_package.calculate_cape_cin(pressure, temperature, dewpoint)
    backend_result = cape_core.cape_profile(pressure, temperature, dewpoint)

    assert len(public_result) == 4
    assert_allclose(public_result, backend_result, rtol=1e-12, atol=1e-12)


def test_cape_public_api_matches_backend_grid_and_xarray():
    pressure_3d, temperature_3d, dewpoint_3d = _sample_grid()

    grid_result = cape_package.calculate_cape_cin(
        pressure_3d, temperature_3d, dewpoint_3d
    )
    backend_grid = cape_core.cape_grid(pressure_3d, temperature_3d, dewpoint_3d)
    for public, backend in zip(grid_result, backend_grid):
        assert_allclose(public, backend, rtol=1e-12, atol=1e-12)

    pressure_xr = xr.DataArray(
        pressure_3d,
        dims=["level", "lat", "lon"],
        coords={
            "level": pressure_3d[:, 0, 0],
            "lat": [0.0, 10.0],
            "lon": [100.0, 110.0, 120.0],
        },
        attrs={"units": "hPa"},
    )
    temperature_xr = xr.DataArray(temperature_3d, dims=["level", "lat", "lon"])
    dewpoint_xr = xr.DataArray(dewpoint_3d, dims=["level", "lat", "lon"])

    xarray_result = cape_package.calculate_cape_cin(
        pressure_xr, temperature_xr, dewpoint_xr
    )
    assert all(isinstance(r, xr.DataArray) for r in xarray_result)
    assert xarray_result[0].dims == ("lat", "lon")
    assert xarray_result[0].coords["lat"].values.tolist() == [0.0, 10.0]
    for i, backend in enumerate(backend_grid):
        assert_allclose(xarray_result[i].values, backend, rtol=1e-12, atol=1e-12)


def test_cape_backend_is_consistent_between_profile_and_grid():
    pressure, temperature, dewpoint = _sample_profile()
    pressure_3d, temperature_3d, dewpoint_3d = _sample_grid()

    profile_result = cape_core.cape_profile(pressure, temperature, dewpoint)
    grid_result = cape_core.cape_grid(pressure_3d, temperature_3d, dewpoint_3d)

    for i, value in enumerate(profile_result):
        assert_allclose(
            grid_result[i],
            np.full((2, 3), value, dtype=np.float64),
            rtol=1e-12,
            atol=1e-12,
        )


# ---------------------------------------------------------------------------
# Correctness against the metpy 1.7.1 reference sounding
# ---------------------------------------------------------------------------
def test_cape_matches_metpy_reference():
    pressure = _SOUNDING_P
    temperature = _SOUNDING_T
    dewpoint = _sounding_dewpoint()

    cape, cin, lfc, el = cape_package.calculate_cape_cin(
        pressure, temperature, dewpoint
    )

    # Bolton LCL + RK4 vs Romps LCL + LSODA: expect agreement within 5%
    assert_allclose(cape, _METPY_CAPE, rtol=0.05, atol=50.0)
    assert_allclose(cin, _METPY_CIN, atol=1.0)
    assert_allclose(lfc, _METPY_LFC, rtol=0.02)
    assert_allclose(el, _METPY_EL, rtol=0.02)


def test_cape_cin_sign_convention():
    """CIN is reported <= 0 (metpy convention); positive values are clipped."""
    pressure, temperature, dewpoint = _sample_profile()
    cape, cin, lfc, el = cape_package.calculate_cape_cin(
        pressure, temperature, dewpoint
    )
    assert cape > 0.0
    assert cin <= 0.0
    assert 0.0 < lfc < pressure[0]
    assert el < lfc


def test_cape_handles_nan_levels():
    pressure, temperature, dewpoint = _sample_profile()
    p2 = pressure.copy()
    t2 = temperature.copy()
    td2 = dewpoint.copy()
    p2[3] = np.nan
    t2[3] = np.nan
    td2[3] = np.nan

    cape, cin, _, _ = cape_package.calculate_cape_cin(p2, t2, td2)
    assert np.isfinite(cape)
    assert np.isfinite(cin)


def test_cape_stable_profile_returns_zero():
    """A fully stable profile has no LFC: CAPE and CIN are zero."""
    pressure = np.array([1000.0, 900.0, 800.0, 700.0, 600.0, 500.0])
    temperature = np.array([20.0, 18.0, 16.0, 14.0, 12.0, 10.0])
    dewpoint = np.array([5.0, 3.0, 1.0, -1.0, -3.0, -5.0])

    cape, cin, lfc, el = cape_package.calculate_cape_cin(
        pressure, temperature, dewpoint
    )
    assert cape == 0.0
    assert cin == 0.0
    assert np.isnan(lfc)
    assert np.isnan(el)


# ---------------------------------------------------------------------------
# Parcel profile
# ---------------------------------------------------------------------------
def test_parcel_profile_matches_metpy_reference():
    pressure = _SOUNDING_P
    temperature = _SOUNDING_T
    dewpoint = _sounding_dewpoint()

    profile = cape_package.calculate_parcel_profile(pressure, temperature, dewpoint)

    # metpy 1.7.1 parcel profile at the first levels (degC)
    metpy_first = np.array(
        [
            29.3,
            28.61221951574845,
            25.174081108824907,
            23.410446413478496,
            21.5304966888134,
        ]
    )
    assert_allclose(profile[:5], metpy_first, rtol=0.02, atol=0.5)


def test_parcel_profile_grid_matches_backend():
    pressure_3d, temperature_3d, dewpoint_3d = _sample_grid()

    public = cape_package.calculate_parcel_profile(
        pressure_3d, temperature_3d, dewpoint_3d
    )
    backend = cape_core.parcel_profile_grid(pressure_3d, temperature_3d, dewpoint_3d)
    assert_allclose(public, backend, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
def test_cape_input_validation():
    with pytest.raises(ValueError, match="pressure must be 1D or 3D"):
        cape_package.calculate_cape_cin(
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
        )

    with pytest.raises(ValueError, match="pressure must be 1D"):
        cape_core.cape_profile(
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
        )

    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.cape_profile(
            np.ones(5, dtype=float),
            np.ones(6, dtype=float),
            np.ones(5, dtype=float),
        )


def test_cape_grid_validation_3d():
    """Test validation for cape_grid: must be 3D and shapes must match."""
    # Test non-3D input
    with pytest.raises(ValueError, match="pressure must be 3D"):
        cape_core.cape_grid(
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.cape_grid(
            np.ones((5, 3, 4), dtype=float),
            np.ones((5, 3, 5), dtype=float),
            np.ones((5, 3, 4), dtype=float),
        )


def test_parcel_profile_validation():
    """Test validation for parcel_profile: must be 1D and shapes must match."""
    # Test non-1D input
    with pytest.raises(ValueError, match="pressure must be 1D"):
        cape_core.parcel_profile(
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.parcel_profile(
            np.ones(5, dtype=float),
            np.ones(6, dtype=float),
            np.ones(5, dtype=float),
        )


def test_parcel_profile_grid_validation():
    """Test validation for parcel_profile_grid: must be 3D and shapes must match."""
    # Test non-3D input
    with pytest.raises(ValueError, match="pressure must be 3D"):
        cape_core.parcel_profile_grid(
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.parcel_profile_grid(
            np.ones((5, 3, 4), dtype=float),
            np.ones((5, 3, 5), dtype=float),
            np.ones((5, 3, 4), dtype=float),
        )


def test_calculate_parcel_profile_invalid_dimension():
    """Test that calculate_parcel_profile rejects 2D input."""
    with pytest.raises(ValueError, match="pressure must be 1D or 3D"):
        cape_package.calculate_parcel_profile(
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
        )


def test_most_unstable_parcel_validation():
    """Test validation for most_unstable_parcel: must be 1D and shapes must match."""
    # Test non-1D input
    with pytest.raises(ValueError, match="pressure must be 1D"):
        cape_core.most_unstable_parcel(
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.most_unstable_parcel(
            np.ones(5, dtype=float),
            np.ones(6, dtype=float),
            np.ones(5, dtype=float),
        )


def test_most_unstable_parcel_grid_validation():
    """Test validation for most_unstable_parcel_grid: must be 3D and shapes must match."""
    # Test non-3D input
    with pytest.raises(ValueError, match="pressure must be 3D"):
        cape_core.most_unstable_parcel_grid(
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.most_unstable_parcel_grid(
            np.ones((5, 3, 4), dtype=float),
            np.ones((5, 3, 5), dtype=float),
            np.ones((5, 3, 4), dtype=float),
        )


def test_most_unstable_cape_cin_validation():
    """Test validation for most_unstable_cape_cin: must be 1D and shapes must match."""
    # Test non-1D input
    with pytest.raises(ValueError, match="pressure must be 1D"):
        cape_core.most_unstable_cape_cin(
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.most_unstable_cape_cin(
            np.ones(5, dtype=float),
            np.ones(6, dtype=float),
            np.ones(5, dtype=float),
        )


def test_most_unstable_cape_cin_grid_validation():
    """Test validation for most_unstable_cape_cin_grid: must be 3D and shapes must match."""
    # Test non-3D input
    with pytest.raises(ValueError, match="pressure must be 3D"):
        cape_core.most_unstable_cape_cin_grid(
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
            np.ones((5,), dtype=float),
        )

    # Test shape mismatch
    with pytest.raises(ValueError, match="must share a shape"):
        cape_core.most_unstable_cape_cin_grid(
            np.ones((5, 3, 4), dtype=float),
            np.ones((5, 3, 5), dtype=float),
            np.ones((5, 3, 4), dtype=float),
        )


def test_calculate_most_unstable_parcel_invalid_dimension():
    """Test that calculate_most_unstable_parcel rejects 2D input."""
    with pytest.raises(ValueError, match="pressure must be 1D or 3D"):
        cape_package.calculate_most_unstable_parcel(
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
        )


def test_calculate_most_unstable_cape_cin_invalid_dimension():
    """Test that calculate_most_unstable_cape_cin rejects 2D input."""
    with pytest.raises(ValueError, match="pressure must be 1D or 3D"):
        cape_package.calculate_most_unstable_cape_cin(
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
            np.ones((5, 3), dtype=float),
        )


def test_xarray_3d_wrapping_with_missing_coords():
    """Test _wrap_xarray_3d when some dimensions have no coordinates."""
    pressure_3d, temperature_3d, dewpoint_3d = _sample_grid()

    # Create xarray without coordinates for one dimension
    p_xr = xr.DataArray(
        pressure_3d,
        dims=["level", "lat", "lon"],
        coords={
            "level": np.arange(pressure_3d.shape[0]),
            # 'lat' intentionally missing
            "lon": np.arange(pressure_3d.shape[2]),
        },
    )
    t_xr = xr.DataArray(temperature_3d, dims=["level", "lat", "lon"])
    td_xr = xr.DataArray(dewpoint_3d, dims=["level", "lat", "lon"])

    # This should work without error - missing coords are just not copied
    profile = cape_package.calculate_parcel_profile(p_xr, t_xr, td_xr)
    assert isinstance(profile, xr.DataArray)
    assert "level" in profile.coords
    assert "lon" in profile.coords
