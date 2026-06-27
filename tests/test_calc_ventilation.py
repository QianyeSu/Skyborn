"""Tests for skyborn.calc.ventilation module.

Tests the ventilated potential intensity (vPI) framework and the
ventilated genesis potential index (GPIv) following Chavas et al. (2025).
"""

import numpy as np
import pytest
import xarray as xr

from skyborn.calc.ventilation import (
    absolute_vorticity_850,
    entropy_deficit,
    genesis_potential_index,
    ventilated_pi,
    ventilation_index,
    vertical_wind_shear,
)
from skyborn.calc.ventilation.ventilation import VI_MAX

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_atmospheric_dataset():
    """Create a small synthetic ERA5-like dataset for testing."""
    np.random.seed(42)

    levels = [200, 600, 850, 1000]
    lats = np.arange(10, 41, 2.5)
    lons = np.arange(100, 181, 2.5)
    times = np.arange(4)

    shape = (len(times), len(levels), len(lats), len(lons))

    # Realistic-ish temperature profile (K) — reshape for broadcasting
    T_base = np.array([220, 275, 285, 290]).reshape(1, -1, 1, 1)
    T = np.broadcast_to(T_base, shape).astype(float).copy()
    T += np.random.randn(*shape) * 2.0

    # Specific humidity (kg/kg)
    Q_base = np.array([1e-5, 0.003, 0.008, 0.012]).reshape(1, -1, 1, 1)
    Q = np.broadcast_to(Q_base, shape).astype(float).copy()
    Q += np.random.randn(*shape) * 0.0005

    # Relative humidity (%)
    R_base = np.array([10, 40, 60, 80]).reshape(1, -1, 1, 1)
    R = np.broadcast_to(R_base, shape).astype(float).copy()
    R += np.random.randn(*shape) * 3.0

    # Zonal wind (m/s)
    U_base = np.array([5, 2, -5, -3]).reshape(1, -1, 1, 1)
    U = np.broadcast_to(U_base, shape).astype(float).copy()
    U += np.random.randn(*shape) * 1.0

    # Meridional wind (m/s)
    V_base = np.array([2, 1, -2, -1]).reshape(1, -1, 1, 1)
    V = np.broadcast_to(V_base, shape).astype(float).copy()
    V += np.random.randn(*shape) * 1.0

    ds = xr.Dataset(
        data_vars={
            "T": (["time", "level", "latitude", "longitude"], T),
            "Q": (["time", "level", "latitude", "longitude"], Q),
            "R": (["time", "level", "latitude", "longitude"], R),
            "U": (["time", "level", "latitude", "longitude"], U),
            "V": (["time", "level", "latitude", "longitude"], V),
        },
        coords={
            "time": times,
            "level": levels,
            "latitude": lats,
            "longitude": lons,
        },
    )
    return ds


@pytest.fixture
def sample_pi_field(sample_atmospheric_dataset):
    """Create a synthetic PI field matching the dataset grid."""
    ds = sample_atmospheric_dataset
    shape = (ds.sizes["time"], ds.sizes["latitude"], ds.sizes["longitude"])
    pi = xr.DataArray(
        np.random.uniform(40, 80, shape),
        dims=["time", "latitude", "longitude"],
        coords={
            "time": ds.time,
            "latitude": ds.latitude,
            "longitude": ds.longitude,
        },
        attrs={"units": "m/s", "long_name": "Potential Intensity"},
    )
    return pi


# ---------------------------------------------------------------------------
# vertical_wind_shear
# ---------------------------------------------------------------------------


class TestVerticalWindShear:
    """Tests for vertical_wind_shear."""

    def test_returns_dataarray(self, sample_atmospheric_dataset):
        vws = vertical_wind_shear(sample_atmospheric_dataset)
        assert isinstance(vws, xr.DataArray)

    def test_shape_no_level_dim(self, sample_atmospheric_dataset):
        ds = sample_atmospheric_dataset
        vws = vertical_wind_shear(ds)
        assert "level" not in vws.dims
        assert vws.sizes["time"] == ds.sizes["time"]
        assert vws.sizes["latitude"] == ds.sizes["latitude"]
        assert vws.sizes["longitude"] == ds.sizes["longitude"]

    def test_values_nonnegative(self, sample_atmospheric_dataset):
        vws = vertical_wind_shear(sample_atmospheric_dataset)
        valid = vws.values[~np.isnan(vws.values)]
        assert np.all(valid >= 0)

    def test_attrs_set(self, sample_atmospheric_dataset):
        vws = vertical_wind_shear(sample_atmospheric_dataset)
        assert "units" in vws.attrs
        assert vws.attrs["units"] == "m/s"

    def test_known_values(self):
        """Test with known wind values."""
        ds = xr.Dataset(
            data_vars={
                "U": (
                    ["level", "latitude", "longitude"],
                    np.array([[[10, 10], [10, 10]], [[5, 5], [5, 5]]]),
                ),
                "V": (
                    ["level", "latitude", "longitude"],
                    np.array([[[0, 0], [0, 0]], [[0, 0], [0, 0]]]),
                ),
            },
            coords={
                "level": [200, 850],
                "latitude": [0, 10],
                "longitude": [0, 10],
            },
        )
        vws = vertical_wind_shear(ds)
        expected = np.full((2, 2), 5.0)
        np.testing.assert_array_almost_equal(vws.values, expected)


# ---------------------------------------------------------------------------
# entropy_deficit
# ---------------------------------------------------------------------------


class TestEntropyDeficit:
    """Tests for entropy_deficit."""

    def test_returns_dataarray(self, sample_atmospheric_dataset):
        chi = entropy_deficit(sample_atmospheric_dataset)
        assert isinstance(chi, xr.DataArray)

    def test_shape_no_level_dim(self, sample_atmospheric_dataset):
        ds = sample_atmospheric_dataset
        chi = entropy_deficit(ds)
        assert "level" not in chi.dims

    def test_values_in_range(self, sample_atmospheric_dataset):
        chi = entropy_deficit(sample_atmospheric_dataset)
        valid = chi.values[~np.isnan(chi.values)]
        assert np.all(valid >= 0.02)
        assert np.all(valid <= 1.5)

    def test_rh_based_calculation(self):
        """Test that chi = 1 - RH/100 when R is available."""
        ds = xr.Dataset(
            data_vars={
                "R": (
                    ["level", "latitude", "longitude"],
                    np.array([[[50, 60], [70, 80]]]),
                ),
            },
            coords={
                "level": [600],
                "latitude": [0, 10],
                "longitude": [0, 10],
            },
        )
        chi = entropy_deficit(ds)
        expected = np.array([[0.5, 0.4], [0.3, 0.2]])
        np.testing.assert_array_almost_equal(chi.values, expected)

    def test_q_fallback(self):
        """Test specific humidity fallback when R is not available."""
        ds = xr.Dataset(
            data_vars={
                "Q": (
                    ["level", "latitude", "longitude"],
                    np.array([[[0.005, 0.010], [0.015, 0.020]]]),
                ),
            },
            coords={
                "level": [600],
                "latitude": [0, 10],
                "longitude": [0, 10],
            },
        )
        chi = entropy_deficit(ds, rh_var="R_missing")
        valid = chi.values[~np.isnan(chi.values)]
        assert np.all(valid >= 0.02)
        assert np.all(valid <= 1.5)


# ---------------------------------------------------------------------------
# ventilation_index
# ---------------------------------------------------------------------------


class TestVentilationIndex:
    """Tests for ventilation_index."""

    def test_returns_dataarray(self, sample_atmospheric_dataset, sample_pi_field):
        ds = sample_atmospheric_dataset
        vws = vertical_wind_shear(ds)
        chi = entropy_deficit(ds)
        vi = ventilation_index(vws, chi, sample_pi_field)
        assert isinstance(vi, xr.DataArray)

    def test_formula(self):
        """Test VI = VWS * Chi / PI with known values."""
        vws = xr.DataArray(np.array([10.0, 20.0]), dims=["x"])
        chi = xr.DataArray(np.array([0.5, 0.3]), dims=["x"])
        pi = xr.DataArray(np.array([50.0, 80.0]), dims=["x"])
        vi = ventilation_index(vws, chi, pi)
        expected = np.array([0.1, 0.075])
        np.testing.assert_array_almost_equal(vi.values, expected)

    def test_negative_vi_is_nan(self):
        """VI should be NaN where the formula gives non-positive values."""
        vws = xr.DataArray(np.array([0.0, -1.0]), dims=["x"])
        chi = xr.DataArray(np.array([0.5, 0.3]), dims=["x"])
        pi = xr.DataArray(np.array([50.0, 80.0]), dims=["x"])
        vi = ventilation_index(vws, chi, pi)
        assert np.isnan(vi.values[0])
        assert np.isnan(vi.values[1])


# ---------------------------------------------------------------------------
# ventilated_pi
# ---------------------------------------------------------------------------


class TestVentilatedPI:
    """Tests for ventilated_pi — the critical analytic cubic solution."""

    def test_returns_dataarray(self):
        pi = xr.DataArray(np.array([70.0, 80.0]), dims=["x"])
        vi = xr.DataArray(np.array([0.05, 0.10]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        assert isinstance(vpi, xr.DataArray)

    def test_vi_near_zero_vpi_equals_pi(self):
        """At VI → 0, vPI should approach PI (ratio → 1)."""
        pi = xr.DataArray(np.array([70.0]), dims=["x"])
        vi = xr.DataArray(np.array([0.001]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        assert abs(vpi.values[0] - 70.0) < 1.0  # Within 1 m/s

    def test_vi_at_vimax_ratio_is_1_over_sqrt3(self):
        """At VI = VI_MAX, ratio should be 1/sqrt(3) ≈ 0.5774."""
        pi_val = 70.0
        pi = xr.DataArray(np.array([pi_val]), dims=["x"])
        vi = xr.DataArray(np.array([VI_MAX]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        expected_ratio = 1.0 / np.sqrt(3.0)
        np.testing.assert_almost_equal(
            vpi.values[0], pi_val * expected_ratio, decimal=3
        )

    def test_vi_above_vimax_is_nan(self):
        """VI > VI_MAX should produce NaN."""
        pi = xr.DataArray(np.array([70.0]), dims=["x"])
        vi = xr.DataArray(np.array([0.146]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        assert np.isnan(vpi.values[0])

    def test_vi_negative_is_nan(self):
        """Negative VI should produce NaN."""
        pi = xr.DataArray(np.array([70.0]), dims=["x"])
        vi = xr.DataArray(np.array([-0.05]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        assert np.isnan(vpi.values[0])

    def test_vi_nan_is_nan(self):
        """NaN VI should produce NaN."""
        pi = xr.DataArray(np.array([70.0]), dims=["x"])
        vi = xr.DataArray(np.array([np.nan]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        assert np.isnan(vpi.values[0])

    def test_monotonic_decrease(self):
        """vPI should decrease monotonically as VI increases."""
        pi_val = 80.0
        vi_values = np.linspace(0.01, VI_MAX, 20)
        pi = xr.DataArray(np.full(20, pi_val), dims=["x"])
        vi = xr.DataArray(vi_values, dims=["x"])
        vpi = ventilated_pi(pi, vi)
        diffs = np.diff(vpi.values)
        assert np.all(diffs < 0), "vPI should decrease as VI increases"

    def test_attrs_preserved(self):
        pi = xr.DataArray(np.array([70.0]), dims=["x"])
        vi = xr.DataArray(np.array([0.05]), dims=["x"])
        vpi = ventilated_pi(pi, vi)
        assert "units" in vpi.attrs
        assert vpi.attrs["units"] == "m/s"

    def test_custom_vi_max(self):
        """Test with a custom vi_max parameter."""
        pi = xr.DataArray(np.array([70.0]), dims=["x"])
        vi = xr.DataArray(np.array([0.20]), dims=["x"])
        # With default vi_max=0.145, VI=0.20 > 0.145 → NaN
        vpi_default = ventilated_pi(pi, vi)
        assert np.isnan(vpi_default.values[0])
        # With vi_max=0.25, VI=0.20 is valid
        vpi_custom = ventilated_pi(pi, vi, vi_max=0.25)
        assert np.isfinite(vpi_custom.values[0])
        assert vpi_custom.values[0] < 70.0  # Should be reduced


# ---------------------------------------------------------------------------
# absolute_vorticity_850
# ---------------------------------------------------------------------------


class TestAbsoluteVorticity850:
    """Tests for absolute_vorticity_850."""

    def test_returns_dataarray(self, sample_atmospheric_dataset):
        eta_c = absolute_vorticity_850(sample_atmospheric_dataset)
        assert isinstance(eta_c, xr.DataArray)

    def test_shape_no_level_dim(self, sample_atmospheric_dataset):
        ds = sample_atmospheric_dataset
        eta_c = absolute_vorticity_850(ds)
        assert "level" not in eta_c.dims

    def test_values_capped(self, sample_atmospheric_dataset):
        """All non-NaN values should be <= cap."""
        eta_c = absolute_vorticity_850(sample_atmospheric_dataset)
        valid = eta_c.values[~np.isnan(eta_c.values)]
        assert np.all(np.abs(valid) <= 3.7e-5 + 1e-15)

    def test_values_positive(self, sample_atmospheric_dataset):
        """All non-NaN values should be positive."""
        eta_c = absolute_vorticity_850(sample_atmospheric_dataset)
        valid = eta_c.values[~np.isnan(eta_c.values)]
        assert np.all(valid > 0)

    def test_attrs_set(self, sample_atmospheric_dataset):
        eta_c = absolute_vorticity_850(sample_atmospheric_dataset)
        assert "units" in eta_c.attrs
        assert eta_c.attrs["units"] == "s^-1"

    def test_custom_cap(self, sample_atmospheric_dataset):
        """Test with a custom cap value."""
        eta_c = absolute_vorticity_850(sample_atmospheric_dataset, cap=5.0e-5)
        valid = eta_c.values[~np.isnan(eta_c.values)]
        assert np.all(np.abs(valid) <= 5.0e-5 + 1e-15)

    def test_coriolis_at_equator(self):
        """At the equator, Coriolis should be ~0."""
        ds = xr.Dataset(
            data_vars={
                "U": (
                    ["level", "latitude", "longitude"],
                    np.zeros((1, 3, 3)),
                ),
                "V": (
                    ["level", "latitude", "longitude"],
                    np.zeros((1, 3, 3)),
                ),
            },
            coords={
                "level": [850],
                "latitude": [0.0, 5.0, 10.0],
                "longitude": [0.0, 5.0, 10.0],
            },
        )
        eta_c = absolute_vorticity_850(ds)
        # At equator (latitude=0), Coriolis=0 and vorticity≈0 → should be NaN (filtered)
        assert np.isnan(eta_c.sel(latitude=0.0).values).all() or np.all(
            eta_c.sel(latitude=0.0).values >= 0
        )


# ---------------------------------------------------------------------------
# genesis_potential_index
# ---------------------------------------------------------------------------


class TestGenesisPotentialIndex:
    """Tests for genesis_potential_index."""

    def test_returns_dataarray(self):
        vpi = xr.DataArray(
            np.array([[50.0, 60.0], [40.0, 45.0]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [10.0, 20.0], "longitude": [100.0, 110.0]},
        )
        eta_c = xr.DataArray(
            np.array([[3e-5, 2e-5], [1e-5, 3.7e-5]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [10.0, 20.0], "longitude": [100.0, 110.0]},
        )
        gpi_v = genesis_potential_index(vpi, eta_c)
        assert isinstance(gpi_v, xr.DataArray)

    def test_values_positive(self):
        vpi = xr.DataArray(
            np.array([[50.0, 60.0]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [10.0], "longitude": [100.0, 110.0]},
        )
        eta_c = xr.DataArray(
            np.array([[3e-5, 2e-5]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [10.0], "longitude": [100.0, 110.0]},
        )
        gpi_v = genesis_potential_index(vpi, eta_c)
        valid = gpi_v.values[~np.isnan(gpi_v.values)]
        assert np.all(valid > 0)

    def test_custom_grid_spacing(self):
        """Test with custom dx/dy."""
        vpi = xr.DataArray(
            np.array([[50.0]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [15.0], "longitude": [120.0]},
        )
        eta_c = xr.DataArray(
            np.array([[3e-5]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [15.0], "longitude": [120.0]},
        )
        gpi_default = genesis_potential_index(vpi, eta_c)
        gpi_custom = genesis_potential_index(vpi, eta_c, dx=1.0, dy=1.0)
        # Custom dx/dy=1 should be 1/4 of default (2*2)
        ratio = gpi_custom.values[0] / gpi_default.values[0]
        np.testing.assert_almost_equal(ratio, 0.25, decimal=5)

    def test_attrs_set(self):
        vpi = xr.DataArray(
            np.array([[50.0]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [15.0], "longitude": [120.0]},
        )
        eta_c = xr.DataArray(
            np.array([[3e-5]]),
            dims=["latitude", "longitude"],
            coords={"latitude": [15.0], "longitude": [120.0]},
        )
        gpi_v = genesis_potential_index(vpi, eta_c)
        assert "long_name" in gpi_v.attrs


# ---------------------------------------------------------------------------
# Integration test
# ---------------------------------------------------------------------------


class TestVentilationIntegration:
    """Integration tests for the full ventilation pipeline."""

    def test_full_pipeline(self, sample_atmospheric_dataset, sample_pi_field):
        """Test the complete calculation chain: VWS → χ → VI → vPI → η_c → GPIv."""
        ds = sample_atmospheric_dataset
        pi = sample_pi_field

        # Step 1: VWS
        vws = vertical_wind_shear(ds)
        assert "level" not in vws.dims

        # Step 2: Entropy deficit
        chi = entropy_deficit(ds)
        assert "level" not in chi.dims

        # Step 3: Ventilation index
        vi = ventilation_index(vws, chi, pi)
        valid_vi = vi.values[~np.isnan(vi.values)]
        assert np.all(valid_vi > 0)

        # Step 4: Ventilated PI
        vpi = ventilated_pi(pi, vi)
        valid_vpi = vpi.values[~np.isnan(vpi.values)]
        assert np.all(valid_vpi > 0)
        # vPI should be <= PI everywhere
        valid_mask = np.isfinite(vpi.values) & np.isfinite(pi.values)
        assert np.all(vpi.values[valid_mask] <= pi.values[valid_mask] + 1e-10)

        # Step 5: Absolute vorticity
        eta_c = absolute_vorticity_850(ds)
        valid_eta = eta_c.values[~np.isnan(eta_c.values)]
        assert np.all(valid_eta > 0)
        assert np.all(np.abs(valid_eta) <= 3.7e-5 + 1e-15)

        # Step 6: GPIv
        gpi_v = genesis_potential_index(vpi, eta_c)
        valid_gpi = gpi_v.values[~np.isnan(gpi_v.values)]
        assert np.all(valid_gpi > 0)

    def test_vpi_decreases_with_increasing_vi(self):
        """vPI should decrease as VI increases (the core physics)."""
        pi_val = 75.0
        vi_low = xr.DataArray(np.array([0.02]), dims=["x"])
        vi_high = xr.DataArray(np.array([0.12]), dims=["x"])
        pi = xr.DataArray(np.array([pi_val]), dims=["x"])

        vpi_low = ventilated_pi(pi, vi_low)
        vpi_high = ventilated_pi(pi, vi_high)

        assert vpi_high.values[0] < vpi_low.values[0]
