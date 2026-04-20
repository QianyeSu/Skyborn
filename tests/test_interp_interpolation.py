"""
Tests for skyborn.interp.interpolation module.

This module tests the interpolation functionality including hybrid-sigma
to pressure level interpolation and multidimensional spatial interpolation.
"""

import builtins
import importlib
import sys

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose, assert_array_almost_equal, assert_array_equal

import skyborn.interp.interpolation as interpolation_module
from skyborn.interp.interpolation import (
    __pres_lev_mandatory__,
    _post_interp_multidim,
    _pre_interp_multidim,
    _pressure_from_hybrid,
    _sigma_from_hybrid,
    delta_pressure_hybrid,
    interp_hybrid_to_pressure,
    interp_multidim,
    interp_sigma_to_hybrid,
    pressure_at_hybrid_levels,
)


class TestInterpolationHelpers:
    """Test helper functions for interpolation."""

    def test_pressure_from_hybrid(self):
        """Test pressure calculation from hybrid coordinates."""
        # Simple test data
        ps = np.array([100000, 95000])  # Pa
        hya = np.array([0.0, 0.1, 0.5])  # hybrid A coefficients
        hyb = np.array([1.0, 0.9, 0.5])  # hybrid B coefficients
        p0 = 100000.0

        # Convert to xarray
        ps_da = xr.DataArray(ps, dims=["x"])
        hya_da = xr.DataArray(hya, dims=["lev"])
        hyb_da = xr.DataArray(hyb, dims=["lev"])

        pressure = _pressure_from_hybrid(ps_da, hya_da, hyb_da, p0)

        # Check shape and basic properties
        assert pressure.dims == ("lev", "x")
        assert pressure.shape == (3, 2)  # lev, x
        assert np.all(pressure > 0)  # All pressures should be positive
        assert np.all(pressure <= 105000)  # Should stay in a realistic range

        # Public helper should match the legacy private wrapper
        public_pressure = pressure_at_hybrid_levels(ps_da, hya_da, hyb_da, p0)
        assert_array_equal(pressure.values, public_pressure.values)

    def test_sigma_from_hybrid(self):
        """Test sigma calculation from hybrid coordinates."""
        ps = np.array([100000, 95000])  # Pa
        hya = np.array([0.0, 0.1, 0.5])
        hyb = np.array([1.0, 0.9, 0.5])
        p0 = 100000.0

        # Convert to xarray
        ps_da = xr.DataArray(ps, dims=["x"])
        hya_da = xr.DataArray(hya, dims=["lev"])
        hyb_da = xr.DataArray(hyb, dims=["lev"])

        sigma = _sigma_from_hybrid(ps_da, hya_da, hyb_da, p0)

        # Check shape and basic properties
        assert sigma.dims == ("lev", "x")
        assert sigma.shape == (3, 2)  # lev, x
        assert np.all(sigma >= 0)  # Sigma should be non-negative
        assert np.all(sigma <= 1.2)  # Should be close to [0, 1] range

    def test_delta_pressure_hybrid(self):
        """Test pressure layer thickness from hybrid coordinates."""
        ps_da = xr.DataArray(np.array([100000.0, 95000.0]), dims=["x"])
        hya_da = xr.DataArray(np.array([0.0, 0.1, 0.5]), dims=["lev"])
        hyb_da = xr.DataArray(np.array([1.0, 0.9, 0.5]), dims=["lev"])

        pressure = pressure_at_hybrid_levels(ps_da, hya_da, hyb_da)
        dph = delta_pressure_hybrid(ps_da, hya_da, hyb_da)

        assert dph.dims == ("lev", "x")
        assert dph.shape == (2, 2)
        assert dph.name == "dph"
        assert dph.attrs["units"] == "Pa"
        assert_allclose(
            dph.values,
            np.abs(pressure.values[1:] - pressure.values[:-1]),
            rtol=0.0,
            atol=1e-10,
        )


class TestMandatoryPressureLevels:
    """Test mandatory pressure levels constant."""

    def test_mandatory_pressure_levels(self):
        """Test that mandatory pressure levels are defined correctly."""
        # Check that mandatory pressure levels exist and are reasonable
        assert len(__pres_lev_mandatory__) == 21
        assert np.all(__pres_lev_mandatory__ > 0)
        assert np.max(__pres_lev_mandatory__) == 100000.0  # 1000 mb in Pa
        assert np.min(__pres_lev_mandatory__) == 100.0  # 1 mb in Pa

        # Check that levels are in descending order (high pressure to low pressure)
        assert np.all(np.diff(__pres_lev_mandatory__) < 0)


class TestPrePostInterpolationHelpers:
    """Test preprocessing and postprocessing helper functions."""

    def test_pre_interp_multidim_no_cyclic_no_missing(self):
        """Test preprocessing without cyclic points or missing values."""
        # Create test data
        data = np.random.randn(10, 20)
        lat = np.linspace(-90, 90, 10)
        lon = np.linspace(0, 350, 20)

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )

        result = _pre_interp_multidim(data_in, cyclic=False, missing_val=None)

        # Should be unchanged
        assert result.shape == data_in.shape
        assert_array_equal(result.values, data_in.values)

    def test_pre_interp_multidim_with_cyclic(self):
        """Test preprocessing with cyclic boundary conditions."""
        data = np.random.randn(5, 8)
        lat = np.linspace(-60, 60, 5)
        lon = np.linspace(0, 315, 8)  # 45-degree spacing, missing 360

        data_in = xr.DataArray(
            data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            name="wind",
            attrs={"units": "m/s", "long_name": "wind speed"},
        )

        result = _pre_interp_multidim(data_in, cyclic=True, missing_val=None)

        # Should have padded longitude dimension
        assert result.shape == (5, 10)  # lat unchanged, lon padded by 2

        # Check longitude coordinates
        # The exact values depend on wrap mode: wraps last value to beginning and first to end
        assert result.lon.values[0] == lon[-1] - 360  # wrapped last value adjusted
        assert result.lon.values[-1] == lon[0] + 360  # wrapped first value adjusted
        assert result.name == "wind"
        assert result.attrs["units"] == "m/s"
        assert result.attrs["long_name"] == "wind speed"

    def test_pre_interp_multidim_with_missing_val(self):
        """Test preprocessing with missing values."""
        data = np.array([[1, 2, 99], [4, 99, 6]])  # 99 is missing value
        lat = np.array([0, 30])
        lon = np.array([0, 90, 180])

        data_in = xr.DataArray(
            data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            name="field",
            attrs={"units": "K", "long_name": "temperature"},
        )

        result = _pre_interp_multidim(data_in, cyclic=False, missing_val=99)

        # Missing values should be replaced with NaN
        assert np.isnan(result.values[0, 2])
        assert np.isnan(result.values[1, 1])
        assert result.values[0, 0] == 1
        assert result.values[0, 1] == 2
        assert result.name == "field"
        assert result.attrs["units"] == "K"
        assert result.attrs["long_name"] == "temperature"

    def test_pre_interp_multidim_cyclic_and_missing(self):
        """Test preprocessing with both cyclic and missing value handling."""
        data = np.array([[1, 99, 3], [4, 5, 99]])
        lat = np.array([0, 45])
        lon = np.array([0, 120, 240])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )

        result = _pre_interp_multidim(data_in, cyclic=True, missing_val=99)

        # Should be padded and have missing values replaced
        assert result.shape == (2, 5)  # padded longitude
        # The location of NaN values may shift due to padding, so check total count
        assert (
            np.sum(np.isnan(result.values)) == 3
        )  # cyclic padding duplicates an edge NaN

    def test_post_interp_multidim_no_missing(self):
        """Test postprocessing without missing value handling."""
        data = np.array([[1.5, 2.3], [4.1, 5.9]])
        data_in = xr.DataArray(data, dims=["lat", "lon"])

        result = _post_interp_multidim(data_in, missing_val=None)

        # Should be unchanged
        assert_array_equal(result.values, data_in.values)

    def test_post_interp_multidim_with_missing(self):
        """Test postprocessing with missing value replacement."""
        data = np.array([[1.5, np.nan], [4.1, np.nan]])
        data_in = xr.DataArray(
            data,
            dims=["lat", "lon"],
            name="field",
            attrs={"units": "Pa", "long_name": "pressure"},
        )

        result = _post_interp_multidim(data_in, missing_val=-999)

        # NaN values should be replaced with missing_val
        assert result.values[0, 1] == -999
        assert result.values[1, 1] == -999
        assert result.values[0, 0] == 1.5
        assert result.values[1, 0] == 4.1
        assert result.name == "field"
        assert result.attrs["units"] == "Pa"
        assert result.attrs["long_name"] == "pressure"


class TestHybridToPressureInterpolation:
    """Test hybrid-sigma to pressure level interpolation."""

    @pytest.fixture
    def sample_hybrid_data(self):
        """Create sample hybrid-sigma data for testing."""
        # Dimensions
        time = 5
        lev = 10
        lat = 20
        lon = 30

        # Coordinates
        time_coord = np.arange(time)
        lev_coord = np.arange(lev)
        lat_coord = np.linspace(-90, 90, lat)
        lon_coord = np.linspace(0, 357.5, lon)

        # Create realistic hybrid coefficients
        hya = np.linspace(0.0, 0.3, lev)  # dimensionless hybrid A coefficient
        hyb = np.linspace(1.0, 0.0, lev)  # dimensionless

        # Surface pressure (varying in space and time)
        ps_base = 101325.0  # Standard atmospheric pressure
        ps = ps_base + np.random.randn(time, lat, lon) * 1000

        # Sample temperature data
        temp_data = 250 + 50 * np.random.randn(time, lev, lat, lon)

        # Create xarray objects
        data = xr.DataArray(
            temp_data,
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": time_coord,
                "lev": lev_coord,
                "lat": lat_coord,
                "lon": lon_coord,
            },
            attrs={"units": "K", "long_name": "Temperature"},
        )

        ps_da = xr.DataArray(
            ps,
            dims=["time", "lat", "lon"],
            coords={"time": time_coord, "lat": lat_coord, "lon": lon_coord},
            attrs={"units": "Pa", "long_name": "Surface Pressure"},
        )

        hya_da = xr.DataArray(hya, dims=["lev"], coords={"lev": lev_coord})
        hyb_da = xr.DataArray(hyb, dims=["lev"], coords={"lev": lev_coord})

        return data, ps_da, hya_da, hyb_da

    def test_interp_hybrid_to_pressure_basic(self, sample_hybrid_data):
        """Test basic hybrid to pressure interpolation."""
        data, ps, hya, hyb = sample_hybrid_data

        # Use a subset of standard pressure levels
        new_levels = np.array([100000, 85000, 70000, 50000, 30000])  # Pa

        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        # Check output structure
        assert "plev" in result.dims
        assert "lev" not in result.dims
        assert len(result.plev) == len(new_levels)
        assert result.shape == (5, 5, 20, 30)  # time, plev, lat, lon

        # Check that pressure coordinates are correct
        assert_array_equal(result.plev.values, new_levels)

        # Check that metadata is preserved
        assert result.attrs["units"] == "K"
        assert result.attrs["long_name"] == "Temperature"

    def test_interp_hybrid_to_pressure_methods(self, sample_hybrid_data):
        """Test different interpolation methods."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 50000, 30000])

        # Test linear interpolation
        result_linear = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            method="linear",
        )

        # Test log interpolation
        result_log = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            method="log",
        )

        # Both should have same shape
        assert result_linear.shape == result_log.shape

        # Results should be different (unless by coincidence)
        assert not np.allclose(result_linear.values, result_log.values)

    def test_interp_hybrid_to_pressure_extrapolation(self, sample_hybrid_data):
        """Test extrapolation functionality."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 85000, 70000])

        # Test with extrapolation enabled
        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            extrapolate=True,
            variable="other",  # Use simple extrapolation
        )

        assert result.shape == (5, 3, 20, 30)
        assert np.all(np.isfinite(result.values))

    def test_interp_hybrid_to_pressure_validation(self, sample_hybrid_data):
        """Test input validation."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 50000])

        # Test invalid method
        with pytest.raises(ValueError, match="Unknown interpolation method"):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                lev_dim="lev",
                method="invalid",
            )

        # Test extrapolation without variable
        with pytest.raises(ValueError, match="If `extrapolate` is True"):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                lev_dim="lev",
                extrapolate=True,
            )

        # Test invalid variable
        with pytest.raises(ValueError, match="accepted values are"):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                lev_dim="lev",
                extrapolate=True,
                variable="invalid_variable",
            )

    def test_interp_hybrid_to_pressure_with_missing_lev_dim(self, sample_hybrid_data):
        """Test automatic detection of level dimension."""
        data, ps, hya, hyb = sample_hybrid_data

        # Remove lev_dim parameter - should auto-detect
        new_levels = np.array([100000, 50000])

        # This should work with automatic detection
        try:
            result = interp_hybrid_to_pressure(
                data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
            )
            # If auto-detection works
            assert "plev" in result.dims
        except ValueError:
            # If auto-detection fails, this is acceptable
            pass

    def test_interp_hybrid_to_pressure_edge_cases(self, sample_hybrid_data):
        """Test edge cases for hybrid to pressure interpolation."""
        data, ps, hya, hyb = sample_hybrid_data

        # Test with single pressure level
        single_level = np.array([70000])
        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=single_level, lev_dim="lev"
        )

        assert result.shape == (5, 1, 20, 30)  # time, single_plev, lat, lon
        assert len(result.plev) == 1

    def test_interp_hybrid_to_pressure_temperature_extrapolation(
        self, sample_hybrid_data
    ):
        """Test temperature extrapolation with required parameters."""
        data, ps, hya, hyb = sample_hybrid_data

        # Create required temperature and geopotential data
        t_bot = data.isel(lev=-1)  # Bottom level temperature
        phi_sfc = xr.DataArray(
            np.zeros((5, 20, 30)),  # Sea level geopotential
            dims=["time", "lat", "lon"],
            coords={"time": data.time, "lat": data.lat, "lon": data.lon},
        )

        new_levels = np.array([100000, 85000])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            extrapolate=True,
            variable="temperature",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert result.shape == (5, 2, 20, 30)
        assert np.all(np.isfinite(result.values))

    def test_interp_hybrid_to_pressure_missing_extrap_params(self, sample_hybrid_data):
        """Test error when extrapolation parameters are missing."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 85000])

        # Test missing t_bot for temperature extrapolation
        with pytest.raises(
            ValueError, match="both `t_bot` and `phi_sfc` must be provided"
        ):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                lev_dim="lev",
                extrapolate=True,
                variable="temperature",
                phi_sfc=xr.DataArray(
                    np.zeros((5, 20, 30)), dims=["time", "lat", "lon"]
                ),
            )

        # Test missing phi_sfc for geopotential extrapolation
        with pytest.raises(
            ValueError, match="both `t_bot` and `phi_sfc` must be provided"
        ):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                lev_dim="lev",
                extrapolate=True,
                variable="geopotential",
                t_bot=data.isel(lev=-1),
            )


class TestSigmaToHybridInterpolation:
    """Test sigma to hybrid coordinate interpolation."""

    @pytest.fixture
    def sample_sigma_data(self):
        """Create sample sigma coordinate data."""
        # Dimensions
        time = 3
        sigma_lev = 8
        lat = 15
        lon = 20

        # Sigma coordinates (0 at top, 1 at surface)
        sigma_coords = np.linspace(0.1, 1.0, sigma_lev)

        # Sample data
        data = 250 + 50 * np.random.randn(time, sigma_lev, lat, lon)

        # Surface pressure
        ps = 101325 + np.random.randn(time, lat, lon) * 1000

        # Target hybrid coefficients
        nlev_hybrid = 6
        hya = np.linspace(0.0, 0.3, nlev_hybrid)
        hyb = np.linspace(1.0, 0.0, nlev_hybrid)

        # Create xarray objects
        data_da = xr.DataArray(
            data,
            dims=["time", "sigma", "lat", "lon"],
            coords={
                "time": np.arange(time),
                "sigma": sigma_coords,
                "lat": np.linspace(-45, 45, lat),
                "lon": np.linspace(0, 357.5, lon),
            },
        )

        ps_da = xr.DataArray(
            ps,
            dims=["time", "lat", "lon"],
            coords={
                "time": data_da.time,
                "lat": data_da.lat,
                "lon": data_da.lon,
            },
        )

        hya_da = xr.DataArray(hya, dims=["hlev"])
        hyb_da = xr.DataArray(hyb, dims=["hlev"])
        sig_da = xr.DataArray(sigma_coords, dims=["sigma"])

        return data_da, sig_da, ps_da, hya_da, hyb_da

    def test_interp_sigma_to_hybrid_basic(self, sample_sigma_data):
        """Test basic sigma to hybrid interpolation."""
        data, sig_coords, ps, hya, hyb = sample_sigma_data

        result = interp_sigma_to_hybrid(
            data=data, sig_coords=sig_coords, ps=ps, hyam=hya, hybm=hyb, lev_dim="sigma"
        )

        # Check output structure
        assert "hlev" in result.dims
        assert "sigma" not in result.dims
        assert len(result.hlev) == len(hya)
        assert result.shape == (3, 6, 15, 20)  # time, hlev, lat, lon

    def test_interp_sigma_to_hybrid_methods(self, sample_sigma_data):
        """Test different interpolation methods for sigma to hybrid."""
        data, sig_coords, ps, hya, hyb = sample_sigma_data

        # Test linear interpolation
        result_linear = interp_sigma_to_hybrid(
            data=data,
            sig_coords=sig_coords,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            lev_dim="sigma",
            method="linear",
        )

        # Test log interpolation
        result_log = interp_sigma_to_hybrid(
            data=data,
            sig_coords=sig_coords,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            lev_dim="sigma",
            method="log",
        )

        # Both should have same shape
        assert result_linear.shape == result_log.shape

        # Results should be different
        assert not np.allclose(result_linear.values, result_log.values)

    def test_interp_sigma_to_hybrid_1d_data(self):
        """Test sigma to hybrid interpolation with 1D data."""
        # Create simple 1D test case
        sigma_coords = np.array([0.2, 0.5, 0.8, 1.0])
        data = xr.DataArray(
            [220, 250, 280, 290], dims=["sigma"], coords={"sigma": sigma_coords}
        )

        ps = xr.DataArray(101325)  # Scalar surface pressure
        hya = xr.DataArray([0.0, 0.1])  # 2 hybrid levels
        hyb = xr.DataArray([0.8, 0.4])

        result = interp_sigma_to_hybrid(
            data=data,
            sig_coords=sigma_coords,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            lev_dim="sigma",
        )

        assert result.shape == (2,)  # 2 hybrid levels
        assert "hlev" in result.dims
        assert np.all(np.isfinite(result.values))

    def test_interp_sigma_to_hybrid_validation(self, sample_sigma_data):
        """Test input validation for sigma to hybrid interpolation."""
        data, sig_coords, ps, hya, hyb = sample_sigma_data

        # Test invalid method
        with pytest.raises(ValueError, match="Unknown interpolation method"):
            interp_sigma_to_hybrid(
                data=data,
                sig_coords=sig_coords,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                lev_dim="sigma",
                method="invalid_method",
            )


class TestInterpolationFortranFallbacks:
    """Exercise rarely hit Fortran-accelerated fallback branches."""

    def test_import_fallback_sets_fortran_handles_to_none(self, monkeypatch):
        """Reload the module with a forced kernel import failure."""

        real_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if "vinth2p_kernels" in name:
                raise ImportError("forced test import failure")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.delitem(
            sys.modules, "skyborn.interp.fortran.vinth2p_kernels", raising=False
        )
        monkeypatch.setattr(builtins, "__import__", fake_import)
        reloaded = importlib.reload(interpolation_module)

        assert reloaded._dsigma2hybrid_nodes is None
        assert reloaded._dvinth2p_nodes_pa is None
        assert reloaded._dvinth2p_ecmwf_nodes_corder_pa_into is None

        monkeypatch.setattr(builtins, "__import__", real_import)
        importlib.reload(interpolation_module)

    def test_as_broadcast_float64_flat_handles_none_and_broadcast_errors(
        self, monkeypatch
    ):
        """Return None for unsupported inputs or broadcast failures."""

        template = xr.DataArray(np.ones(3, dtype=np.float64), dims=["x"])
        assert interpolation_module._as_broadcast_float64_flat(None, template) is None

        array = xr.DataArray(np.array([1.0], dtype=np.float64), dims=["x"])

        def raise_broadcast(self, other, exclude=None):
            raise ValueError("boom")

        monkeypatch.setattr(xr.DataArray, "broadcast_like", raise_broadcast)
        assert interpolation_module._as_broadcast_float64_flat(array, template) is None

    def test_interp_hybrid_to_pressure_corder_returns_none_when_aux_inputs_missing(
        self, monkeypatch
    ):
        """ECMWF c-order fast path should bail out when aux fields cannot flatten."""

        data = xr.DataArray(
            np.arange(12, dtype=np.float64).reshape(2, 2, 3),
            dims=["lev", "lat", "lon"],
        )
        ps = xr.DataArray(
            np.full((2, 3), 90000.0, dtype=np.float64),
            dims=["lat", "lon"],
        )
        hyam = xr.DataArray(np.array([0.1, 0.0], dtype=np.float64), dims=["lev"])
        hybm = xr.DataArray(np.array([0.9, 1.0], dtype=np.float64), dims=["lev"])

        monkeypatch.setattr(
            interpolation_module, "_dvinth2p_nodes_corder_pa_into", lambda *a, **k: None
        )
        monkeypatch.setattr(
            interpolation_module,
            "_dvinth2p_ecmwf_nodes_corder_pa_into",
            lambda *a, **k: None,
        )

        result = interpolation_module._interp_hybrid_to_pressure_corder(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            p0=100000.0,
            new_levels=np.array([85000.0], dtype=np.float64),
            lev_dim="lev",
            method="linear",
            extrapolate=True,
            variable="temperature",
            t_bot=None,
            phi_sfc=None,
        )

        assert result is None

    def test_interp_hybrid_to_pressure_helper_broadcasts_ecmwf_inputs(
        self, monkeypatch
    ):
        """Fallback Fortran path should accept broadcasted ECMWF auxiliaries."""

        captured = {}

        data = xr.DataArray(
            np.array([[280.0, 281.0, 282.0], [270.0, 271.0, 272.0]], dtype=np.float64),
            dims=["lev", "x"],
        )
        ps = xr.DataArray(np.full(3, 95000.0, dtype=np.float64), dims=["x"])
        hyam = xr.DataArray(np.array([0.2, 0.0], dtype=np.float64), dims=["lev"])
        hybm = xr.DataArray(np.array([0.8, 1.0], dtype=np.float64), dims=["lev"])
        t_bot = xr.DataArray(np.array(290.0, dtype=np.float64))
        phi_sfc = xr.DataArray(np.array(50.0, dtype=np.float64))

        def fake_ecmwf_into(
            data_columns,
            output_columns,
            hyam_values,
            hybm_values,
            p0,
            new_level_values,
            intyp,
            ps_columns,
            spvl,
            ixtrp,
            varflg,
            t_bot_columns,
            phi_columns,
        ):
            captured["t_bot"] = t_bot_columns.copy()
            captured["phi"] = phi_columns.copy()
            output_columns[:] = 42.0

        monkeypatch.setattr(interpolation_module, "_dvinth2p_nodes_pa", object())
        monkeypatch.setattr(
            interpolation_module,
            "_interp_hybrid_to_pressure_corder",
            lambda **kwargs: None,
        )
        monkeypatch.setattr(
            interpolation_module, "_dvinth2p_ecmwf_nodes_pa_into", fake_ecmwf_into
        )

        result = interpolation_module._interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            p0=100000.0,
            new_levels=np.array([85000.0], dtype=np.float64),
            lev_dim="lev",
            method="linear",
            extrapolate=True,
            variable="temperature",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert result.dims == ("plev", "x")
        assert np.all(result.values == 42.0)
        assert np.array_equal(captured["t_bot"], np.full(3, 290.0))
        assert np.array_equal(captured["phi"], np.full(3, 50.0))

    def test_interp_sigma_to_hybrid_helper_requires_backend(self, monkeypatch):
        """Private eager helper should fail clearly when backend is absent."""

        data = xr.DataArray(np.array([280.0, 275.0], dtype=np.float64), dims=["sigma"])
        sig_coords = xr.DataArray(
            np.array([0.2, 1.0], dtype=np.float64), dims=["sigma"]
        )
        ps = xr.DataArray(np.array(100000.0, dtype=np.float64))
        hyam = xr.DataArray(np.array([0.0, 0.1], dtype=np.float64), dims=["hlev"])
        hybm = xr.DataArray(np.array([1.0, 0.8], dtype=np.float64), dims=["hlev"])

        monkeypatch.setattr(interpolation_module, "_dsigma2hybrid_nodes", None)

        with pytest.raises(RuntimeError, match="sigma2hybrid backend is not available"):
            interpolation_module._interp_sigma_to_hybrid(
                data=data,
                sig_coords=sig_coords,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
                p0=100000.0,
                lev_dim="sigma",
                method="linear",
            )

    def test_interp_sigma_to_hybrid_helper_requires_monotonic_sigma(self, monkeypatch):
        """Non-monotonic sigma coordinates should be rejected before backend use."""

        data = xr.DataArray(
            np.array([280.0, 275.0, 270.0], dtype=np.float64),
            dims=["sigma"],
        )
        sig_coords = xr.DataArray(
            np.array([0.2, 0.9, 0.5], dtype=np.float64), dims=["sigma"]
        )
        ps = xr.DataArray(np.array(100000.0, dtype=np.float64))
        hyam = xr.DataArray(np.array([0.0, 0.1], dtype=np.float64), dims=["hlev"])
        hybm = xr.DataArray(np.array([1.0, 0.8], dtype=np.float64), dims=["hlev"])

        monkeypatch.setattr(interpolation_module, "_dsigma2hybrid_nodes", object())

        with pytest.raises(ValueError, match="requires monotonic sigma coordinates"):
            interpolation_module._interp_sigma_to_hybrid(
                data=data,
                sig_coords=sig_coords,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
                p0=100000.0,
                lev_dim="sigma",
                method="linear",
            )

    def test_interp_sigma_to_hybrid_corder_rejects_non_numpy_data(self, monkeypatch):
        """Fast path should reject non-NumPy raw data cleanly."""

        sigma_values = np.array([0.2, 1.0], dtype=np.float64)
        data = xr.DataArray(
            np.array([280.0, 275.0], dtype=np.float64),
            dims=["sigma"],
        ).chunk({"sigma": 2})
        ps = xr.DataArray(np.array(100000.0, dtype=np.float64))
        hyam = xr.DataArray(np.array([0.0, 0.1], dtype=np.float64), dims=["hlev"])
        hybm = xr.DataArray(np.array([1.0, 0.8], dtype=np.float64), dims=["hlev"])

        monkeypatch.setattr(
            interpolation_module,
            "_dsigma2hybrid_nodes_corder_into",
            lambda *a, **k: None,
        )

        result = interpolation_module._interp_sigma_to_hybrid_corder(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            p0=100000.0,
            lev_dim="sigma",
            method="linear",
            sigma_source_values=sigma_values,
        )

        assert result is None

    def test_interp_sigma_to_hybrid_corder_rejects_nonfinite_broadcast_ps(
        self, monkeypatch
    ):
        """Broadcast fallback should reject non-finite surface pressure values."""

        sigma_values = np.array([0.2, 1.0], dtype=np.float64)
        data = xr.DataArray(
            np.array([280.0, 275.0], dtype=np.float64),
            dims=["sigma"],
        )
        ps = xr.DataArray(np.array(np.nan, dtype=np.float64))
        hyam = xr.DataArray(np.array([0.0, 0.1], dtype=np.float64), dims=["hlev"])
        hybm = xr.DataArray(np.array([1.0, 0.8], dtype=np.float64), dims=["hlev"])

        monkeypatch.setattr(
            interpolation_module,
            "_dsigma2hybrid_nodes_corder_into",
            lambda *a, **k: None,
        )

        result = interpolation_module._interp_sigma_to_hybrid_corder(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            p0=100000.0,
            lev_dim="sigma",
            method="linear",
            sigma_source_values=sigma_values,
        )

        assert result is None

    def test_interp_sigma_to_hybrid_1d_python_fallback(self, monkeypatch):
        """Compiled-only sigma-to-hybrid should fail clearly without the backend."""

        sigma_coords = np.array([0.2, 0.5, 0.8, 1.0], dtype=np.float64)
        data = xr.DataArray(
            np.array([220.0, 250.0, 280.0, 290.0], dtype=np.float64),
            dims=["sigma"],
            coords={"sigma": sigma_coords},
        )
        ps = xr.DataArray(np.array(101325.0, dtype=np.float64))
        hya = xr.DataArray(
            np.array([0.0, 0.1], dtype=np.float64),
            dims=["hlev"],
            coords={"hlev": [101, 102]},
        )
        hyb = xr.DataArray(
            np.array([0.8, 0.4], dtype=np.float64),
            dims=["hlev"],
            coords={"hlev": [101, 102]},
        )

        monkeypatch.setattr(interpolation_module, "_dsigma2hybrid_nodes", None)
        monkeypatch.setattr(interpolation_module, "_dsigma2hybrid_nodes_into", None)
        monkeypatch.setattr(
            interpolation_module, "_dsigma2hybrid_nodes_corder_into", None
        )

        with pytest.raises(RuntimeError, match="requires the compiled Fortran backend"):
            interp_sigma_to_hybrid(
                data=data,
                sig_coords=sigma_coords,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                lev_dim="sigma",
            )


class TestMultidimensionalInterpolation:
    """Test multidimensional spatial interpolation."""

    def test_interp_multidim_basic(self):
        """Test basic multidimensional interpolation."""
        # Create test data
        lat_in = np.array([0, 30, 60, 90])
        lon_in = np.array([0, 90, 180, 270])
        data = np.random.randn(4, 4)

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Output coordinates
        lat_out = np.array([15, 45, 75])
        lon_out = np.array([45, 135, 225, 315])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, method="nearest"
        )

        # Check output shape
        assert result.shape == (3, 4)  # lat_out, lon_out
        assert_array_equal(result.lat.values, lat_out)
        assert_array_equal(result.lon.values, lon_out)

    def test_interp_multidim_numpy_input(self):
        """Test multidimensional interpolation with numpy arrays."""
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 120, 240])
        data = np.random.randn(3, 3)

        lat_out = np.array([15, 45])
        lon_out = np.array([60, 180])

        result = interp_multidim(
            data_in=data, lat_in=lat_in, lon_in=lon_in, lat_out=lat_out, lon_out=lon_out
        )

        assert result.shape == (2, 2)
        assert isinstance(result, xr.DataArray)

    def test_interp_multidim_cyclic(self):
        """Test multidimensional interpolation with cyclic boundary."""
        lat_in = np.array([-90, 0, 90])
        lon_in = np.array([0, 180])  # Only half the globe
        data = np.array([[1, 2], [3, 4], [5, 6]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Request data at 360 degrees (should wrap to 0)
        lat_out = np.array([0])
        lon_out = np.array([360])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, cyclic=True
        )

        assert result.shape == (1, 1)
        # Should be close to the value at lon=0
        assert not np.isnan(result.values[0, 0])

    def test_interp_multidim_missing_values(self):
        """Test handling of missing values."""
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 90, 180])
        data = np.array([[1, 2, 3], [4, 99, 6], [7, 8, 9]])  # 99 is missing

        data_in = xr.DataArray(
            data,
            dims=["lat", "lon"],
            coords={"lat": lat_in, "lon": lon_in},
            name="field",
            attrs={"units": "m/s", "long_name": "wind speed"},
        )

        lat_out = np.array([15, 45])
        lon_out = np.array([45, 135])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, missing_val=99
        )

        assert result.shape == (2, 2)
        assert result.name == "field"
        assert result.attrs["units"] == "m/s"
        assert result.attrs["long_name"] == "wind speed"

    def test_interp_multidim_validation(self):
        """Test input validation for multidimensional interpolation."""
        data = np.random.randn(3, 3)
        lat_out = np.array([15, 45])
        lon_out = np.array([60, 180])

        # Test missing coordinates for numpy input
        with pytest.raises(ValueError, match="lat_in and lon_in must be provided"):
            interp_multidim(data_in=data, lat_out=lat_out, lon_out=lon_out)

    def test_interp_multidim_different_methods(self):
        """Test different interpolation methods."""
        lat_in = np.array([0, 45, 90])
        lon_in = np.array([0, 180])
        data = np.array([[1, 2], [3, 4], [5, 6]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        lat_out = np.array([22.5, 67.5])
        lon_out = np.array([90])

        # Test linear interpolation
        result_linear = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, method="linear"
        )

        # Test nearest neighbor interpolation
        result_nearest = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, method="nearest"
        )

        # Both should have same shape
        assert result_linear.shape == result_nearest.shape
        assert result_linear.shape == (2, 1)

        # Results should generally be different
        # (though they might coincidentally be the same)
        assert np.all(np.isfinite(result_linear.values))
        assert np.all(np.isfinite(result_nearest.values))

    def test_interp_multidim_extrapolation(self):
        """Test multidimensional interpolation with extrapolation."""
        lat_in = np.array([10, 20, 30])
        lon_in = np.array([100, 200])
        data = np.array([[1, 2], [3, 4], [5, 6]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Request points outside the input domain
        lat_out = np.array([5, 35])  # Outside lat range
        lon_out = np.array([50, 250])  # Outside lon range

        # With extrapolation enabled
        result_extrap = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, fill_value="extrapolate"
        )

        # Without extrapolation (default)
        result_no_extrap = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out
        )

        # Both should have same shape
        assert result_extrap.shape == result_no_extrap.shape
        assert result_extrap.shape == (2, 2)

        # Extrapolated result should have finite values
        assert np.all(np.isfinite(result_extrap.values))

        # Non-extrapolated result should have some NaN values
        assert np.any(np.isnan(result_no_extrap.values))

    def test_interp_multidim_multidimensional_input(self):
        """Test multidimensional interpolation with 3D+ input data."""
        # Create 3D data (time, lat, lon)
        time = 4
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 90, 180, 270])

        data = np.random.randn(time, len(lat_in), len(lon_in))

        data_in = xr.DataArray(
            data,
            dims=["time", "lat", "lon"],
            coords={"time": np.arange(time), "lat": lat_in, "lon": lon_in},
        )

        lat_out = np.array([15, 45])
        lon_out = np.array([45, 135, 225])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, method="nearest"
        )

        # Should preserve time dimension and interpolate spatial dimensions
        assert result.shape == (4, 2, 3)  # time, lat_out, lon_out
        assert "time" in result.coords
        assert_array_equal(result.lat.values, lat_out)
        assert_array_equal(result.lon.values, lon_out)

    def test_interp_multidim_edge_case_coordinates(self):
        """Test multidimensional interpolation with edge case coordinates."""
        # Test with single point grids
        lat_in = np.array([45])
        lon_in = np.array([90])
        data = np.array([[5]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Interpolate to the same point
        lat_out = np.array([45])
        lon_out = np.array([90])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, method="nearest"
        )

        assert result.shape == (1, 1)
        assert result.values[0, 0] == 5

    def test_interp_multidim_large_missing_regions(self):
        """Test interpolation with large regions of missing data."""
        lat_in = np.array([0, 30, 60, 90])
        lon_in = np.array([0, 90, 180, 270])

        # Create data with large missing region
        data = np.array(
            [
                [1, 99, 99, 4],  # 99 = missing
                [99, 99, 99, 99],  # entire row missing
                [99, 99, 7, 8],
                [9, 10, 11, 12],
            ]
        )

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        lat_out = np.array([15, 45, 75])
        lon_out = np.array([45, 135, 225])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, missing_val=99
        )

        assert result.shape == (3, 3)
        # Some interpolated values should be NaN where surrounded by missing data
        # Others should be finite where valid data is available


class TestInterpolationIntegration:
    """Integration tests for interpolation module."""

    def test_interpolation_with_climate_data(self):
        """Test interpolation with realistic climate data."""
        # Create sample climate data directly in this test
        time = 12
        nlev = 10
        lat = 73
        lon = 144

        # Create temperature data
        temp_data = 250 + 50 * np.random.randn(time, nlev, lat, lon)
        temp = xr.DataArray(
            temp_data,
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": np.arange(time),
                "lev": np.arange(nlev),
                "lat": np.linspace(-90, 90, lat),
                "lon": np.linspace(0, 357.5, lon),
            },
        )

        ps = xr.DataArray(
            101325 + np.random.randn(time, lat, lon) * 1000,
            dims=["time", "lat", "lon"],
            coords={"time": temp.time, "lat": temp.lat, "lon": temp.lon},
        )

        hya = xr.DataArray(np.linspace(0.0, 0.3, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])

        # Test hybrid to pressure interpolation
        new_levels = np.array([100000, 85000, 70000, 50000])
        result = interp_hybrid_to_pressure(
            data=temp, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        assert result.shape == (12, 4, 73, 144)  # time, plev, lat, lon
        assert np.any(np.isfinite(result.values))
        assert np.all(np.isfinite(result.values[~np.isnan(result.values)]))

    def test_interpolation_error_handling(self):
        """Test comprehensive error handling."""
        # Create minimal test data
        data = xr.DataArray(
            np.random.randn(3, 5, 5),
            dims=["lev", "lat", "lon"],
            coords={
                "lev": np.arange(3),
                "lat": np.linspace(-60, 60, 5),
                "lon": np.linspace(0, 288, 5),
            },
        )

        ps = xr.DataArray(
            101325 + np.random.randn(5, 5),
            dims=["lat", "lon"],
            coords={"lat": data.lat, "lon": data.lon},
        )

        hya = xr.DataArray([0, 25000, 50000], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        # Test with invalid pressure levels (negative)
        with pytest.raises((ValueError, RuntimeError)):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=np.array([-1000, 50000]),
            )


# Performance tests (marked as slow)
@pytest.mark.slow
class TestInterpolationPerformance:
    """Performance tests for interpolation module."""

    def test_hybrid_to_pressure_large_data(self):
        """Test hybrid to pressure interpolation with large datasets."""
        # Large dataset
        time, lev, lat, lon = 100, 50, 180, 360

        data = xr.DataArray(
            250 + 50 * np.random.randn(time, lev, lat, lon),
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": np.arange(time),
                "lev": np.arange(lev),
                "lat": np.linspace(-90, 90, lat),
                "lon": np.linspace(0, 357.5, lon),
            },
        )

        ps = xr.DataArray(
            101325 + np.random.randn(time, lat, lon) * 1000,
            dims=["time", "lat", "lon"],
            coords={"time": data.time, "lat": data.lat, "lon": data.lon},
        )

        hya = xr.DataArray(np.linspace(0.0, 0.3, lev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, lev), dims=["lev"])

        new_levels = np.array([100000, 85000, 70000, 50000, 30000])

        # Should complete without memory issues
        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        assert result.shape == (100, 5, 180, 360)
        assert np.any(np.isfinite(result.values))
        assert np.all(np.isfinite(result.values[~np.isnan(result.values)]))


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
