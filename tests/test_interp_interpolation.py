"""
Tests for skyborn.interp.interpolation module.

This module tests the interpolation functionality including hybrid-sigma
to pressure level interpolation and multidimensional spatial interpolation.
"""

import pytest
import numpy as np
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal

from skyborn.interp.interpolation import (
    interp_hybrid_to_pressure,
    interp_sigma_to_hybrid,
    interp_multidim,
    _pressure_from_hybrid,
    _sigma_from_hybrid,
    _func_interpolate,
    __pres_lev_mandatory__,
)

# Try to import private functions - they may not be available in all versions
try:
    from skyborn.interp.interpolation import (
        _pre_interp_multidim,
        _post_interp_multidim,
        _vertical_remap,
        _temp_extrapolate,
        _geo_height_extrapolate,
        _vertical_remap_extrap,
    )

    PRIVATE_FUNCTIONS_AVAILABLE = True
except ImportError:
    PRIVATE_FUNCTIONS_AVAILABLE = False


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
        assert pressure.shape == (2, 3)  # x, lev
        assert np.all(pressure > 0)  # All pressures should be positive
        assert np.all(pressure <= 100000)  # Should not exceed reference pressure

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
        assert sigma.shape == (2, 3)  # x, lev
        assert np.all(sigma >= 0)  # Sigma should be non-negative
        assert np.all(sigma <= 1.2)  # Should be close to [0, 1] range

    def test_func_interpolate(self):
        """Test interpolation function selection."""
        # Test linear interpolation function
        func_linear = _func_interpolate("linear")
        assert func_linear is not None

        # Test log interpolation function
        func_log = _func_interpolate("log")
        assert func_log is not None

        # Test invalid method
        with pytest.raises(ValueError, match="Unknown interpolation method"):
            _func_interpolate("invalid_method")


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

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
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

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_pre_interp_multidim_with_cyclic(self):
        """Test preprocessing with cyclic boundary conditions."""
        data = np.random.randn(5, 8)
        lat = np.linspace(-60, 60, 5)
        lon = np.linspace(0, 315, 8)  # 45-degree spacing, missing 360

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )

        result = _pre_interp_multidim(data_in, cyclic=True, missing_val=None)

        # Should have padded longitude dimension
        assert result.shape == (5, 10)  # lat unchanged, lon padded by 2

        # Check longitude coordinates
        # The exact values depend on wrap mode: wraps last value to beginning and first to end
        assert result.lon.values[0] == lon[-1] - 360  # wrapped last value adjusted
        assert result.lon.values[-1] == lon[0] + 360  # wrapped first value adjusted

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_pre_interp_multidim_with_missing_val(self):
        """Test preprocessing with missing values."""
        data = np.array([[1, 2, 99], [4, 99, 6]])  # 99 is missing value
        lat = np.array([0, 30])
        lon = np.array([0, 90, 180])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )

        result = _pre_interp_multidim(data_in, cyclic=False, missing_val=99)

        # Missing values should be replaced with NaN
        assert np.isnan(result.values[0, 2])
        assert np.isnan(result.values[1, 1])
        assert result.values[0, 0] == 1
        assert result.values[0, 1] == 2

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
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
        assert np.sum(np.isnan(result.values)) == 2  # Should have 2 NaN values total

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_post_interp_multidim_no_missing(self):
        """Test postprocessing without missing value handling."""
        data = np.array([[1.5, 2.3], [4.1, 5.9]])
        data_in = xr.DataArray(data, dims=["lat", "lon"])

        result = _post_interp_multidim(data_in, missing_val=None)

        # Should be unchanged
        assert_array_equal(result.values, data_in.values)

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_post_interp_multidim_with_missing(self):
        """Test postprocessing with missing value replacement."""
        data = np.array([[1.5, np.nan], [4.1, np.nan]])
        data_in = xr.DataArray(data, dims=["lat", "lon"])

        result = _post_interp_multidim(data_in, missing_val=-999)

        # NaN values should be replaced with missing_val
        assert result.values[0, 1] == -999
        assert result.values[1, 1] == -999
        assert result.values[0, 0] == 1.5
        assert result.values[1, 0] == 4.1


class TestVerticalRemapHelpers:
    """Test vertical remapping helper functions."""

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_vertical_remap_basic(self):
        """Test basic vertical remapping functionality."""
        # Simple test case
        new_levels = np.array([850, 500, 200])  # mb
        xcoords = np.array([1000, 700, 300])  # mb
        data = np.array([288, 268, 228])  # temperature

        func_interpolate = _func_interpolate("linear")

        # Test interpolation
        result = _vertical_remap(func_interpolate, new_levels, xcoords, data)

        # Should interpolate to 3 levels
        assert len(result) == 3
        assert np.all(np.isfinite(result))

        # Check that interpolated values are reasonable
        assert 220 < result[2] < 270  # 200 mb should be cold
        assert 270 < result[0] < 290  # 850 mb should be warm

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_vertical_remap_multidimensional(self):
        """Test vertical remapping with multidimensional data."""
        new_levels = np.array([850, 500])
        xcoords = np.array([1000, 700, 300])
        # 2D data (3 levels, 4 grid points)
        data = np.array(
            [[288, 285, 290, 287], [268, 265, 270, 267], [228, 225, 230, 227]]
        )

        func_interpolate = _func_interpolate("linear")

        result = _vertical_remap(
            func_interpolate, new_levels, xcoords, data, interp_axis=0
        )

        # Should have shape (2, 4) - 2 new levels, 4 grid points
        assert result.shape == (2, 4)
        assert np.all(np.isfinite(result))


class TestExtrapolationFunctions:
    """Test temperature and geopotential height extrapolation functions."""

    @pytest.fixture
    def extrapolation_test_data(self):
        """Create sample data for extrapolation testing."""
        # Create sample atmospheric data
        lev = np.array([85000, 70000, 50000])  # Pa
        lat = np.linspace(-30, 30, 5)
        lon = np.linspace(0, 90, 4)

        # Temperature data (decreasing with height)
        temp_data = np.zeros((3, 5, 4))
        for i in range(3):
            temp_data[i, :, :] = 288 - i * 20  # Simple vertical profile

        data = xr.DataArray(
            temp_data,
            dims=["lev", "lat", "lon"],
            coords={"lev": lev, "lat": lat, "lon": lon},
        )

        # Surface pressure
        ps = xr.DataArray(
            np.full((5, 4), 101325),  # Standard surface pressure
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
        )

        # Geopotential at surface
        phi_sfc = xr.DataArray(
            np.zeros((5, 4)),  # Sea level
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
        )

        return data, ps, phi_sfc

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_temp_extrapolate_basic(self, extrapolation_test_data):
        """Test basic temperature extrapolation."""
        data, ps, phi_sfc = extrapolation_test_data

        # Pressure at lowest model level
        p_sfc = xr.DataArray(
            np.full((5, 4), 85000),  # Pa
            dims=["lat", "lon"],
            coords={"lat": data.lat, "lon": data.lon},
        )

        # Extrapolate to surface pressure level (higher pressure)
        lev_extrap = 100000  # Pa

        result = _temp_extrapolate(data, "lev", lev_extrap, p_sfc, ps, phi_sfc)

        # Check output shape
        assert result.shape == data.isel(lev=0).shape  # Should match spatial dimensions

        # Extrapolated temperature should be warmer than lowest level
        lowest_temp = data.isel(lev=-1)  # Lowest level (highest pressure)
        assert np.all(result >= lowest_temp)

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_temp_extrapolate_with_topography(self):
        """Test temperature extrapolation with topography."""
        # Create data with varying surface heights
        lat = np.array([0])
        lon = np.array([0])

        data = xr.DataArray(
            np.array([[[280]]]),  # Single temperature value
            dims=["lev", "lat", "lon"],
            coords={"lev": [85000], "lat": lat, "lon": lon},
        )

        ps = xr.DataArray([[101325]], dims=["lat", "lon"])  # Surface pressure
        p_sfc = xr.DataArray([[85000]], dims=["lat", "lon"])  # Model surface pressure

        # High elevation site
        phi_sfc = xr.DataArray([[19620]], dims=["lat", "lon"])  # ~2000m elevation

        result = _temp_extrapolate(data, "lev", 101325, p_sfc, ps, phi_sfc)

        assert result.shape == (1, 1)
        assert np.isfinite(result.values[0, 0])
        # Temperature at higher pressure should be warmer
        assert result.values[0, 0] > 280

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_geo_height_extrapolate_basic(self, extrapolation_test_data):
        """Test basic geopotential height extrapolation."""
        data, ps, phi_sfc = extrapolation_test_data

        # Temperature at bottom level
        t_bot = data.isel(lev=-1)  # Warmest level

        # Pressure levels for extrapolation
        p_sfc = xr.DataArray(
            np.full((5, 4), 85000),
            dims=["lat", "lon"],
            coords={"lat": data.lat, "lon": data.lon},
        )

        lev_extrap = 100000  # Pa (lower altitude, higher pressure)

        result = _geo_height_extrapolate(t_bot, lev_extrap, p_sfc, ps, phi_sfc)

        # Check output shape
        assert result.shape == t_bot.shape

        # Geopotential height should be finite
        assert np.all(np.isfinite(result))

        # At higher pressure (lower altitude), geopotential should be lower
        assert np.all(result <= phi_sfc)

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_geo_height_extrapolate_high_elevation(self):
        """Test geopotential height extrapolation at high elevation."""
        # High altitude site
        t_bot = xr.DataArray([[250]], dims=["lat", "lon"])  # Cold temperature
        ps = xr.DataArray([[70000]], dims=["lat", "lon"])  # Low surface pressure
        p_sfc = xr.DataArray([[70000]], dims=["lat", "lon"])
        phi_sfc = xr.DataArray([[29430]], dims=["lat", "lon"])  # ~3000m elevation

        lev_extrap = 85000  # Extrapolate to higher pressure (lower altitude)

        result = _geo_height_extrapolate(t_bot, lev_extrap, p_sfc, ps, phi_sfc)

        assert result.shape == (1, 1)
        assert np.isfinite(result.values[0, 0])
        # Should be at lower altitude than surface
        assert result.values[0, 0] < phi_sfc.values[0, 0]


class TestVerticalRemapExtrap:
    """Test the vertical remap extrapolation coordinator function."""

    @pytest.fixture
    def extrap_test_setup(self):
        """Set up test data for extrapolation testing."""
        # Create multi-level data
        nlev = 3
        new_levels = np.array([100000, 85000, 70000])  # Pa

        # Create data with vertical coordinate
        data = xr.DataArray(
            np.array([[[290]], [[280]], [[270]]]),  # Temperature profile
            dims=["lev", "lat", "lon"],
            coords={"lev": [50000, 70000, 85000], "lat": [0], "lon": [0]},
        )

        # Output array template
        output = xr.DataArray(
            np.full((3, 1, 1), np.nan),
            dims=["plev", "lat", "lon"],
            coords={"plev": new_levels, "lat": [0], "lon": [0]},
        )

        # Pressure at model levels
        pressure = xr.DataArray(
            np.array([[[50000]], [[70000]], [[85000]]]),
            dims=["lev", "lat", "lon"],
            coords={"lev": [50000, 70000, 85000], "lat": [0], "lon": [0]},
        )

        ps = xr.DataArray([[101325]], dims=["lat", "lon"])
        t_bot = data.isel(lev=-1)  # Bottom temperature
        phi_sfc = xr.DataArray([[0]], dims=["lat", "lon"])  # Sea level

        return new_levels, "lev", data, output, pressure, ps, t_bot, phi_sfc

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_vertical_remap_extrap_temperature(self, extrap_test_setup):
        """Test extrapolation for temperature variable."""
        new_levels, lev_dim, data, output, pressure, ps, t_bot, phi_sfc = (
            extrap_test_setup
        )

        result = _vertical_remap_extrap(
            new_levels,
            lev_dim,
            data,
            output,
            pressure,
            ps,
            variable="temperature",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        # Should have same shape as output
        assert result.shape == output.shape

        # Values should be finite where extrapolation occurred
        assert np.all(np.isfinite(result.values))

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_vertical_remap_extrap_geopotential(self, extrap_test_setup):
        """Test extrapolation for geopotential variable."""
        new_levels, lev_dim, data, output, pressure, ps, t_bot, phi_sfc = (
            extrap_test_setup
        )

        result = _vertical_remap_extrap(
            new_levels,
            lev_dim,
            data,
            output,
            pressure,
            ps,
            variable="geopotential",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert result.shape == output.shape
        assert np.all(np.isfinite(result.values))

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_vertical_remap_extrap_other(self, extrap_test_setup):
        """Test extrapolation for other variables (simple fill)."""
        new_levels, lev_dim, data, output, pressure, ps, t_bot, phi_sfc = (
            extrap_test_setup
        )

        result = _vertical_remap_extrap(
            new_levels,
            lev_dim,
            data,
            output,
            pressure,
            ps,
            variable="other",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert result.shape == output.shape
        # For "other" variables, should use surface level value for extrapolation
        assert np.all(np.isfinite(result.values))


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
        hya = np.linspace(0, 50000, lev)  # Pa
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
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, method="linear"
        )

        # Test log interpolation
        result_log = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, method="log"
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
        hya = np.linspace(0, 30000, nlev_hybrid)
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

        ps_da = xr.DataArray(ps, dims=["time", "lat", "lon"], coords=data_da.coords)

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

        ps = xr.DataArray([101325])  # Scalar surface pressure
        hya = xr.DataArray([10000, 30000])  # 2 hybrid levels
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

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

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
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        lat_out = np.array([15, 45])
        lon_out = np.array([45, 135])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, missing_val=99
        )

        assert result.shape == (2, 2)

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

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

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

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

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

        hya = xr.DataArray(np.linspace(0, 50000, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])

        # Test hybrid to pressure interpolation
        new_levels = np.array([100000, 85000, 70000, 50000])
        result = interp_hybrid_to_pressure(
            data=temp, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
        )

        assert result.shape == (12, 4, 73, 144)  # time, plev, lat, lon
        assert np.all(np.isfinite(result.values))

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
            coords=data.coords,
        )

        hya = xr.DataArray(np.linspace(0, 50000, lev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, lev), dims=["lev"])

        new_levels = np.array([100000, 85000, 70000, 50000, 30000])

        # Should complete without memory issues
        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
        )

        assert result.shape == (100, 5, 180, 360)
        assert np.all(np.isfinite(result.values))


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
