"""
Comprehensive tests for GPI xarray interface to improve code coverage.
"""

import warnings

import numpy as np
import pytest

try:
    import xarray as xr

    _has_xarray = True
except ImportError:
    _has_xarray = False


@pytest.fixture
def sample_3d_data():
    """Create sample 3D atmospheric data for testing."""
    nlat, nlon = 10, 12
    nlevs = 8

    # Pressure levels (mb)
    levels = np.array([1000, 925, 850, 700, 500, 300, 200, 100])

    # Create realistic data
    np.random.seed(42)

    # Temperature decreases with height
    temp_base = 300 - np.arange(nlevs) * 10
    temperature = np.zeros((nlevs, nlat, nlon))
    for i in range(nlevs):
        temperature[i] = temp_base[i] + np.random.randn(nlat, nlon) * 2

    # Mixing ratio decreases with height
    mixr_base = 0.015 * np.exp(-np.arange(nlevs) / 3)
    mixing_ratio = np.zeros((nlevs, nlat, nlon))
    for i in range(nlevs):
        mixing_ratio[i] = mixr_base[i] + np.random.randn(nlat, nlon) * 0.001
        mixing_ratio[i] = np.maximum(mixing_ratio[i], 0.00001)  # Ensure positive

    # Surface fields
    sst = 298 + np.random.randn(nlat, nlon) * 2  # K
    psl = 101325 + np.random.randn(nlat, nlon) * 500  # Pa

    return {
        "levels": levels,
        "temperature": temperature,
        "mixing_ratio": mixing_ratio,
        "sst": sst,
        "psl": psl,
        "nlat": nlat,
        "nlon": nlon,
        "nlevs": nlevs,
    }


@pytest.fixture
def sample_4d_data():
    """Create sample 4D atmospheric data for testing."""
    ntimes = 5
    nlat, nlon = 8, 10
    nlevs = 6

    # Pressure levels (mb)
    levels = np.array([1000, 850, 700, 500, 300, 200])

    # Create realistic data
    np.random.seed(42)

    # Temperature
    temp_base = 300 - np.arange(nlevs) * 12
    temperature = np.zeros((ntimes, nlevs, nlat, nlon))
    for t in range(ntimes):
        for i in range(nlevs):
            temperature[t, i] = temp_base[i] + np.random.randn(nlat, nlon) * 2

    # Mixing ratio
    mixr_base = 0.015 * np.exp(-np.arange(nlevs) / 3)
    mixing_ratio = np.zeros((ntimes, nlevs, nlat, nlon))
    for t in range(ntimes):
        for i in range(nlevs):
            mixing_ratio[t, i] = mixr_base[i] + np.random.randn(nlat, nlon) * 0.001
            mixing_ratio[t, i] = np.maximum(mixing_ratio[t, i], 0.00001)

    # Surface fields
    sst = np.zeros((ntimes, nlat, nlon))
    psl = np.zeros((ntimes, nlat, nlon))
    for t in range(ntimes):
        sst[t] = 298 + np.random.randn(nlat, nlon) * 2
        psl[t] = 101325 + np.random.randn(nlat, nlon) * 500

    return {
        "levels": levels,
        "temperature": temperature,
        "mixing_ratio": mixing_ratio,
        "sst": sst,
        "psl": psl,
        "ntimes": ntimes,
        "nlat": nlat,
        "nlon": nlon,
        "nlevs": nlevs,
    }


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestUnitConversions:
    """Test unit conversion functionality."""

    def test_temperature_celsius_to_kelvin(self):
        """Test temperature conversion from Celsius to Kelvin."""
        from skyborn.calc.GPI.xarray import _check_units

        # Test various Celsius unit formats
        celsius_units = ["C", "°C", "celsius", "Celsius", "degC", "deg C"]

        for unit in celsius_units:
            temp_c = xr.DataArray(25.0, attrs={"units": unit})
            temp_k, converted = _check_units(temp_c, "temperature", ["K"])

            assert converted == True
            assert np.isclose(temp_k.values, 298.15)
            assert temp_k.attrs["units"] == "K"

    def test_temperature_fahrenheit_to_kelvin(self):
        """Test temperature conversion from Fahrenheit to Kelvin."""
        from skyborn.calc.GPI.xarray import _check_units

        # Test various Fahrenheit unit formats
        fahrenheit_units = ["F", "°F", "fahrenheit", "Fahrenheit", "degF"]

        for unit in fahrenheit_units:
            temp_f = xr.DataArray(77.0, attrs={"units": unit})  # 77°F = 25°C = 298.15K
            temp_k, converted = _check_units(temp_f, "temperature", ["K"])

            assert converted == True
            assert np.isclose(temp_k.values, 298.15, rtol=0.01)
            assert temp_k.attrs["units"] == "K"

    def test_temperature_kelvin_standardization(self):
        """Test Kelvin unit standardization."""
        from skyborn.calc.GPI.xarray import _check_units

        kelvin_units = ["kelvin", "Kelvin", "KELVIN", "degK"]

        for unit in kelvin_units:
            temp = xr.DataArray(298.15, attrs={"units": unit})
            temp_std, converted = _check_units(temp, "temperature", ["K"])

            assert converted == False or unit != "K"  # Only convert if not already "K"
            assert temp_std.values == 298.15
            if unit != "K":
                assert temp_std.attrs["units"] == "K"

    def test_pressure_pa_to_hpa(self):
        """Test pressure conversion from Pa to hPa."""
        from skyborn.calc.GPI.xarray import _check_units

        pascal_units = ["Pa", "pa", "pascal", "Pascal", "N/m^2"]

        for unit in pascal_units:
            pres_pa = xr.DataArray(101325.0, attrs={"units": unit})
            pres_hpa, converted = _check_units(pres_pa, "pressure", ["hPa"])

            assert converted == True
            assert np.isclose(pres_hpa.values, 1013.25)
            assert pres_hpa.attrs["units"] == "hPa"

    def test_pressure_hpa_to_pa(self):
        """Test pressure conversion from hPa to Pa."""
        from skyborn.calc.GPI.xarray import _check_units

        hpa_units = ["hPa", "mb", "mbar", "millibar"]

        for unit in hpa_units:
            pres_hpa = xr.DataArray(1013.25, attrs={"units": unit})
            pres_pa, converted = _check_units(pres_hpa, "pressure", ["Pa"])

            assert converted == True
            assert np.isclose(pres_pa.values, 101325.0)
            assert pres_pa.attrs["units"] == "Pa"

    def test_pressure_kpa_conversion(self):
        """Test pressure conversion from kPa."""
        from skyborn.calc.GPI.xarray import _check_units

        # kPa to Pa
        pres_kpa = xr.DataArray(101.325, attrs={"units": "kPa"})
        pres_pa, converted = _check_units(pres_kpa, "pressure", ["Pa"])
        assert converted == True
        assert np.isclose(pres_pa.values, 101325.0)

        # kPa to hPa
        pres_kpa = xr.DataArray(101.325, attrs={"units": "kPa"})
        pres_hpa, converted = _check_units(pres_kpa, "pressure", ["hPa"])
        assert converted == True
        assert np.isclose(pres_hpa.values, 1013.25)

    def test_pressure_atm_conversion(self):
        """Test pressure conversion from atmospheres."""
        from skyborn.calc.GPI.xarray import _check_units

        # atm to Pa
        pres_atm = xr.DataArray(1.0, attrs={"units": "atm"})
        pres_pa, converted = _check_units(pres_atm, "pressure", ["Pa"])
        assert converted == True
        assert np.isclose(pres_pa.values, 101325.0)

        # atm to hPa
        pres_atm = xr.DataArray(1.0, attrs={"units": "atmosphere"})
        pres_hpa, converted = _check_units(pres_atm, "pressure", ["hPa"])
        assert converted == True
        assert np.isclose(pres_hpa.values, 1013.25)

    def test_mixing_ratio_g_per_kg(self):
        """Test mixing ratio conversion from g/kg to kg/kg."""
        from skyborn.calc.GPI.xarray import _check_units

        g_kg_units = ["g/kg", "g kg-1", "g kg^-1", "g.kg-1", "G/KG"]

        for unit in g_kg_units:
            mixr_g = xr.DataArray(10.0, attrs={"units": unit})
            mixr_kg, converted = _check_units(mixr_g, "mixing_ratio", ["kg/kg"])

            assert converted == True
            assert np.isclose(mixr_kg.values, 0.01)
            assert mixr_kg.attrs["units"] == "kg/kg"

    def test_mixing_ratio_mg_per_kg(self):
        """Test mixing ratio conversion from mg/kg to kg/kg."""
        from skyborn.calc.GPI.xarray import _check_units

        mixr_mg = xr.DataArray(10000.0, attrs={"units": "mg/kg"})
        mixr_kg, converted = _check_units(mixr_mg, "mixing_ratio", ["kg/kg"])

        assert converted == True
        assert np.isclose(mixr_kg.values, 0.01)
        assert mixr_kg.attrs["units"] == "kg/kg"

    def test_mixing_ratio_dimensionless(self):
        """Test mixing ratio with dimensionless units."""
        from skyborn.calc.GPI.xarray import _check_units

        dimensionless_units = ["dimensionless", "1", "fraction"]

        for unit in dimensionless_units:
            mixr = xr.DataArray(0.01, attrs={"units": unit})
            mixr_std, converted = _check_units(mixr, "mixing_ratio", ["kg/kg"])

            assert mixr_std.values == 0.01
            if unit != "kg/kg":
                assert mixr_std.attrs["units"] == "kg/kg"

    def test_mixing_ratio_percentage(self):
        """Test mixing ratio conversion from percentage."""
        from skyborn.calc.GPI.xarray import _check_units

        mixr_pct = xr.DataArray(1.0, attrs={"units": "%"})
        mixr_kg, converted = _check_units(mixr_pct, "mixing_ratio", ["kg/kg"])

        assert converted == True
        assert np.isclose(mixr_kg.values, 0.01)
        assert mixr_kg.attrs["units"] == "kg/kg"

    def test_mixing_ratio_ppmv(self):
        """Test mixing ratio conversion from ppmv."""
        from skyborn.calc.GPI.xarray import _check_units

        ppmv_units = ["ppmv", "ppm", "PPM"]

        for unit in ppmv_units:
            mixr_ppmv = xr.DataArray(16077.0, attrs={"units": unit})  # ~0.01 kg/kg
            mixr_kg, converted = _check_units(mixr_ppmv, "mixing_ratio", ["kg/kg"])

            assert converted == True
            assert np.isclose(mixr_kg.values, 0.01, rtol=0.01)
            assert mixr_kg.attrs["units"] == "kg/kg"

    def test_no_units_warning(self):
        """Test warning when no units attribute exists."""
        from skyborn.calc.GPI.xarray import _check_units

        data = xr.DataArray(298.15)  # No units attribute

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result, converted = _check_units(data, "temperature", ["K"])

            assert len(w) == 1
            assert "has no 'units' attribute" in str(w[0].message)
            assert converted == False

    def test_invalid_units_error(self):
        """Test error for unrecognized units."""
        from skyborn.calc.GPI.xarray import _check_units

        data = xr.DataArray(100.0, attrs={"units": "invalid_unit"})

        with pytest.raises(ValueError, match="has units 'invalid_unit' but expected"):
            _check_units(data, "temperature", ["K"])


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestDimensionDetection:
    """Test dimension detection functionality."""

    def test_detect_standard_dimensions(self):
        """Test detection of standard dimension names."""
        from skyborn.calc.GPI.xarray import _detect_atmospheric_dimensions

        # Test various dimension name combinations
        test_cases = [
            (["lon", "lat", "level"], (0, 1, 2, None)),
            (["longitude", "latitude", "pressure"], (0, 1, 2, None)),
            (["x", "y", "z"], (0, 1, 2, None)),
            (["time", "plev", "lat", "lon"], (3, 2, 1, 0)),
            (["t", "height", "nlat", "mlon"], (3, 2, 1, 0)),
        ]

        for dims, expected in test_cases:
            data = xr.DataArray(np.zeros([2] * len(dims)), dims=dims)
            result = _detect_atmospheric_dimensions(data)
            assert result == expected

    def test_detect_pressure_variations(self):
        """Test detection of various pressure level names."""
        from skyborn.calc.GPI.xarray import _detect_atmospheric_dimensions

        pressure_names = [
            "level",
            "lev",
            "plev",
            "pressure",
            "p",
            "height",
            "altitude",
            "isobaric",
            "model_level",
            "PRESSURE_LEVEL",
        ]

        for pname in pressure_names:
            data = xr.DataArray(np.zeros((3, 4, 5)), dims=[pname, "lat", "lon"])
            xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(data)
            assert levdim == 0

    def test_missing_level_dimension_error(self):
        """Test error when level dimension cannot be detected."""
        from skyborn.calc.GPI.xarray import _detect_atmospheric_dimensions

        data = xr.DataArray(np.zeros((10, 20)), dims=["unknown1", "unknown2"])

        with pytest.raises(ValueError, match="Could not auto-detect level dimension"):
            _detect_atmospheric_dimensions(data)

    def test_3d_missing_spatial_dims_error(self):
        """Test error for 3D data missing spatial dimensions."""
        from skyborn.calc.GPI.xarray import _detect_atmospheric_dimensions

        data = xr.DataArray(
            np.zeros((5, 10, 15)), dims=["level", "unknown1", "unknown2"]
        )

        with pytest.raises(ValueError, match="need both lat and lon dimensions"):
            _detect_atmospheric_dimensions(data)

    def test_4d_missing_time_error(self):
        """Test error for 4D data missing time dimension."""
        from skyborn.calc.GPI.xarray import _detect_atmospheric_dimensions

        data = xr.DataArray(
            np.zeros((5, 10, 15, 20)), dims=["unknown", "level", "lat", "lon"]
        )

        with pytest.raises(ValueError, match="need time, lat, and lon dimensions"):
            _detect_atmospheric_dimensions(data)


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class Test3DCalculations:
    """Test 3D gridded calculations."""

    def test_3d_calculation(self, sample_3d_data):
        """Test 3D potential intensity calculation."""
        from skyborn.calc.GPI.xarray import potential_intensity

        data = sample_3d_data

        # Create xarray objects
        levels = xr.DataArray(data["levels"], dims=["level"], attrs={"units": "mb"})

        temperature = xr.DataArray(
            data["temperature"],
            dims=["level", "lat", "lon"],
            coords={
                "level": levels,
                "lat": np.linspace(-90, 90, data["nlat"]),
                "lon": np.linspace(0, 360, data["nlon"]),
            },
            attrs={"units": "K"},
        )

        mixing_ratio = xr.DataArray(
            data["mixing_ratio"],
            dims=["level", "lat", "lon"],
            coords=temperature.coords,
            attrs={"units": "kg/kg"},
        )

        sst = xr.DataArray(
            data["sst"],
            dims=["lat", "lon"],
            coords={"lat": temperature.lat, "lon": temperature.lon},
            attrs={"units": "K"},
        )

        psl = xr.DataArray(
            data["psl"],
            dims=["lat", "lon"],
            coords={"lat": temperature.lat, "lon": temperature.lon},
            attrs={"units": "Pa"},
        )

        # Perform calculation
        result = potential_intensity(
            sst=sst,
            psl=psl,
            pressure_levels=levels,
            temperature=temperature,
            mixing_ratio=mixing_ratio,
        )

        # Check results
        assert isinstance(result, xr.Dataset)
        assert "min_pressure" in result.data_vars
        assert "pi" in result.data_vars
        assert result.min_pressure.shape == (data["nlat"], data["nlon"])
        assert result.pi.shape == (data["nlat"], data["nlon"])
        assert result.error_flag.values == 1  # Success

        # Check metadata
        assert (
            result.attrs["description"] == "3D gridded potential intensity calculation"
        )
        assert result.attrs["vertical_levels"] == data["nlevs"]

    def test_3d_with_unit_conversion(self, sample_3d_data):
        """Test 3D calculation with unit conversions."""
        from skyborn.calc.GPI.xarray import potential_intensity

        data = sample_3d_data

        # Create data with different units
        levels = xr.DataArray(
            data["levels"] * 100, dims=["level"], attrs={"units": "Pa"}  # Convert to Pa
        )

        temperature = xr.DataArray(
            data["temperature"] - 273.15,  # Convert to Celsius
            dims=["level", "lat", "lon"],
            coords={"level": data["levels"]},
            attrs={"units": "°C"},
        )

        mixing_ratio = xr.DataArray(
            data["mixing_ratio"] * 1000,  # Convert to g/kg
            dims=["level", "lat", "lon"],
            coords={"level": data["levels"]},
            attrs={"units": "g/kg"},
        )

        sst = xr.DataArray(data["sst"], dims=["lat", "lon"], attrs={"units": "K"})

        psl = xr.DataArray(
            data["psl"] / 100,  # Convert to hPa
            dims=["lat", "lon"],
            attrs={"units": "hPa"},
        )

        # Should handle conversions automatically
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Ignore conversion warnings
            result = potential_intensity(
                sst=sst,
                psl=psl,
                pressure_levels=levels,
                temperature=temperature,
                mixing_ratio=mixing_ratio,
            )

        assert isinstance(result, xr.Dataset)
        assert result.error_flag.values == 1

    def test_3d_dimension_validation(self, sample_3d_data):
        """Test dimension validation for 3D data."""
        from skyborn.calc.GPI.xarray import _potential_intensity_3d

        data = sample_3d_data

        # Create valid xarray objects
        levels = xr.DataArray(data["levels"], dims=["level"])
        temp = xr.DataArray(data["temperature"], dims=["level", "lat", "lon"])
        mixr = xr.DataArray(data["mixing_ratio"], dims=["level", "lat", "lon"])
        sst = xr.DataArray(data["sst"], dims=["lat", "lon"])
        psl = xr.DataArray(data["psl"], dims=["lat", "lon"])

        # Test with wrong dimensions for temperature
        temp_wrong = xr.DataArray(
            data["temperature"][:, :, 0], dims=["level", "lat"]  # 2D instead of 3D
        )

        with pytest.raises(ValueError, match="must be 3D arrays"):
            _potential_intensity_3d(sst, psl, levels, temp_wrong, mixr)

        # Test with wrong dimensions for SST
        sst_wrong = xr.DataArray(data["sst"][0, :], dims=["lat"])  # 1D instead of 2D

        with pytest.raises(ValueError, match="must be 2D arrays"):
            _potential_intensity_3d(sst_wrong, psl, levels, temp, mixr)


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class Test4DCalculations:
    """Test 4D time series calculations."""

    def test_4d_calculation(self, sample_4d_data):
        """Test 4D potential intensity calculation."""
        from skyborn.calc.GPI.xarray import potential_intensity

        data = sample_4d_data

        # Create xarray objects
        levels = xr.DataArray(data["levels"], dims=["level"], attrs={"units": "mb"})

        times = np.arange(data["ntimes"])
        lats = np.linspace(-90, 90, data["nlat"])
        lons = np.linspace(0, 360, data["nlon"])

        temperature = xr.DataArray(
            data["temperature"],
            dims=["time", "level", "lat", "lon"],
            coords={"time": times, "level": levels, "lat": lats, "lon": lons},
            attrs={"units": "K"},
        )

        mixing_ratio = xr.DataArray(
            data["mixing_ratio"],
            dims=["time", "level", "lat", "lon"],
            coords=temperature.coords,
            attrs={"units": "kg/kg"},
        )

        sst = xr.DataArray(
            data["sst"],
            dims=["time", "lat", "lon"],
            coords={"time": times, "lat": lats, "lon": lons},
            attrs={"units": "K"},
        )

        psl = xr.DataArray(
            data["psl"],
            dims=["time", "lat", "lon"],
            coords={"time": times, "lat": lats, "lon": lons},
            attrs={"units": "Pa"},
        )

        # Perform calculation
        result = potential_intensity(
            sst=sst,
            psl=psl,
            pressure_levels=levels,
            temperature=temperature,
            mixing_ratio=mixing_ratio,
        )

        # Check results
        assert isinstance(result, xr.Dataset)
        assert result.min_pressure.shape == (data["ntimes"], data["nlat"], data["nlon"])
        assert result.pi.shape == (data["ntimes"], data["nlat"], data["nlon"])
        assert result.error_flag.values == 1

        # Check coordinates
        assert "time" in result.coords
        assert "lat" in result.coords
        assert "lon" in result.coords

        # Check metadata
        assert (
            result.attrs["description"]
            == "4D time series potential intensity calculation"
        )

    def test_4d_dimension_ordering(self, sample_4d_data):
        """Test different dimension orderings for 4D data."""
        from skyborn.calc.GPI.xarray import potential_intensity

        data = sample_4d_data

        # Create data with different dimension order
        levels = xr.DataArray(data["levels"], dims=["pressure_level"])

        # Different order: (level, time, lat, lon) instead of (time, level, lat, lon)
        temp_reordered = np.transpose(data["temperature"], (1, 0, 2, 3))
        temperature = xr.DataArray(
            temp_reordered,
            dims=["pressure_level", "t", "latitude", "longitude"],
            attrs={"units": "K"},
        )

        mixr_reordered = np.transpose(data["mixing_ratio"], (1, 0, 2, 3))
        mixing_ratio = xr.DataArray(
            mixr_reordered,
            dims=["pressure_level", "t", "latitude", "longitude"],
            attrs={"units": "kg/kg"},
        )

        sst = xr.DataArray(
            data["sst"], dims=["t", "latitude", "longitude"], attrs={"units": "K"}
        )

        psl = xr.DataArray(
            data["psl"], dims=["t", "latitude", "longitude"], attrs={"units": "Pa"}
        )

        # Should handle different dimension names and orders
        result = potential_intensity(
            sst=sst,
            psl=psl,
            pressure_levels=levels,
            temperature=temperature,
            mixing_ratio=mixing_ratio,
        )

        assert isinstance(result, xr.Dataset)
        assert result.error_flag.values == 1
        assert "t" in result.dims
        assert "latitude" in result.dims
        assert "longitude" in result.dims

    def test_4d_validation_errors(self, sample_4d_data):
        """Test validation errors for 4D data."""
        from skyborn.calc.GPI.xarray import potential_intensity

        data = sample_4d_data

        # Create valid xarray objects
        levels = xr.DataArray(data["levels"], dims=["level"], attrs={"units": "mb"})
        temp = xr.DataArray(
            data["temperature"],
            dims=["time", "level", "lat", "lon"],
            attrs={"units": "K"},
        )
        mixr = xr.DataArray(
            data["mixing_ratio"],
            dims=["time", "level", "lat", "lon"],
            attrs={"units": "kg/kg"},
        )
        sst = xr.DataArray(
            data["sst"], dims=["time", "lat", "lon"], attrs={"units": "K"}
        )
        psl = xr.DataArray(
            data["psl"], dims=["time", "lat", "lon"], attrs={"units": "Pa"}
        )

        # Test with wrong number of dimensions for temperature (3D instead of 4D)
        temp_3d = xr.DataArray(
            data["temperature"][0],  # Remove time dimension
            dims=["level", "lat", "lon"],
            attrs={"units": "K"},
        )

        # This should raise an error because dimension mismatch
        with pytest.raises(ValueError, match="Unsupported number of dimensions"):
            # Using potential_intensity which routes to the appropriate function
            potential_intensity(sst[0], psl[0], levels, temp_3d, mixr[0])

        # Test with wrong SST dimensions (2D when expecting 3D for 4D data)
        sst_2d = xr.DataArray(
            data["sst"][0],  # Remove time dimension
            dims=["lat", "lon"],
            attrs={"units": "K"},
        )

        # This should work because it detects as 3D calculation based on temp dimensions
        # So we need to test differently - use 4D temp with 2D sst which should fail
        with pytest.raises(ValueError):
            potential_intensity(sst_2d, psl[0], levels, temp, mixr)


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestOutputDataset:
    """Test output dataset creation."""

    def test_create_output_dataset_profile(self):
        """Test output dataset creation for profile data."""
        from skyborn.calc.GPI.xarray import _create_output_dataset

        min_p = 920.5
        pi = 75.3
        error = 1

        result = _create_output_dataset(
            min_p,
            pi,
            error,
            sst_val=301.0,
            psl_val=101325.0,
            vertical_levels=8,
            vertical_dim="level",
            data_type="profile",
        )

        assert isinstance(result, xr.Dataset)
        assert float(result.min_pressure.values) == min_p
        assert float(result.pi.values) == pi
        assert int(result.error_flag.values) == error

        # Check attributes
        assert result.attrs["sst_input"] == 301.0
        assert result.attrs["psl_input"] == 101325.0
        assert result.attrs["vertical_levels"] == 8
        assert result.attrs["description"] == "profile potential intensity calculation"

    def test_create_output_dataset_3d(self):
        """Test output dataset creation for 3D data."""
        from skyborn.calc.GPI.xarray import _create_output_dataset

        nlat, nlon = 5, 6
        min_p = np.random.rand(nlat, nlon) * 100 + 900
        pi = np.random.rand(nlat, nlon) * 50 + 50
        error = 1

        coords = {"lat": np.linspace(-90, 90, nlat), "lon": np.linspace(0, 360, nlon)}

        result = _create_output_dataset(
            min_p,
            pi,
            error,
            input_coords=coords,
            vertical_levels=10,
            vertical_dim="pressure",
            data_type="3D",
        )

        assert isinstance(result, xr.Dataset)
        assert result.min_pressure.shape == (nlat, nlon)
        assert result.pi.shape == (nlat, nlon)
        assert "lat" in result.coords
        assert "lon" in result.coords

    def test_create_output_dataset_4d(self):
        """Test output dataset creation for 4D data."""
        from skyborn.calc.GPI.xarray import _create_output_dataset

        ntimes, nlat, nlon = 3, 4, 5
        min_p = np.random.rand(ntimes, nlat, nlon) * 100 + 900
        pi = np.random.rand(ntimes, nlat, nlon) * 50 + 50
        error = 1

        coords = {
            "time": np.arange(ntimes),
            "lat": np.linspace(-90, 90, nlat),
            "lon": np.linspace(0, 360, nlon),
        }

        result = _create_output_dataset(
            min_p,
            pi,
            error,
            input_coords=coords,
            vertical_levels=8,
            vertical_dim="level",
            data_type="4D time series",
        )

        assert isinstance(result, xr.Dataset)
        assert result.min_pressure.shape == (ntimes, nlat, nlon)
        assert result.pi.shape == (ntimes, nlat, nlon)
        assert "time" in result.coords
        assert (
            result.attrs["description"]
            == "4D time series potential intensity calculation"
        )


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestEdgeCases:
    """Test edge cases and special scenarios."""

    def test_sorted_pressure_levels(self):
        """Test handling of already sorted pressure levels."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Create data with pressure levels already sorted (high to low)
        levels = xr.DataArray(
            [1000, 850, 700, 500], dims=["level"], attrs={"units": "mb"}
        )
        temp = xr.DataArray([300, 290, 280, 270], dims=["level"], attrs={"units": "K"})
        mixr = xr.DataArray(
            [0.015, 0.010, 0.005, 0.001], dims=["level"], attrs={"units": "kg/kg"}
        )

        result = potential_intensity(
            sst=302.0,
            psl=101325.0,
            pressure_levels=levels,
            temperature=temp,
            mixing_ratio=mixr,
        )

        assert isinstance(result, xr.Dataset)

    def test_reverse_sorted_pressure_levels(self):
        """Test handling of reverse sorted pressure levels."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Create data with pressure levels sorted low to high (needs reversal)
        levels = xr.DataArray(
            [500, 700, 850, 1000], dims=["level"], attrs={"units": "mb"}
        )
        temp = xr.DataArray([270, 280, 290, 300], dims=["level"], attrs={"units": "K"})
        mixr = xr.DataArray(
            [0.001, 0.005, 0.010, 0.015], dims=["level"], attrs={"units": "kg/kg"}
        )

        result = potential_intensity(
            sst=302.0,
            psl=101325.0,
            pressure_levels=levels,
            temperature=temp,
            mixing_ratio=mixr,
        )

        assert isinstance(result, xr.Dataset)

    def test_nan_handling(self):
        """Test handling of NaN values in input."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Create data with some NaN values
        levels = xr.DataArray(
            [1000, 850, 700, 500], dims=["level"], attrs={"units": "mb"}
        )
        temp = xr.DataArray(
            [300, np.nan, 280, 270], dims=["level"], attrs={"units": "K"}
        )
        mixr = xr.DataArray(
            [0.015, 0.010, 0.005, 0.001], dims=["level"], attrs={"units": "kg/kg"}
        )

        # Should handle NaN appropriately (likely return error or NaN result)
        result = potential_intensity(
            sst=302.0,
            psl=101325.0,
            pressure_levels=levels,
            temperature=temp,
            mixing_ratio=mixr,
        )

        assert isinstance(result, xr.Dataset)
        # Check that error flag indicates a problem or result contains NaN
        assert result.error_flag.values != 1 or np.isnan(result.pi.values)

    def test_extreme_values(self):
        """Test calculation with extreme but valid values."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Very warm SST and low pressure (strong storm conditions)
        levels = xr.DataArray(
            [1000, 850, 700, 500], dims=["level"], attrs={"units": "mb"}
        )
        temp = xr.DataArray([305, 295, 285, 270], dims=["level"], attrs={"units": "K"})
        mixr = xr.DataArray(
            [0.020, 0.015, 0.008, 0.002], dims=["level"], attrs={"units": "kg/kg"}
        )

        result = potential_intensity(
            sst=305.0,  # Very warm SST
            psl=100000.0,  # Low pressure
            pressure_levels=levels,
            temperature=temp,
            mixing_ratio=mixr,
        )

        assert isinstance(result, xr.Dataset)
        # Should produce high potential intensity
        if result.error_flag.values == 1:
            assert result.pi.values > 0  # Should have positive intensity

    def test_minimal_vertical_levels(self):
        """Test with minimal number of vertical levels."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Only 2 levels (minimum for a profile)
        levels = xr.DataArray([1000, 500], dims=["level"], attrs={"units": "mb"})
        temp = xr.DataArray([300, 270], dims=["level"], attrs={"units": "K"})
        mixr = xr.DataArray([0.015, 0.002], dims=["level"], attrs={"units": "kg/kg"})

        result = potential_intensity(
            sst=301.0,
            psl=101325.0,
            pressure_levels=levels,
            temperature=temp,
            mixing_ratio=mixr,
        )

        assert isinstance(result, xr.Dataset)
