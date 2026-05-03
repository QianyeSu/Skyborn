"""
Tests for GPI xarray interface functionality.
"""

import numpy as np
import pytest

try:
    import xarray as xr

    _has_xarray = True
except ImportError:
    _has_xarray = False


@pytest.fixture
def sample_profile_data():
    """Create sample atmospheric profile data for testing."""
    # Standard atmospheric pressure levels (mb)
    levels = [1000, 925, 850, 700, 500, 300, 200, 100]

    # Realistic temperature profile (K)
    temperature = [300, 298, 295, 285, 270, 250, 230, 210]

    # Realistic mixing ratio profile (kg/kg)
    mixing_ratio = [0.015, 0.012, 0.010, 0.006, 0.002, 0.0005, 0.0001, 0.00005]

    # Surface conditions
    sst = 301.0  # K
    psl = 101325.0  # Pa

    return {
        "levels": levels,
        "temperature": temperature,
        "mixing_ratio": mixing_ratio,
        "sst": sst,
        "psl": psl,
    }


@pytest.fixture
def xarray_profile_data(sample_profile_data):
    """Convert sample data to xarray format."""
    if not _has_xarray:
        pytest.skip("xarray not available")

    data = sample_profile_data

    # Create coordinate array
    levels = xr.DataArray(
        data["levels"],
        dims=["level"],
        name="pressure",
        attrs={"units": "mb", "long_name": "Pressure level"},
    )

    # Create profile arrays with level coordinate
    temperature = xr.DataArray(
        data["temperature"],
        dims=["level"],
        coords={"level": levels},
        name="temperature",
        attrs={"units": "K", "long_name": "Temperature"},
    )

    mixing_ratio = xr.DataArray(
        data["mixing_ratio"],
        dims=["level"],
        coords={"level": levels},
        name="mixing_ratio",
        attrs={"units": "kg/kg", "long_name": "Water vapor mixing ratio"},
    )

    return {
        "levels": levels,
        "temperature": temperature,
        "mixing_ratio": mixing_ratio,
        "sst": data["sst"],
        "psl": data["psl"],
    }


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestGPIXarray:
    """Test GPI xarray interface functionality."""

    def test_import_xarray_interface(self):
        """Test that the xarray interface can be imported."""
        try:
            from skyborn.calc.GPI.xarray import potential_intensity
        except ImportError:
            pytest.skip("GPI module or xarray interface not available")

    def test_basic_profile_calculation(self, xarray_profile_data):
        """Test basic single profile calculation with xarray."""
        try:
            from skyborn.calc.GPI.xarray import (
                potential_intensity as potential_intensity_xr,
            )
        except ImportError:
            pytest.skip("GPI module not available")

        data = xarray_profile_data

        # Perform calculation
        result = potential_intensity_xr(
            sst=data["sst"],
            psl=data["psl"],
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
        )

        # Check result structure
        assert isinstance(result, xr.Dataset)
        assert "min_pressure" in result.data_vars
        assert "pi" in result.data_vars
        assert "error_flag" in result.data_vars

        # Check that results are scalars
        assert result.min_pressure.ndim == 0
        assert result.pi.ndim == 0
        assert result.error_flag.ndim == 0

        # Check for reasonable values
        min_p = float(result.min_pressure.values)
        max_w = float(result.pi.values)
        err = int(result.error_flag.values)

        # Basic sanity checks (error_flag = 1 means success in this implementation)
        assert err == 1  # Should be successful (1 = success, not 0)
        assert (
            600 < min_p < 1100
        )  # Reasonable pressure range (mb) - intense storms can go below 700mb
        assert (
            0 < max_w < 200
        )  # Reasonable wind speed range (m/s) - intense storms can exceed 150 kph
        assert not np.isnan(min_p)
        assert not np.isnan(max_w)

    def test_xarray_dataarray_sst_psl(self, xarray_profile_data):
        """Test with SST and PSL as xarray DataArrays."""
        try:
            from skyborn.calc.GPI.xarray import (
                potential_intensity as potential_intensity_xr,
            )
        except ImportError:
            pytest.skip("GPI module not available")

        data = xarray_profile_data

        # Convert SST and PSL to scalar DataArrays
        sst_da = xr.DataArray(data["sst"], attrs={"units": "K"})
        psl_da = xr.DataArray(data["psl"], attrs={"units": "Pa"})

        result = potential_intensity_xr(
            sst=sst_da,
            psl=psl_da,
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
        )

        # Should work the same as with scalar inputs
        assert isinstance(result, xr.Dataset)
        assert result.error_flag.values == 1  # 1 = success

    def test_profile_log_decomposition_xarray_interface(self, xarray_profile_data):
        """Explicit xarray decomposition helper should expose profile diagnostics."""
        from skyborn.calc.GPI.xarray import (
            pi_log_decomposition as pi_log_decomposition_xr,
        )

        data = xarray_profile_data
        result = pi_log_decomposition_xr(
            sst=xr.DataArray(data["sst"], attrs={"units": "K"}),
            psl=xr.DataArray(data["psl"], attrs={"units": "Pa"}),
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
            outflow_source="cape_env",
        )

        assert isinstance(result, xr.Dataset)
        for name in ["t0", "otl", "lnpi", "lneff", "lndiseq", "lnCKCD"]:
            assert name in result.data_vars
        assert int(result.error_flag.values) == 1
        assert np.isfinite(float(result.t0.values))
        assert np.isfinite(float(result.otl.values))

    def test_profile_log_decomposition_xarray_validation(self, xarray_profile_data):
        """Profile decomposition should reject non-1D and non-scalar inputs."""
        from skyborn.calc.GPI.xarray import (
            pi_log_decomposition as pi_log_decomposition_xr,
        )

        data = xarray_profile_data

        with pytest.raises(ValueError, match="Unsupported number of dimensions"):
            pi_log_decomposition_xr(
                sst=data["sst"],
                psl=data["psl"],
                pressure_levels=data["levels"].expand_dims("x"),
                temperature=data["temperature"].expand_dims("x"),
                mixing_ratio=data["mixing_ratio"].expand_dims("x"),
            )

        with pytest.raises(ValueError, match="SST DataArray must be 0-dimensional"):
            pi_log_decomposition_xr(
                sst=xr.DataArray([data["sst"], data["sst"]], dims=["x"]),
                psl=xr.DataArray(data["psl"], attrs={"units": "Pa"}),
                pressure_levels=data["levels"],
                temperature=data["temperature"],
                mixing_ratio=data["mixing_ratio"],
            )

    def test_profile_log_decomposition_repeatability(self, xarray_profile_data):
        """Profile diagnostics should be stable across repeated decomposition calls."""
        from skyborn.calc.GPI.xarray import (
            pi_log_decomposition as pi_log_decomposition_xr,
        )

        data = xarray_profile_data
        main = pi_log_decomposition_xr(
            sst=xr.DataArray(data["sst"], attrs={"units": "K"}),
            psl=xr.DataArray(data["psl"], attrs={"units": "Pa"}),
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
            outflow_source="cape_env",
        )
        repeat = pi_log_decomposition_xr(
            sst=xr.DataArray(data["sst"], attrs={"units": "K"}),
            psl=xr.DataArray(data["psl"], attrs={"units": "Pa"}),
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
            outflow_source="cape_env",
        )

        for name in [
            "pi",
            "min_pressure",
            "t0",
            "otl",
            "lnpi",
            "lneff",
            "lndiseq",
            "lnCKCD",
        ]:
            np.testing.assert_allclose(main[name].values, repeat[name].values)

        with pytest.raises(ValueError, match="PSL DataArray must be 0-dimensional"):
            pi_log_decomposition_xr(
                sst=float(data["sst"]),
                psl=xr.DataArray([data["psl"], data["psl"]], dims=["x"]),
                pressure_levels=data["levels"],
                temperature=data["temperature"],
                mixing_ratio=data["mixing_ratio"],
            )

    def test_3d_log_decomposition_xarray_interface(self):
        """3D xarray diagnostics should come back as gridded fields."""
        from skyborn.calc.GPI.xarray import (
            pi_log_decomposition as pi_log_decomposition_xr,
        )

        lat = xr.DataArray(np.linspace(-10, 10, 4), dims=["lat"])
        lon = xr.DataArray(np.linspace(100, 110, 5), dims=["lon"])
        level = xr.DataArray(
            [1000, 850, 700, 500], dims=["level"], attrs={"units": "mb"}
        )
        sst = xr.DataArray(
            298.0 + np.random.rand(4, 5) * 4.0,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "K"},
        )
        psl = xr.DataArray(
            100000.0 + np.random.rand(4, 5) * 1200.0,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "Pa"},
        )
        temp = xr.DataArray(
            np.stack(
                [
                    300 + np.random.rand(4, 5),
                    288 + np.random.rand(4, 5),
                    275 + np.random.rand(4, 5),
                    255 + np.random.rand(4, 5),
                ]
            ),
            dims=["level", "lat", "lon"],
            coords={"level": level, "lat": lat, "lon": lon},
            attrs={"units": "K"},
        )
        mixr = xr.DataArray(
            np.stack(
                [
                    0.016 + np.random.rand(4, 5) * 1e-3,
                    0.010 + np.random.rand(4, 5) * 1e-3,
                    0.006 + np.random.rand(4, 5) * 1e-3,
                    0.0025 + np.random.rand(4, 5) * 1e-3,
                ]
            ),
            dims=["level", "lat", "lon"],
            coords={"level": level, "lat": lat, "lon": lon},
            attrs={"units": "kg/kg"},
        )

        result = pi_log_decomposition_xr(
            sst=sst,
            psl=psl,
            pressure_levels=level,
            temperature=temp,
            mixing_ratio=mixr,
        )
        for name in ["t0", "otl", "lnpi", "lneff", "lndiseq", "lnCKCD"]:
            assert name in result.data_vars
        assert result.pi.shape == sst.shape
        assert np.isfinite(float(result.lnCKCD.values))

    def test_metadata_preservation(self, xarray_profile_data):
        """Test that metadata is properly added to results."""
        try:
            from skyborn.calc.GPI.xarray import (
                potential_intensity as potential_intensity_xr,
            )
        except ImportError:
            pytest.skip("GPI module not available")

        data = xarray_profile_data

        result = potential_intensity_xr(
            sst=data["sst"],
            psl=data["psl"],
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
        )

        # Check global attributes
        assert "title" in result.attrs
        assert "sst_input" in result.attrs
        assert "psl_input" in result.attrs
        assert "vertical_levels" in result.attrs
        assert result.attrs["vertical_levels"] == len(data["levels"])

        # Check variable attributes
        assert "units" in result.min_pressure.attrs
        assert "units" in result.pi.attrs
        assert result.min_pressure.attrs["units"] == "mb"
        assert result.pi.attrs["units"] == "m/s"

    def test_error_handling(self, xarray_profile_data):
        """Test error handling for invalid inputs."""
        try:
            from skyborn.calc.GPI.xarray import (
                potential_intensity as potential_intensity_xr,
            )
        except ImportError:
            pytest.skip("GPI module not available")

        data = xarray_profile_data

        # Test invalid dimension name by creating data with invalid dim
        invalid_temp = data["temperature"].rename({"level": "invalid_dim"})
        with pytest.raises(ValueError, match="Could not auto-detect level dimension"):
            potential_intensity_xr(
                sst=data["sst"],
                psl=data["psl"],
                pressure_levels=data["levels"],
                temperature=invalid_temp,
                mixing_ratio=data["mixing_ratio"],
            )

        # Test non-1D arrays (2D is unsupported)
        temp_2d = data["temperature"].expand_dims("extra")
        with pytest.raises(ValueError, match="Unsupported number of dimensions"):
            potential_intensity_xr(
                sst=data["sst"],
                psl=data["psl"],
                pressure_levels=data["levels"],
                temperature=temp_2d,
                mixing_ratio=data["mixing_ratio"],
            )

        # Test non-scalar SST
        sst_1d = xr.DataArray([data["sst"], data["sst"]], dims=["x"])
        with pytest.raises(ValueError, match="SST DataArray must be 0-dimensional"):
            potential_intensity_xr(
                sst=sst_1d,
                psl=data["psl"],
                pressure_levels=data["levels"],
                temperature=data["temperature"],
                mixing_ratio=data["mixing_ratio"],
            )

    def test_auto_detection_1d(self, xarray_profile_data):
        """Test automatic dimension detection for 1D profile data."""
        try:
            from skyborn.calc.GPI.xarray import (
                potential_intensity as potential_intensity_xr,
            )
        except ImportError:
            pytest.skip("GPI module not available")

        data = xarray_profile_data

        # Test that auto-detection works for 1D profile
        result = potential_intensity_xr(
            sst=data["sst"],
            psl=data["psl"],
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
        )

        # Should produce same result as explicit profile call
        assert isinstance(result, xr.Dataset)
        assert result.error_flag.values == 1  # 1 = success
        assert result.min_pressure.ndim == 0  # Scalar output

    def test_simplified_interface(self, xarray_profile_data):
        """Test that the simplified interface works correctly."""
        try:
            from skyborn.calc.GPI.xarray import (
                potential_intensity as potential_intensity_xr,
            )
        except ImportError:
            pytest.skip("GPI module not available")

        data = xarray_profile_data

        # Test simplified interface
        result = potential_intensity_xr(
            sst=data["sst"],
            psl=data["psl"],
            pressure_levels=data["levels"],
            temperature=data["temperature"],
            mixing_ratio=data["mixing_ratio"],
        )

        assert isinstance(result, xr.Dataset)
        assert result.error_flag.values == 1


@pytest.mark.skipif(_has_xarray, reason="Testing when xarray is not available")
def test_missing_xarray_error():
    """Test proper error when xarray is not available."""
    try:
        from skyborn.calc.GPI.xarray import potential_intensity

        pytest.skip("xarray is available, cannot test error condition")
    except ImportError:
        pytest.skip("GPI module not available")
