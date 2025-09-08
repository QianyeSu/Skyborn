"""
Additional tests for GPI xarray interface to achieve complete coverage.
Focus on edge cases and uncovered unit conversion scenarios.
"""

import warnings

import numpy as np
import pytest

try:
    import xarray as xr

    _has_xarray = True
except ImportError:
    _has_xarray = False


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestAdditionalUnitConversions:
    """Test additional unit conversion scenarios for complete coverage."""

    def test_pressure_pa_standardization(self):
        """Test Pascal unit standardization when already in Pa but with different format."""
        from skyborn.calc.GPI.xarray import _check_units

        # Test Pa to Pa standardization (line 339-345)
        pres = xr.DataArray(101325.0, attrs={"units": "pascal"})
        pres_std, converted = _check_units(pres, "pressure", ["Pa"])

        # Should standardize the unit string but not convert value
        assert converted == False
        assert pres_std.values == 101325.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pres = xr.DataArray(101325.0, attrs={"units": "Pascal"})
            pres_std, converted = _check_units(pres, "pressure", ["Pa"])
            assert len(w) == 1
            assert "Standardized pressure units from Pascal to Pa" in str(w[0].message)

    def test_pressure_hpa_mb_standardization(self):
        """Test hPa/mb standardization when target is hPa or mb."""
        from skyborn.calc.GPI.xarray import _check_units

        # Test standardization from various hPa formats to standard hPa (lines 375-381)
        hpa_variations = ["hectopascal", "hpascal", "MILLIBAR", "mbars"]

        for unit in hpa_variations:
            pres = xr.DataArray(1013.25, attrs={"units": unit})
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                pres_std, converted = _check_units(pres, "pressure", ["hPa"])

                # Should standardize but not convert
                assert pres_std.values == 1013.25
                if unit not in ["hPa", "mb"]:
                    assert len(w) == 1
                    assert f"Standardized pressure units from {unit} to hPa" in str(
                        w[0].message
                    )
                    assert converted == False

    def test_pressure_mb_to_mb_no_conversion(self):
        """Test that mb to mb doesn't trigger conversion."""
        from skyborn.calc.GPI.xarray import _check_units

        # mb is already acceptable for expected mb
        pres = xr.DataArray(1013.25, attrs={"units": "mb"})
        pres_result, converted = _check_units(pres, "pressure", ["mb"])

        assert converted == False
        assert pres_result.values == 1013.25
        assert pres_result.attrs["units"] == "mb"

    def test_pressure_hpa_to_hpa_no_conversion(self):
        """Test that hPa to hPa doesn't trigger conversion."""
        from skyborn.calc.GPI.xarray import _check_units

        # hPa is already acceptable for expected hPa
        pres = xr.DataArray(1013.25, attrs={"units": "hPa"})
        pres_result, converted = _check_units(pres, "pressure", ["hPa"])

        assert converted == False
        assert pres_result.values == 1013.25
        assert pres_result.attrs["units"] == "hPa"


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestProfileValidationEdgeCases:
    """Test edge cases in profile validation."""

    def test_profile_1d_validation_fails(self):
        """Test that profile validation correctly identifies dimension mismatches."""
        from skyborn.calc.GPI.xarray import _potential_intensity_profile

        # Create 2D arrays that should fail 1D validation (line 700)
        levels = xr.DataArray(
            [[1000, 850], [1000, 850]], dims=["x", "level"], attrs={"units": "mb"}
        )
        temp = xr.DataArray(
            [[300, 290], [301, 291]], dims=["x", "level"], attrs={"units": "K"}
        )
        mixr = xr.DataArray(
            [[0.015, 0.010], [0.016, 0.011]],
            dims=["x", "level"],
            attrs={"units": "kg/kg"},
        )

        with pytest.raises(ValueError, match="Only 1D profiles are supported"):
            _potential_intensity_profile(
                sst=302.0,
                psl=101325.0,
                pressure_levels=levels,
                temperature=temp,
                mixing_ratio=mixr,
            )

    def test_profile_scalar_extraction_fails(self):
        """Test scalar extraction with non-scalar DataArray."""
        from skyborn.calc.GPI.xarray import _potential_intensity_profile

        # Valid 1D arrays
        levels = xr.DataArray([1000, 850, 700], dims=["level"], attrs={"units": "mb"})
        temp = xr.DataArray([300, 290, 280], dims=["level"], attrs={"units": "K"})
        mixr = xr.DataArray(
            [0.015, 0.010, 0.005], dims=["level"], attrs={"units": "kg/kg"}
        )

        # Non-scalar SST (lines 707-710)
        sst_array = xr.DataArray([301, 302], dims=["x"], attrs={"units": "K"})

        with pytest.raises(ValueError, match="SST DataArray must be 0-dimensional"):
            _potential_intensity_profile(
                sst=sst_array,
                psl=101325.0,
                pressure_levels=levels,
                temperature=temp,
                mixing_ratio=mixr,
            )

        # Non-scalar PSL
        psl_array = xr.DataArray([101325, 101300], dims=["x"], attrs={"units": "Pa"})

        with pytest.raises(ValueError, match="PSL DataArray must be 0-dimensional"):
            _potential_intensity_profile(
                sst=302.0,
                psl=psl_array,
                pressure_levels=levels,
                temperature=temp,
                mixing_ratio=mixr,
            )


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class Test3DValidationEdgeCases:
    """Test edge cases in 3D data validation."""

    def test_3d_incorrect_spatial_dims(self):
        """Test 3D calculation with incorrect number of spatial dimensions."""
        from skyborn.calc.GPI.xarray import _potential_intensity_3d

        # Create data with only 1 spatial dimension (should have 2)
        levels = xr.DataArray([1000, 850, 700], dims=["level"], attrs={"units": "mb"})

        # 2D data instead of 3D (level + 1 spatial instead of level + 2 spatial)
        temp = xr.DataArray(
            np.random.rand(3, 10), dims=["level", "x"], attrs={"units": "K"}
        )
        mixr = xr.DataArray(
            np.random.rand(3, 10), dims=["level", "x"], attrs={"units": "kg/kg"}
        )

        # SST/PSL with wrong dimensions
        sst = xr.DataArray(np.random.rand(10), dims=["x"], attrs={"units": "K"})
        psl = xr.DataArray(np.random.rand(10), dims=["x"], attrs={"units": "Pa"})

        # This should fail because not enough spatial dimensions (line 791)
        with pytest.raises(ValueError, match="Expected exactly 2 spatial dimensions"):
            _potential_intensity_3d(sst, psl, levels, temp, mixr)


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class Test4DSpecificEdgeCases:
    """Test specific edge cases in 4D calculations."""

    def test_4d_incorrect_output_dims(self):
        """Test 4D calculation with incorrect dimension setup."""
        from skyborn.calc.GPI.xarray import _potential_intensity_4d

        # Create data with wrong number of output dimensions (line 874-875)
        levels = xr.DataArray([1000, 850], dims=["level"], attrs={"units": "mb"})

        # Create 3D data (missing time) but try to use as 4D
        temp_3d = xr.DataArray(
            np.random.rand(2, 5, 6), dims=["level", "lat", "lon"], attrs={"units": "K"}
        )
        mixr_3d = xr.DataArray(
            np.random.rand(2, 5, 6),
            dims=["level", "lat", "lon"],
            attrs={"units": "kg/kg"},
        )
        sst_2d = xr.DataArray(
            np.random.rand(5, 6), dims=["lat", "lon"], attrs={"units": "K"}
        )
        psl_2d = xr.DataArray(
            np.random.rand(5, 6), dims=["lat", "lon"], attrs={"units": "Pa"}
        )

        # This should fail validation (line 864, 869)
        with pytest.raises(ValueError, match="must be 4D arrays"):
            _potential_intensity_4d(sst_2d, psl_2d, levels, temp_3d, mixr_3d)

    def test_4d_mismatched_sst_psl_dims(self):
        """Test 4D calculation with mismatched SST/PSL dimensions."""
        from skyborn.calc.GPI.xarray import _potential_intensity_4d

        # Create valid 4D temperature and mixing ratio
        temp = xr.DataArray(
            np.random.rand(3, 2, 5, 6),
            dims=["time", "level", "lat", "lon"],
            attrs={"units": "K"},
        )
        mixr = xr.DataArray(
            np.random.rand(3, 2, 5, 6),
            dims=["time", "level", "lat", "lon"],
            attrs={"units": "kg/kg"},
        )
        levels = xr.DataArray([1000, 850], dims=["level"], attrs={"units": "mb"})

        # SST with wrong dimensions (2D instead of 3D)
        sst_wrong = xr.DataArray(
            np.random.rand(5, 6), dims=["lat", "lon"], attrs={"units": "K"}
        )
        psl_correct = xr.DataArray(
            np.random.rand(3, 5, 6), dims=["time", "lat", "lon"], attrs={"units": "Pa"}
        )

        # Should fail because SST doesn't have time dimension (line 869)
        with pytest.raises(ValueError, match="must be 3D arrays"):
            _potential_intensity_4d(sst_wrong, psl_correct, levels, temp, mixr)


@pytest.mark.skipif(not _has_xarray, reason="xarray not available")
class TestMainInterfaceEdgeCases:
    """Test edge cases in the main potential_intensity function."""

    def test_unsupported_2d_dimensions(self):
        """Test with 2D data which is not supported."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Create 2D temperature data (no vertical dimension)
        temp_2d = xr.DataArray(
            np.random.rand(10, 12), dims=["lat", "lon"], attrs={"units": "K"}
        )
        mixr_2d = xr.DataArray(
            np.random.rand(10, 12), dims=["lat", "lon"], attrs={"units": "kg/kg"}
        )
        levels = xr.DataArray([1000], dims=["level"], attrs={"units": "mb"})
        sst = xr.DataArray(
            np.random.rand(10, 12), dims=["lat", "lon"], attrs={"units": "K"}
        )
        psl = xr.DataArray(
            np.random.rand(10, 12), dims=["lat", "lon"], attrs={"units": "Pa"}
        )

        # Should fail because 2D is not supported (line 967)
        with pytest.raises(ValueError, match="Could not auto-detect level dimension"):
            potential_intensity(sst, psl, levels, temp_2d, mixr_2d)

    def test_unsupported_5d_dimensions(self):
        """Test with 5D data which is not supported."""
        from skyborn.calc.GPI.xarray import potential_intensity

        # Create 5D temperature data
        temp_5d = xr.DataArray(
            np.random.rand(2, 3, 4, 5, 6),
            dims=["time", "level", "lat", "lon", "extra"],
            attrs={"units": "K"},
        )
        mixr_5d = xr.DataArray(
            np.random.rand(2, 3, 4, 5, 6),
            dims=["time", "level", "lat", "lon", "extra"],
            attrs={"units": "kg/kg"},
        )
        levels = xr.DataArray([1000, 850, 700], dims=["level"], attrs={"units": "mb"})
        sst = xr.DataArray(
            np.random.rand(2, 4, 5, 6),
            dims=["time", "lat", "lon", "extra"],
            attrs={"units": "K"},
        )
        psl = xr.DataArray(
            np.random.rand(2, 4, 5, 6),
            dims=["time", "lat", "lon", "extra"],
            attrs={"units": "Pa"},
        )

        # Should fail because 5D is not supported (line 967-969)
        with pytest.raises(ValueError, match="Unsupported number of dimensions"):
            potential_intensity(sst, psl, levels, temp_5d, mixr_5d)
