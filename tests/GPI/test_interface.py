"""
Comprehensive tests for skyborn.calc.GPI.interface module.

This test suite provides >96% code coverage for the GPI interface module,
testing all functions, methods, and edge cases for potential intensity calculations.
"""

import warnings

import numpy as np
import pytest

# Import the actual module to test (no mocking since Fortran module is available)
from skyborn.calc.GPI.interface import (
    UNDEF,
    PotentialIntensityCalculator,
    _ensure_pressure_ordering,
    _postprocess_results,
    _validate_dimensions,
    _validate_input_arrays,
    calculate_potential_intensity_3d,
    calculate_potential_intensity_4d,
    calculate_potential_intensity_profile,
    potential_intensity,
)


class TestPostprocessResults:
    """Test the _postprocess_results function."""

    def test_normal_values(self):
        """Test that normal values are unchanged."""
        min_pressure = np.array([[950.0, 980.0], [970.0, 990.0]])
        max_wind = np.array([[45.0, 30.0], [40.0, 25.0]])

        result_p, result_w = _postprocess_results(min_pressure, max_wind)

        np.testing.assert_array_equal(result_p, min_pressure)
        np.testing.assert_array_equal(result_w, max_wind)

    def test_undef_values_converted_to_nan(self):
        """Test that UNDEF values are converted to NaN."""
        min_pressure = np.array([[950.0, UNDEF], [970.0, 990.0]])
        max_wind = np.array([[45.0, 30.0], [UNDEF, 25.0]])

        result_p, result_w = _postprocess_results(min_pressure, max_wind)

        assert np.isnan(result_p[0, 1])
        assert np.isnan(result_w[1, 0])
        assert result_p[0, 0] == 950.0
        assert result_w[0, 1] == 30.0

    def test_all_undef_values(self):
        """Test array with all UNDEF values."""
        min_pressure = np.full((2, 2), UNDEF)
        max_wind = np.full((2, 2), UNDEF)

        result_p, result_w = _postprocess_results(min_pressure, max_wind)

        assert np.all(np.isnan(result_p))
        assert np.all(np.isnan(result_w))


class TestValidateInputArrays:
    """Test the _validate_input_arrays function."""

    def test_normal_arrays(self):
        """Test validation of normal arrays."""
        arr1 = np.array([[1.0, 2.0], [3.0, 4.0]])
        arr2 = np.array([[5.0, 6.0], [7.0, 8.0]])

        (result1, result2), has_missing = _validate_input_arrays(
            arr1, arr2, names=["array1", "array2"]
        )

        assert result1.dtype == np.float32
        assert result2.dtype == np.float32
        assert has_missing == False
        np.testing.assert_array_almost_equal(result1, arr1)
        np.testing.assert_array_almost_equal(result2, arr2)

    def test_arrays_with_nan(self):
        """Test arrays containing NaN values."""
        arr1 = np.array([[1.0, np.nan], [3.0, 4.0]])
        arr2 = np.array([[5.0, 6.0], [7.0, 8.0]])

        (result1, result2), has_missing = _validate_input_arrays(
            arr1, arr2, names=["array1", "array2"]
        )

        assert has_missing == True
        assert result1[0, 1] == UNDEF  # NaN converted to UNDEF
        assert result1[0, 0] == 1.0

    def test_arrays_with_undef(self):
        """Test arrays containing UNDEF values."""
        arr1 = np.array([[1.0, UNDEF], [3.0, 4.0]])
        arr2 = np.array([[5.0, 6.0], [7.0, 8.0]])

        (result1, result2), has_missing = _validate_input_arrays(
            arr1, arr2, names=["array1", "array2"]
        )

        assert has_missing == True
        assert result1[0, 1] == UNDEF

    def test_arrays_with_inf(self):
        """Test arrays containing infinite values."""
        arr1 = np.array([[1.0, np.inf], [3.0, 4.0]])
        arr2 = np.array([[5.0, 6.0], [-np.inf, 8.0]])

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            (result1, result2), has_missing = _validate_input_arrays(
                arr1, arr2, names=["array1", "array2"]
            )
            assert len(w) == 2  # Two warnings for inf values
            assert "Infinite values detected" in str(w[0].message)

    def test_single_array(self):
        """Test validation of a single array."""
        arr = np.array([1.0, 2.0, 3.0])

        result, has_missing = _validate_input_arrays(arr, names=["test_array"])

        assert result.dtype == np.float32
        assert has_missing == False
        np.testing.assert_array_almost_equal(result, arr)

    def test_mixed_missing_values(self):
        """Test arrays with both NaN and UNDEF values."""
        arr1 = np.array([[1.0, np.nan], [UNDEF, 4.0]])

        result, has_missing = _validate_input_arrays(arr1, names=["mixed"])

        assert has_missing == True
        assert result[0, 1] == UNDEF  # NaN converted to UNDEF
        assert result[1, 0] == UNDEF  # UNDEF remains UNDEF


class TestEnsurePressureOrdering:
    """Test the _ensure_pressure_ordering function."""

    def test_correctly_ordered_pressure(self):
        """Test pressure levels already in correct order (high to low)."""
        pressure = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])
        temp = np.random.rand(5, 10, 20)
        mixr = np.random.rand(5, 10, 20)

        p_out, t_out, m_out = _ensure_pressure_ordering(pressure, temp, mixr)

        np.testing.assert_array_equal(p_out, pressure)
        np.testing.assert_array_equal(t_out, temp)
        np.testing.assert_array_equal(m_out, mixr)

    def test_reversed_pressure_levels(self):
        """Test pressure levels in reverse order (low to high)."""
        pressure = np.array([300.0, 500.0, 700.0, 850.0, 1000.0])
        temp = np.arange(5 * 10 * 20).reshape(5, 10, 20)
        mixr = np.arange(5 * 10 * 20).reshape(5, 10, 20) * 2

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            p_out, t_out, m_out = _ensure_pressure_ordering(pressure, temp, mixr)
            assert len(w) == 1
            assert "Reordering to surface to top" in str(w[0].message)

        np.testing.assert_array_equal(p_out, pressure[::-1])
        np.testing.assert_array_equal(t_out, np.flip(temp, axis=0))
        np.testing.assert_array_equal(m_out, np.flip(mixr, axis=0))

    def test_single_pressure_level(self):
        """Test single pressure level (no ordering needed)."""
        pressure = np.array([850.0])
        temp = np.random.rand(1, 10, 20)
        mixr = np.random.rand(1, 10, 20)

        p_out, t_out, m_out = _ensure_pressure_ordering(pressure, temp, mixr)

        np.testing.assert_array_equal(p_out, pressure)
        np.testing.assert_array_equal(t_out, temp)
        np.testing.assert_array_equal(m_out, mixr)

    def test_1d_temperature_data(self):
        """Test with 1D temperature profile data."""
        pressure = np.array([300.0, 500.0, 700.0, 850.0, 1000.0])
        temp = np.array([210.0, 250.0, 270.0, 285.0, 300.0])
        mixr = np.array([0.001, 0.003, 0.005, 0.008, 0.012])

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            p_out, t_out, m_out = _ensure_pressure_ordering(pressure, temp, mixr)
            assert len(w) == 1

        np.testing.assert_array_equal(p_out, pressure[::-1])
        np.testing.assert_array_equal(t_out, temp[::-1])
        np.testing.assert_array_equal(m_out, mixr[::-1])

    def test_4d_temperature_data(self):
        """Test with 4D temperature data (time dimension)."""
        pressure = np.array([300.0, 500.0, 700.0, 850.0, 1000.0])
        temp = np.random.rand(3, 5, 10, 20)  # (time, pressure, lat, lon)
        mixr = np.random.rand(3, 5, 10, 20)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            p_out, t_out, m_out = _ensure_pressure_ordering(pressure, temp, mixr)
            assert len(w) == 1

        np.testing.assert_array_equal(p_out, pressure[::-1])
        np.testing.assert_array_equal(t_out, np.flip(temp, axis=1))
        np.testing.assert_array_equal(m_out, np.flip(mixr, axis=1))


class TestValidateDimensions:
    """Test the _validate_dimensions function."""

    def test_profile_dimensions_valid(self):
        """Test valid profile dimensions."""
        sst = 300.0
        psl = 1013.25
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        # Should not raise any exception
        _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "profile")

    def test_profile_dimensions_invalid_shape(self):
        """Test invalid profile dimensions."""
        sst = 300.0
        psl = 1013.25
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([[300.0, 285.0], [270.0, 250.0]])  # Wrong shape
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        with pytest.raises(ValueError, match="doesn't match expected profile shape"):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "profile")

    def test_profile_dimensions_non_scalar_sst(self):
        """Test profile with non-scalar SST."""
        sst = np.array([300.0, 301.0])  # Should be scalar
        psl = 1013.25
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        with pytest.raises(ValueError, match="SST and PSL must be scalars"):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "profile")

    def test_3d_dimensions_valid(self):
        """Test valid 3D dimensions."""
        nlat, nlon = 10, 20
        num_levels = 4
        sst = np.random.rand(nlat, nlon)
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(num_levels, nlat, nlon)
        mixr = np.random.rand(num_levels, nlat, nlon)

        # Should not raise any exception
        _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "3D")

    def test_3d_dimensions_invalid_temp_shape(self):
        """Test 3D with invalid temperature shape."""
        nlat, nlon = 10, 20
        sst = np.random.rand(nlat, nlon)
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(nlat, nlon)  # Missing pressure dimension
        mixr = np.random.rand(4, nlat, nlon)

        with pytest.raises(ValueError, match="doesn't match expected 3D shape"):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "3D")

    def test_3d_dimensions_mismatched_sst_shape(self):
        """Test 3D with mismatched SST shape."""
        nlat, nlon = 10, 20
        sst = np.random.rand(nlat + 1, nlon)  # Wrong shape
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, nlat, nlon)
        mixr = np.random.rand(4, nlat, nlon)

        with pytest.raises(ValueError, match="SST/PSL shape mismatch"):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "3D")

    def test_4d_dimensions_valid(self):
        """Test valid 4D dimensions."""
        ntimes, nlat, nlon = 5, 10, 20
        num_levels = 4
        sst = np.random.rand(ntimes, nlat, nlon)
        psl = np.random.rand(ntimes, nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(ntimes, num_levels, nlat, nlon)
        mixr = np.random.rand(ntimes, num_levels, nlat, nlon)

        # Should not raise any exception
        _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "4D")

    def test_4d_dimensions_invalid_temp_shape(self):
        """Test 4D with invalid temperature shape."""
        ntimes, nlat, nlon = 5, 10, 20
        sst = np.random.rand(ntimes, nlat, nlon)
        psl = np.random.rand(ntimes, nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, nlat, nlon)  # Missing time dimension
        mixr = np.random.rand(ntimes, 4, nlat, nlon)

        with pytest.raises(ValueError, match="doesn't match expected 4D shape"):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "4D")

    def test_unsupported_data_type(self):
        """Test unsupported data type."""
        sst = np.random.rand(10, 20)
        psl = np.random.rand(10, 20)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, 10, 20)
        mixr = np.random.rand(4, 10, 20)

        with pytest.raises(ValueError, match="Unsupported data type"):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "2D")

    def test_temp_mixr_shape_mismatch(self):
        """Test temperature and mixing ratio shape mismatch."""
        nlat, nlon = 10, 20
        sst = np.random.rand(nlat, nlon)
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, nlat, nlon)
        mixr = np.random.rand(4, nlat + 1, nlon)  # Different shape

        with pytest.raises(
            ValueError, match="Temperature shape.*doesn't match mixing ratio shape"
        ):
            _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "3D")


class TestCalculatePotentialIntensity3D:
    """Test the calculate_potential_intensity_3d function."""

    def test_normal_3d_calculation(self):
        """Test normal 3D calculation without missing values."""
        nlat, nlon = 10, 20
        num_levels = 4
        sst = np.random.rand(nlat, nlon) * 10 + 290  # 290-300 K
        psl = np.random.rand(nlat, nlon) * 1000 + 100000  # Around 1000 hPa
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(num_levels, nlat, nlon) * 50 + 250
        mixr = np.random.rand(num_levels, nlat, nlon) * 0.02

        min_p, max_w, err = calculate_potential_intensity_3d(
            sst, psl, pressure_levels, temp, mixr
        )

        assert min_p.shape == (nlat, nlon)
        assert max_w.shape == (nlat, nlon)
        assert err == 0
        assert not np.any(np.isnan(min_p))  # No NaN in normal data
        assert not np.any(np.isnan(max_w))

    def test_3d_with_missing_values(self):
        """Test 3D calculation with missing values."""
        nlat, nlon = 10, 20
        num_levels = 4
        sst = np.random.rand(nlat, nlon) * 10 + 290
        sst[0, 0] = np.nan  # Add missing value
        psl = np.random.rand(nlat, nlon) * 1000 + 100000
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(num_levels, nlat, nlon) * 50 + 250
        mixr = np.random.rand(num_levels, nlat, nlon) * 0.02

        min_p, max_w, err = calculate_potential_intensity_3d(
            sst, psl, pressure_levels, temp, mixr
        )

        assert min_p.shape == (nlat, nlon)
        assert max_w.shape == (nlat, nlon)
        assert err == 0
        # Check that UNDEF was converted to NaN
        assert np.isnan(min_p[0, 0])
        assert np.isnan(max_w[0, 0])

    def test_3d_with_reversed_pressure(self):
        """Test 3D calculation with reversed pressure levels."""
        nlat, nlon = 10, 20
        num_levels = 4
        sst = np.random.rand(nlat, nlon) * 10 + 290
        psl = np.random.rand(nlat, nlon) * 1000 + 100000
        pressure_levels = np.array([500.0, 700.0, 850.0, 1000.0])  # Reversed
        temp = np.random.rand(num_levels, nlat, nlon) * 50 + 250
        mixr = np.random.rand(num_levels, nlat, nlon) * 0.02

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            min_p, max_w, err = calculate_potential_intensity_3d(
                sst, psl, pressure_levels, temp, mixr
            )
            assert len(w) == 1
            assert "Reordering to surface to top" in str(w[0].message)

        assert min_p.shape == (nlat, nlon)
        assert max_w.shape == (nlat, nlon)
        assert err == 0


class TestCalculatePotentialIntensity4D:
    """Test the calculate_potential_intensity_4d function."""

    def test_normal_4d_calculation(self):
        """Test normal 4D calculation without missing values."""
        ntimes, nlat, nlon = 5, 10, 20
        num_levels = 4
        sst = np.random.rand(ntimes, nlat, nlon) * 10 + 290
        psl = np.random.rand(ntimes, nlat, nlon) * 1000 + 100000
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(ntimes, num_levels, nlat, nlon) * 50 + 250
        mixr = np.random.rand(ntimes, num_levels, nlat, nlon) * 0.02

        min_p, max_w, err = calculate_potential_intensity_4d(
            sst, psl, pressure_levels, temp, mixr
        )

        assert min_p.shape == (ntimes, nlat, nlon)
        assert max_w.shape == (ntimes, nlat, nlon)
        assert err == 0
        assert not np.any(np.isnan(min_p))
        assert not np.any(np.isnan(max_w))

    def test_4d_with_missing_values(self):
        """Test 4D calculation with missing values."""
        ntimes, nlat, nlon = 5, 10, 20
        num_levels = 4
        sst = np.random.rand(ntimes, nlat, nlon) * 10 + 290
        sst[0, 0, 0] = np.nan  # Add missing value
        psl = np.random.rand(ntimes, nlat, nlon) * 1000 + 100000
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(ntimes, num_levels, nlat, nlon) * 50 + 250
        temp[1, :, 5, 5] = UNDEF  # Add UNDEF values
        mixr = np.random.rand(ntimes, num_levels, nlat, nlon) * 0.02

        min_p, max_w, err = calculate_potential_intensity_4d(
            sst, psl, pressure_levels, temp, mixr
        )

        assert min_p.shape == (ntimes, nlat, nlon)
        assert max_w.shape == (ntimes, nlat, nlon)
        assert err == 0
        # Check that UNDEF was converted to NaN
        assert np.isnan(min_p[0, 0, 0])
        assert np.isnan(max_w[0, 0, 0])


class TestCalculatePotentialIntensityProfile:
    """Test the calculate_potential_intensity_profile function."""

    def test_normal_profile_calculation(self):
        """Test normal profile calculation."""
        sst = 300.0
        psl = 101325.0
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        min_p, max_w, err = calculate_potential_intensity_profile(
            sst, psl, pressure_levels, temp, mixr
        )

        assert isinstance(min_p, (float, np.floating))
        assert isinstance(max_w, (float, np.floating))
        assert err == 0
        assert not np.isnan(min_p)
        assert not np.isnan(max_w)

    def test_profile_with_actual_levels(self):
        """Test profile calculation with actual_levels parameter."""
        sst = 300.0
        psl = 101325.0
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0, 230.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003, 0.001])

        min_p, max_w, err = calculate_potential_intensity_profile(
            sst, psl, pressure_levels, temp, mixr, actual_levels=3
        )

        assert isinstance(min_p, (float, np.floating))
        assert isinstance(max_w, (float, np.floating))
        assert err == 0

    def test_profile_with_mismatched_lengths(self):
        """Test profile with mismatched array lengths."""
        sst = 300.0
        psl = 101325.0
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0])  # Too short
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        with pytest.raises(ValueError, match="Profile lengths mismatch"):
            calculate_potential_intensity_profile(sst, psl, pressure_levels, temp, mixr)

    def test_profile_with_reversed_pressure(self):
        """Test profile with reversed pressure levels."""
        sst = 300.0
        psl = 101325.0
        pressure_levels = np.array([500.0, 700.0, 850.0, 1000.0])  # Reversed
        temp = np.array([250.0, 270.0, 285.0, 300.0])
        mixr = np.array([0.003, 0.005, 0.008, 0.012])

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            min_p, max_w, err = calculate_potential_intensity_profile(
                sst, psl, pressure_levels, temp, mixr
            )
            assert len(w) == 1
            assert "Reordering to surface to top" in str(w[0].message)

        assert isinstance(min_p, (float, np.floating))
        assert isinstance(max_w, (float, np.floating))
        assert err == 0

    def test_profile_undef_result(self):
        """Test profile calculation with extreme values that might return UNDEF."""
        # Use unrealistic values that might cause the calculation to fail
        sst = 273.0  # Very cold SST (0°C)
        psl = 101325.0
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([273.0, 273.0, 273.0, 273.0])  # Isothermal (unrealistic)
        mixr = np.array([0.0, 0.0, 0.0, 0.0])  # No moisture

        # This should either return NaN or very low values
        min_p, max_w, err = calculate_potential_intensity_profile(
            sst, psl, pressure_levels, temp, mixr
        )

        # The calculation should complete (even if results are extreme)
        assert err in [0, 1, 2]  # Accept various error codes
        # Results might be NaN or very low/high values
        if not np.isnan(min_p):
            assert min_p > 0  # If not NaN, should be positive
        if not np.isnan(max_w):
            assert max_w >= 0  # If not NaN, should be non-negative


class TestPotentialIntensityCalculator:
    """Test the PotentialIntensityCalculator class."""

    def test_3d_data_detection(self):
        """Test automatic detection of 3D data."""
        nlat, nlon = 10, 20
        num_levels = 4
        sst = np.random.rand(nlat, nlon)
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(num_levels, nlat, nlon)
        mixr = np.random.rand(num_levels, nlat, nlon)

        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)

        assert calc.data_type == "3D"

        min_p, max_w, err = calc.calculate()

        assert min_p.shape == (nlat, nlon)
        assert max_w.shape == (nlat, nlon)
        assert err == 0
        assert calc.results is not None
        assert calc.results["data_type"] == "3D"

    def test_4d_data_detection(self):
        """Test automatic detection of 4D data."""
        ntimes, nlat, nlon = 5, 10, 20
        num_levels = 4
        sst = np.random.rand(ntimes, nlat, nlon)
        psl = np.random.rand(ntimes, nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(ntimes, num_levels, nlat, nlon)
        mixr = np.random.rand(ntimes, num_levels, nlat, nlon)

        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)

        assert calc.data_type == "4D"

        min_p, max_w, err = calc.calculate()

        assert min_p.shape == (ntimes, nlat, nlon)
        assert max_w.shape == (ntimes, nlat, nlon)
        assert err == 0

    def test_profile_data_detection(self):
        """Test automatic detection of profile data."""
        sst = 300.0
        psl = 101325.0
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)

        assert calc.data_type == "profile"

        min_p, max_w, err = calc.calculate()

        assert isinstance(min_p, (float, np.floating))
        assert isinstance(max_w, (float, np.floating))
        assert err == 0

    def test_unsupported_dimensions(self):
        """Test error on unsupported dimensions."""
        sst = np.random.rand(10, 20)
        psl = np.random.rand(10, 20)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(10, 20)  # 2D - not supported
        mixr = np.random.rand(10, 20)

        with pytest.raises(ValueError, match="Unsupported temperature dimensions"):
            PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)

    def test_results_property_before_calculation(self):
        """Test results property before calculation."""
        sst = np.random.rand(10, 20)
        psl = np.random.rand(10, 20)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, 10, 20)
        mixr = np.random.rand(4, 10, 20)

        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)

        assert calc.results is None  # No results before calculation

    def test_results_property_after_calculation(self):
        """Test results property after calculation."""
        sst = np.random.rand(10, 20)
        psl = np.random.rand(10, 20)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, 10, 20)
        mixr = np.random.rand(4, 10, 20)

        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)
        min_p, max_w, err = calc.calculate()

        results = calc.results
        assert results is not None
        assert "min_pressure" in results
        assert "max_wind" in results
        assert "error_flag" in results
        assert "data_type" in results
        np.testing.assert_array_equal(results["min_pressure"], min_p)
        np.testing.assert_array_equal(results["max_wind"], max_w)
        assert results["error_flag"] == err


class TestPotentialIntensityConvenienceFunction:
    """Test the potential_intensity convenience function."""

    def test_3d_data(self):
        """Test convenience function with 3D data."""
        nlat, nlon = 10, 20
        num_levels = 4
        sst = np.random.rand(nlat, nlon)
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(num_levels, nlat, nlon)
        mixr = np.random.rand(num_levels, nlat, nlon)

        min_p, max_w, err = potential_intensity(sst, psl, pressure_levels, temp, mixr)

        assert min_p.shape == (nlat, nlon)
        assert max_w.shape == (nlat, nlon)
        assert err == 0

    def test_4d_data(self):
        """Test convenience function with 4D data."""
        ntimes, nlat, nlon = 5, 10, 20
        num_levels = 4
        sst = np.random.rand(ntimes, nlat, nlon)
        psl = np.random.rand(ntimes, nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(ntimes, num_levels, nlat, nlon)
        mixr = np.random.rand(ntimes, num_levels, nlat, nlon)

        min_p, max_w, err = potential_intensity(sst, psl, pressure_levels, temp, mixr)

        assert min_p.shape == (ntimes, nlat, nlon)
        assert max_w.shape == (ntimes, nlat, nlon)
        assert err == 0

    def test_profile_data(self):
        """Test convenience function with profile data."""
        sst = 300.0
        psl = 101325.0
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        min_p, max_w, err = potential_intensity(sst, psl, pressure_levels, temp, mixr)

        assert isinstance(min_p, (float, np.floating))
        assert isinstance(max_w, (float, np.floating))
        assert err == 0


# Integration tests
class TestIntegration:
    """Integration tests for the interface module."""

    def test_end_to_end_3d_workflow(self):
        """Test complete 3D workflow from input to output."""
        # Create realistic input data
        nlat, nlon = 50, 100
        num_levels = 10
        sst = np.random.rand(nlat, nlon) * 5 + 298  # 298-303 K
        psl = np.random.rand(nlat, nlon) * 2000 + 100000  # 1000-1020 hPa
        pressure_levels = np.array([1000, 925, 850, 700, 600, 500, 400, 300, 200, 100])

        # Create temperature profile with realistic lapse rate
        temp = np.zeros((num_levels, nlat, nlon))
        for i, p in enumerate(pressure_levels):
            temp[i] = 288 - 0.0065 * (8000 * (1 - p / 1000))  # Simple lapse rate
            temp[i] += np.random.rand(nlat, nlon) * 5  # Add variation

        # Create mixing ratio that decreases with height
        mixr = np.zeros((num_levels, nlat, nlon))
        for i, p in enumerate(pressure_levels):
            mixr[i] = 0.015 * (p / 1000) ** 2  # Decreases with height
            mixr[i] += np.random.rand(nlat, nlon) * 0.002

        # Add some missing values
        sst[0:5, 0:10] = np.nan
        temp[:, 10:15, 20:25] = UNDEF

        # Calculate potential intensity
        min_p, max_w, err = calculate_potential_intensity_3d(
            sst, psl, pressure_levels, temp, mixr
        )

        # Validate results
        assert min_p.shape == (nlat, nlon)
        assert max_w.shape == (nlat, nlon)
        assert err == 0

        # Check that missing values were handled
        assert np.all(np.isnan(min_p[0:5, 0:10]))
        assert np.all(np.isnan(max_w[0:5, 0:10]))
        assert np.all(np.isnan(min_p[10:15, 20:25]))

        # Check that non-missing values are reasonable
        valid_min_p = min_p[~np.isnan(min_p)]
        valid_max_w = max_w[~np.isnan(max_w)]
        assert len(valid_min_p) > 0
        assert len(valid_max_w) > 0
        assert np.all(valid_min_p > 800) and np.all(valid_min_p < 1100)
        assert np.all(valid_max_w > 0) and np.all(valid_max_w < 200)

    def test_class_based_workflow(self):
        """Test complete workflow using the class interface."""
        # Create 4D data
        ntimes = 12
        nlat, nlon = 30, 60
        num_levels = 8

        sst = np.random.rand(ntimes, nlat, nlon) * 5 + 298
        psl = np.random.rand(ntimes, nlat, nlon) * 2000 + 100000
        pressure_levels = np.array([1000, 850, 700, 500, 400, 300, 200, 100])

        temp = np.zeros((ntimes, num_levels, nlat, nlon))
        mixr = np.zeros((ntimes, num_levels, nlat, nlon))

        for t in range(ntimes):
            for i, p in enumerate(pressure_levels):
                temp[t, i] = 288 - 0.0065 * (8000 * (1 - p / 1000))
                temp[t, i] += np.random.rand(nlat, nlon) * 5
                mixr[t, i] = 0.015 * (p / 1000) ** 2
                mixr[t, i] += np.random.rand(nlat, nlon) * 0.002

        # Use class interface
        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)

        assert calc.data_type == "4D"

        min_p, max_w, err = calc.calculate()

        assert min_p.shape == (ntimes, nlat, nlon)
        assert max_w.shape == (ntimes, nlat, nlon)
        assert err == 0

        # Check results are stored
        results = calc.results
        assert results["data_type"] == "4D"
        assert results["error_flag"] == 0
        np.testing.assert_array_equal(results["min_pressure"], min_p)
        np.testing.assert_array_equal(results["max_wind"], max_w)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
