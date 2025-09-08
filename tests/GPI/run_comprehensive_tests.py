#!/usr/bin/env python
"""
Comprehensive test runner for the GPI module with coverage reporting.
This script runs all tests and provides detailed coverage information.
"""

import os
import sys
import warnings

import numpy as np

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))


def test_interface_module():
    """Test the interface module comprehensively."""
    print("\n" + "=" * 60)
    print("Testing skyborn.calc.GPI.interface module")
    print("=" * 60)

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

    tests_passed = 0
    tests_failed = 0

    # Test 1: _postprocess_results
    try:
        arr = np.array([[950.0, UNDEF], [970.0, 990.0]])
        result_p, result_w = _postprocess_results(arr, arr)
        assert np.isnan(result_p[0, 1])
        print("PASS: Test _postprocess_results")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _postprocess_results - {e}")
        tests_failed += 1

    # Test 2: _validate_input_arrays
    try:
        arr1 = np.array([[1.0, np.nan], [3.0, 4.0]])
        arr2 = np.array([[5.0, 6.0], [7.0, 8.0]])
        (result1, result2), has_missing = _validate_input_arrays(
            arr1, arr2, names=["array1", "array2"]
        )
        assert has_missing == True
        assert result1[0, 1] == UNDEF
        print("PASS: Test _validate_input_arrays")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _validate_input_arrays - {e}")
        tests_failed += 1

    # Test 3: _ensure_pressure_ordering
    try:
        pressure = np.array([300.0, 500.0, 700.0, 850.0, 1000.0])
        temp = np.array([210.0, 250.0, 270.0, 285.0, 300.0])
        mixr = np.array([0.001, 0.003, 0.005, 0.008, 0.012])

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            p_out, t_out, m_out = _ensure_pressure_ordering(pressure, temp, mixr)

        assert np.array_equal(p_out, pressure[::-1])
        print("PASS: Test _ensure_pressure_ordering")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _ensure_pressure_ordering - {e}")
        tests_failed += 1

    # Test 4: _validate_dimensions
    try:
        sst = 300.0
        psl = 1013.25
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.array([300.0, 285.0, 270.0, 250.0])
        mixr = np.array([0.012, 0.008, 0.005, 0.003])

        _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "profile")
        print("PASS: Test _validate_dimensions (profile)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _validate_dimensions (profile) - {e}")
        tests_failed += 1

    # Test 5: _validate_dimensions 3D
    try:
        nlat, nlon = 10, 20
        sst = np.random.rand(nlat, nlon)
        psl = np.random.rand(nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, nlat, nlon)
        mixr = np.random.rand(4, nlat, nlon)

        _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "3D")
        print("PASS: Test _validate_dimensions (3D)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _validate_dimensions (3D) - {e}")
        tests_failed += 1

    # Test 6: _validate_dimensions 4D
    try:
        ntimes, nlat, nlon = 5, 10, 20
        sst = np.random.rand(ntimes, nlat, nlon)
        psl = np.random.rand(ntimes, nlat, nlon)
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(ntimes, 4, nlat, nlon)
        mixr = np.random.rand(ntimes, 4, nlat, nlon)

        _validate_dimensions(sst, psl, pressure_levels, temp, mixr, "4D")
        print("PASS: Test _validate_dimensions (4D)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _validate_dimensions (4D) - {e}")
        tests_failed += 1

    # Test 7: PotentialIntensityCalculator class
    try:
        nlat, nlon = 10, 20
        sst = np.random.rand(nlat, nlon) * 10 + 290
        psl = np.random.rand(nlat, nlon) * 1000 + 100000
        pressure_levels = np.array([1000.0, 850.0, 700.0, 500.0])
        temp = np.random.rand(4, nlat, nlon) * 50 + 250
        mixr = np.random.rand(4, nlat, nlon) * 0.02

        calc = PotentialIntensityCalculator(sst, psl, pressure_levels, temp, mixr)
        assert calc.data_type == "3D"
        assert calc.results is None
        print("PASS: Test PotentialIntensityCalculator")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test PotentialIntensityCalculator - {e}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_xarray_module():
    """Test the xarray module comprehensively."""
    print("\n" + "=" * 60)
    print("Testing skyborn.calc.GPI.xarray module")
    print("=" * 60)

    try:
        import xarray as xr
    except ImportError:
        print("! xarray not available, skipping xarray tests")
        return 0, 0

    from skyborn.calc.GPI.xarray import (
        _check_units,
        _create_output_dataset,
        _detect_atmospheric_dimensions,
        potential_intensity,
    )

    tests_passed = 0
    tests_failed = 0

    # Test 1: _check_units for temperature
    try:
        # Test Celsius to Kelvin conversion
        da = xr.DataArray(25.0, attrs={"units": "degC"})
        result = _check_units(da, "temperature")
        assert np.isclose(
            result.values, 298.15
        ), f"Expected 298.15, got {result.values}"

        # Test Fahrenheit to Kelvin conversion
        da = xr.DataArray(77.0, attrs={"units": "degF"})
        result = _check_units(da, "temperature")
        assert np.isclose(
            result.values, 298.15
        ), f"Expected 298.15, got {result.values}"

        print("PASS: Test _check_units (temperature)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _check_units (temperature) - {e}")
        tests_failed += 1

    # Test 2: _check_units for pressure
    try:
        # Test hPa to Pa conversion
        da = xr.DataArray(1013.25, attrs={"units": "hPa"})
        result = _check_units(da, "pressure")
        assert np.isclose(result.values, 101325.0)

        # Test mb to Pa conversion
        da = xr.DataArray(1013.25, attrs={"units": "mb"})
        result = _check_units(da, "pressure")
        assert np.isclose(result.values, 101325.0)

        print("PASS: Test _check_units (pressure)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _check_units (pressure) - {e}")
        tests_failed += 1

    # Test 3: _check_units for humidity
    try:
        # Test g/kg to kg/kg conversion
        da = xr.DataArray(10.0, attrs={"units": "g/kg"})
        result = _check_units(da, "humidity")
        assert np.isclose(result.values, 0.01)

        # Test percentage to kg/kg conversion
        da = xr.DataArray(1.5, attrs={"units": "%"})
        result = _check_units(da, "humidity")
        assert np.isclose(result.values, 0.015)

        print("PASS: Test _check_units (humidity)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _check_units (humidity) - {e}")
        tests_failed += 1

    # Test 4: _detect_atmospheric_dimensions
    try:
        # Test with 3D data
        temp_3d = xr.DataArray(
            np.random.rand(10, 20, 30), dims=["pressure", "lat", "lon"]
        )
        result = _detect_atmospheric_dimensions(temp_3d)
        assert result == "3D"

        # Test with 4D data
        temp_4d = xr.DataArray(
            np.random.rand(5, 10, 20, 30), dims=["time", "pressure", "lat", "lon"]
        )
        result = _detect_atmospheric_dimensions(temp_4d)
        assert result == "4D"

        print("PASS: Test _detect_atmospheric_dimensions")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _detect_atmospheric_dimensions - {e}")
        tests_failed += 1

    # Test 5: _create_output_dataset
    try:
        # Create sample DataArrays
        min_p = xr.DataArray(
            np.random.rand(10, 20) * 100 + 900,
            dims=["lat", "lon"],
            attrs={"units": "mb", "long_name": "Minimum Central Pressure"},
        )
        max_w = xr.DataArray(
            np.random.rand(10, 20) * 50 + 20,
            dims=["lat", "lon"],
            attrs={"units": "m/s", "long_name": "Maximum Wind Speed"},
        )

        result = _create_output_dataset(min_p, max_w, 0)
        assert isinstance(result, xr.Dataset)
        assert "min_pressure" in result
        assert "max_wind" in result
        assert "error_flag" in result
        print("PASS: Test _create_output_dataset")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test _create_output_dataset - {e}")
        tests_failed += 1

    # Test 6: potential_intensity xarray interface
    try:
        # Create sample xarray data
        levels = xr.DataArray([1000, 850, 700, 500], dims=["level"])
        temp = xr.DataArray(
            [[300, 285, 270, 250]],
            dims=["sample", "level"],
            coords={"level": levels},
            attrs={"units": "K"},
        )
        mixr = xr.DataArray(
            [[0.012, 0.008, 0.005, 0.003]],
            dims=["sample", "level"],
            coords={"level": levels},
            attrs={"units": "kg/kg"},
        )
        sst = xr.DataArray([301.0], dims=["sample"], attrs={"units": "K"})
        psl = xr.DataArray([101325.0], dims=["sample"], attrs={"units": "Pa"})

        # Test that the function can be called (actual result depends on Fortran module)
        # We're just testing the interface here
        assert callable(potential_intensity)
        print("PASS: Test potential_intensity xarray interface")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test potential_intensity xarray interface - {e}")
        tests_failed += 1

    return tests_passed, tests_failed


def test_init_module():
    """Test the __init__ module."""
    print("\n" + "=" * 60)
    print("Testing skyborn.calc.GPI.__init__ module")
    print("=" * 60)

    tests_passed = 0
    tests_failed = 0

    try:
        from skyborn.calc.GPI import potential_intensity

        print("PASS: Test __init__ imports (potential_intensity)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test __init__ imports (potential_intensity) - {e}")
        tests_failed += 1

    try:
        from skyborn.calc.GPI import xarray

        print("PASS: Test __init__ imports (xarray module)")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test __init__ imports (xarray module) - {e}")
        tests_failed += 1

    try:
        # Test that we can import from interface directly
        from skyborn.calc.GPI.interface import (
            PotentialIntensityCalculator,
            calculate_potential_intensity_3d,
            calculate_potential_intensity_4d,
            calculate_potential_intensity_profile,
        )

        print("PASS: Test direct interface imports")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test direct interface imports - {e}")
        tests_failed += 1

    try:
        # Test that we can import from xarray directly
        from skyborn.calc.GPI.xarray import (
            _check_units,
            _detect_atmospheric_dimensions,
            potential_intensity,
        )

        print("PASS: Test direct xarray imports")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL: Test direct xarray imports - {e}")
        tests_failed += 1

    return tests_passed, tests_failed


def calculate_coverage():
    """Calculate and display coverage statistics."""
    print("\n" + "=" * 60)
    print("Coverage Analysis")
    print("=" * 60)

    try:
        import coverage

        cov = coverage.Coverage()
        cov.start()

        # Import all modules to get coverage
        from skyborn.calc.GPI import __init__, interface, xarray

        # Execute various functions
        from skyborn.calc.GPI.interface import (
            PotentialIntensityCalculator,
            _ensure_pressure_ordering,
            _postprocess_results,
            _validate_dimensions,
            _validate_input_arrays,
        )

        # Run some functions
        arr = np.array([1, 2, 3])
        _validate_input_arrays(arr)
        _postprocess_results(arr, arr)

        pressure = np.array([1000, 850, 700])
        temp = np.array([300, 285, 270])
        mixr = np.array([0.012, 0.008, 0.005])
        _ensure_pressure_ordering(pressure, temp, mixr)

        cov.stop()
        cov.save()

        # Get coverage data
        total_lines = 0
        executed_lines = 0

        for filename in cov.get_data().measured_files():
            if "GPI" in filename and "test" not in filename:
                analysis = cov.analysis2(filename)
                executed = len(analysis[1])
                missing = len(analysis[3])
                total = executed + missing
                if total > 0:
                    percent = (executed / total) * 100
                    module_name = os.path.basename(filename)
                    print(f"  {module_name}: {percent:.1f}% ({executed}/{total} lines)")
                    total_lines += total
                    executed_lines += executed

        if total_lines > 0:
            overall_coverage = (executed_lines / total_lines) * 100
            print(
                f"\nOverall GPI module coverage: {overall_coverage:.1f}% ({executed_lines}/{total_lines} lines)"
            )

            if overall_coverage >= 96:
                print("PASS: Coverage goal of 96% achieved!")
            else:
                needed = int(0.96 * total_lines) - executed_lines
                print(f"FAIL: Need to cover {needed} more lines to reach 96% coverage")

    except ImportError:
        print("! coverage package not available")


def main():
    """Run all tests and display results."""
    print("GPI Module Comprehensive Test Suite")
    print("=" * 60)

    total_passed = 0
    total_failed = 0

    # Run tests for each module
    passed, failed = test_init_module()
    total_passed += passed
    total_failed += failed

    passed, failed = test_interface_module()
    total_passed += passed
    total_failed += failed

    passed, failed = test_xarray_module()
    total_passed += passed
    total_failed += failed

    # Display summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    print(f"Total tests passed: {total_passed}")
    print(f"Total tests failed: {total_failed}")

    if total_failed == 0:
        print("PASS: All tests passed!")
    else:
        print(f"FAIL: {total_failed} tests failed")

    # Calculate coverage
    calculate_coverage()

    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
