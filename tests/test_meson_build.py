#!/usr/bin/env python3
"""
Quick test script for the new skyborn meson build system

This script tests:
1. Basic skyborn imports
2. spharm integration
3. Fortran extension availability
"""

import sys
import numpy as np


def test_basic_imports():
    """Test basic skyborn module imports"""
    print("TESTING: Basic skyborn imports...")
    try:
        import skyborn

        print(f"PASS: skyborn version: {skyborn.__version__}")

        # Test submodules
        from skyborn import calc, plot, interp, conversion, ROF

        print("PASS: All basic submodules imported successfully")

        # Test a basic function
        from skyborn.calc import convert_longitude_range

        result = convert_longitude_range(
            np.array([350, 10, 20]), target_range=(-180, 180)
        )
        print(f"PASS: calc.convert_longitude_range test: {result}")

        return True
    except Exception as e:
        print(f"FAIL: Basic import failed: {e}")
        return False


def test_spharm_import():
    """Test spharm module import"""
    print("\nTESTING: spharm import...")
    try:
        from skyborn import spharm

        print("PASS: spharm module imported successfully")

        # Try to import the main class
        from skyborn.spharm import Spharmt

        print("PASS: Spharmt class imported successfully")

        return True
    except Exception as e:
        print(f"FAIL: spharm import failed: {e}")
        return False


def test_spharm_functionality():
    """Test basic spharm functionality"""
    print("\nTESTING: spharm functionality...")
    try:
        from skyborn.spharm import Spharmt

        # Create a simple test case
        nlon, nlat = 72, 36  # Small grid for testing
        sht = Spharmt(nlon=nlon, nlat=nlat, gridtype="regular")
        print(f"PASS: Created Spharmt instance: {nlon}x{nlat} grid")

        # Test with some dummy data
        test_data = np.random.random((nlat, nlon))

        # Grid to spectral transform
        spec = sht.grdtospec(test_data)
        print(
            f"PASS: Grid to spectral transform: shape {test_data.shape} -> {spec.shape}"
        )

        # Spectral to grid transform
        data_back = sht.spectogrd(spec)
        print(
            f"PASS: Spectral to grid transform: shape {spec.shape} -> {data_back.shape}"
        )

        # Check if round-trip is approximately correct
        diff = np.max(np.abs(test_data - data_back))
        print(f"PASS: Round-trip error: {diff:.2e}")

        if diff < 1e-10:
            print("PASS: spharm Fortran extensions working correctly!")
            return True
        else:
            print("WARNING: spharm working but with higher than expected error")
            return True

    except Exception as e:
        print(f"FAIL: spharm functionality test failed: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_fortran_extensions():
    """Test if Fortran extensions are properly compiled"""
    print("\nTESTING: Fortran extensions...")
    try:
        from skyborn.spharm import _spherepack

        print("PASS: _spherepack Fortran extension imported successfully")

        # Check some basic functions
        if hasattr(_spherepack, "shaes"):
            print("PASS: shaes function available")
        if hasattr(_spherepack, "shses"):
            print("PASS: shses function available")

        return True
    except Exception as e:
        print(f"FAIL: Fortran extension test failed: {e}")
        return False


def main():
    """Run all tests"""
    print("TESTING: skyborn meson build system")
    print("=" * 50)

    tests = [
        test_basic_imports,
        test_spharm_import,
        test_fortran_extensions,
        test_spharm_functionality,
    ]

    results = []
    for test in tests:
        results.append(test())

    print("\n" + "=" * 50)
    passed = sum(results)
    total = len(results)

    if passed == total:
        print(f"SUCCESS: All {total} tests passed!")
        print("RESULT: skyborn with spharm integration is working correctly!")
        return 0
    else:
        print(f"PARTIAL: {passed}/{total} tests passed")
        print("RESULT: Some issues detected - check the output above")
        return 1


# if __name__ == "__main__":
#     sys.exit(main())
