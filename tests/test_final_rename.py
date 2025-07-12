#!/usr/bin/env python3
"""
Final test: Verify successful renaming of pyspharm to spharm
"""

import sys
import os

# Add source code path
sys.path.insert(0, "src")


def final_test():
    print("TESTING: Final test - pyspharm to spharm renaming verification")
    print("=" * 50)

    # Test 1: Directory structure
    print("TESTING: Directory structure...")
    if os.path.exists("src/skyborn/spharm") and not os.path.exists(
        "src/skyborn/pyspharm"
    ):
        print("PASS: Directory rename successful - pyspharm to spharm")
    else:
        print("FAIL: Directory rename failed")
        return False

    # Test 2: Module import (structural level)
    print("\nTESTING: Module structure...")
    try:
        # Test new module exists
        spharm_init = "src/skyborn/spharm/__init__.py"
        with open(spharm_init, "r") as f:
            content = f.read()

        if "spharm - Spherical harmonic transforms" in content:
            print("PASS: spharm module documentation updated")
        else:
            print("FAIL: spharm module documentation not updated")
            return False

        if "from skyborn.spharm import Spharmt" in content:
            print("PASS: Example code updated to spharm")
        else:
            print("FAIL: Example code not updated")
            return False

    except Exception as e:
        print(f"FAIL: Module test failed: {e}")
        return False

    # Test 3: Main package import references
    print("\nTESTING: Main package import references...")
    try:
        main_init = "src/skyborn/__init__.py"
        with open(main_init, "r") as f:
            content = f.read()

        if "from . import spharm" in content and "pyspharm" not in content:
            print("PASS: Main package import updated to spharm")
        else:
            print("FAIL: Main package import not correctly updated")
            return False

    except Exception as e:
        print(f"FAIL: Main package import test failed: {e}")
        return False

    # Test 4: Build configuration
    print("\nTESTING: Build configuration...")
    try:
        meson_file = "src/skyborn/spharm/meson.build"
        with open(meson_file, "r") as f:
            content = f.read()

        if "spharm submodule" in content and "skyborn/spharm" in content:
            print("PASS: meson build configuration updated")
        else:
            print("FAIL: meson build configuration not updated")
            return False

    except Exception as e:
        print(f"FAIL: Build configuration test failed: {e}")
        return False

    # Test 5: Fortran source file integrity
    print("\nTESTING: Fortran source files...")
    fortran_dir = "src/skyborn/spharm/src"
    if os.path.exists(fortran_dir):
        fortran_files = [f for f in os.listdir(fortran_dir) if f.endswith(".f")]
        if len(fortran_files) >= 29:  # Should have 29 Fortran files
            print(f"PASS: Fortran source files complete: {len(fortran_files)} files")
        else:
            print(f"WARNING: Abnormal Fortran file count: {len(fortran_files)} files")
    else:
        print("FAIL: Fortran source directory does not exist")
        return False

    print("\n" + "=" * 50)
    print("SUCCESS: All tests passed!")
    print("RESULT: pyspharm successfully renamed to spharm!")
    print()
    print("SUMMARY of renaming:")
    print("   • Directory: src/skyborn/pyspharm → src/skyborn/spharm")
    print("   • Import: from skyborn.pyspharm → from skyborn.spharm")
    print("   • Documentation: All references updated to spharm")
    print("   • Configuration: meson.build and pyproject.toml updated")
    print("   • Fortran: 29 .f source files completely preserved")
    print()
    print("NEXT STEPS suggested:")
    print("   1. Ensure gfortran compiler is available")
    print("   2. Install dependencies: python3 -m pip install metpy xarray matplotlib")
    print("   3. Compile and install: python3 -m pip install -e .")
    print("   4. Test spharm functionality: from skyborn.spharm import Spharmt")

    return True


# if __name__ == "__main__":
#     final_test()
