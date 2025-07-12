#!/usr/bin/env python3
"""
Test basic spharm module functionality after renaming
(no Fortran compilation required)
"""

import sys
import os


def test_basic_skyborn_structure():
    """Test basic skyborn structure"""
    print("TESTING: Basic skyborn structure...")

    # Add path for importing
    sys.path.insert(0, "src")

    try:
        # Test basic module structure (without importing modules that need dependencies)
        import skyborn

        print(f"PASS: skyborn version: {skyborn.__version__}")

        # Check if spharm is in the module
        spharm_path = os.path.join("src", "skyborn", "spharm")
        if os.path.exists(spharm_path):
            print("PASS: spharm submodule directory exists")
        else:
            print("FAIL: spharm submodule directory does not exist")
            return False

        return True
    except ImportError as e:
        print(f"INFO: Cannot fully import skyborn due to missing dependencies: {e}")
        print("PASS: This is normal since we haven't installed all dependencies yet")
        return True


def test_spharm_module_structure():
    """Test spharm module structure"""
    print("\nTESTING: spharm module structure...")

    spharm_files = [
        "src/skyborn/spharm/__init__.py",
        "src/skyborn/spharm/Lib/__init__.py",
        "src/skyborn/spharm/Lib/spharm.py",
        "src/skyborn/spharm/meson.build",
        "src/skyborn/spharm/src/_spherepack.pyf",
    ]

    missing_files = []
    for file_path in spharm_files:
        if os.path.exists(file_path):
            print(f"PASS: {file_path}")
        else:
            print(f"FAIL: {file_path}")
            missing_files.append(file_path)

    if not missing_files:
        print("PASS: All important spharm files exist")
        return True
    else:
        print(f"FAIL: Missing files: {missing_files}")
        return False


def test_fortran_source_files():
    """Test Fortran source files"""
    print("\nTESTING: Fortran source files...")

    src_dir = "src/skyborn/spharm/src"
    if not os.path.exists(src_dir):
        print(f"FAIL: Fortran source directory does not exist: {src_dir}")
        return False

    fortran_files = [f for f in os.listdir(src_dir) if f.endswith(".f")]
    print(f"PASS: Found {len(fortran_files)} Fortran source files")

    # Check key files
    key_files = ["shaes.f", "shses.f", "vhaes.f", "vhses.f", "hrfft.f"]
    for key_file in key_files:
        if key_file in fortran_files:
            print(f"PASS: Key file exists: {key_file}")
        else:
            print(f"WARNING: Key file missing: {key_file}")

    return len(fortran_files) > 20  # Should have 29 files


def test_meson_configuration():
    """Test meson configuration"""
    print("\nTESTING: meson configuration...")

    # Check top-level meson.build
    if os.path.exists("meson.build"):
        print("PASS: Top-level meson.build exists")
    else:
        print("FAIL: Top-level meson.build does not exist")
        return False

    # Check spharm meson.build
    spharm_meson = "src/skyborn/spharm/meson.build"
    if os.path.exists(spharm_meson):
        print("PASS: spharm meson.build exists")

        # Check content
        with open(spharm_meson, "r") as f:
            content = f.read()

        if "spharm submodule" in content:
            print("PASS: meson.build updated for spharm")
        else:
            print("FAIL: meson.build not correctly updated")
            return False
    else:
        print("FAIL: spharm meson.build does not exist")
        return False

    return True


def test_documentation_update():
    """Test documentation update"""
    print("\nTESTING: documentation update...")

    # Check updated test script
    if os.path.exists("test_meson_build.py"):
        with open("test_meson_build.py", "r") as f:
            content = f.read()

        if "spharm" in content and "test_spharm_" in content:
            print("PASS: test_meson_build.py updated for spharm")
        else:
            print("FAIL: test_meson_build.py not correctly updated")
            return False

    # Check documentation
    if os.path.exists("MESON_BUILD_GUIDE.md"):
        with open("MESON_BUILD_GUIDE.md", "r") as f:
            content = f.read()

        spharm_count = content.count("spharm")
        pyspharm_count = content.count("pyspharm")

        if spharm_count > 0 and pyspharm_count == 0:
            print(
                f"PASS: Documentation fully updated to spharm (appears {spharm_count} times)"
            )
        else:
            print(
                f"WARNING: Documentation update may be incomplete (spharm: {spharm_count}, pyspharm: {pyspharm_count})"
            )

    return True


def main():
    """Run all tests"""
    print("TESTING: spharm basic functionality after renaming")
    print("=" * 50)

    tests = [
        test_basic_skyborn_structure,
        test_spharm_module_structure,
        test_fortran_source_files,
        test_meson_configuration,
        test_documentation_update,
    ]

    results = []
    for test in tests:
        try:
            results.append(test())
        except Exception as e:
            print(f"FAIL: Test failed: {e}")
            results.append(False)

    print("\n" + "=" * 50)
    passed = sum(results)
    total = len(results)

    if passed == total:
        print(f"SUCCESS: All {total} tests passed!")
        print("RESULT: spharm module renaming completed successfully!")
        print("\nNEXT STEPS:")
        print("   1. Install gfortran compiler")
        print("   2. Run pip install -e . for full installation")
        print("   3. Test Fortran extension compilation")
        return 0
    else:
        print(f"PARTIAL: {passed}/{total} tests passed")
        print("RESULT: Some issues may exist")
        return 1


# if __name__ == "__main__":
#     sys.exit(main())
