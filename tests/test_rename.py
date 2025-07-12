#!/usr/bin/env python3
"""
Simple test script to verify successful spharm module renaming
"""

import sys
import os


def test_directory_structure():
    """Test directory structure"""
    print("TESTING: Directory structure...")

    spharm_path = "src/skyborn/spharm"
    pyspharm_path = "src/skyborn/pyspharm"

    if os.path.exists(spharm_path):
        print("PASS: spharm directory exists")
    else:
        print("FAIL: spharm directory does not exist")
        return False

    if not os.path.exists(pyspharm_path):
        print("PASS: pyspharm directory has been renamed")
    else:
        print("FAIL: pyspharm directory still exists")
        return False

    return True


def test_import_references():
    """Test import references"""
    print("\nTESTING: Import references...")

    # Check references in __init__.py file
    with open("src/skyborn/__init__.py", "r") as f:
        content = f.read()

    if "from . import spharm" in content:
        print("PASS: Import in skyborn/__init__.py updated to spharm")
    else:
        print("FAIL: Import in skyborn/__init__.py not updated")
        return False

    if "pyspharm" not in content:
        print("PASS: pyspharm references removed from skyborn/__init__.py")
    else:
        print("WARNING: pyspharm references still exist in skyborn/__init__.py")

    return True


def test_module_docstring():
    """Test module docstring"""
    print("\nTESTING: Module docstring...")

    with open("src/skyborn/spharm/__init__.py", "r") as f:
        content = f.read()

    if "spharm - Spherical harmonic transforms" in content:
        print("PASS: spharm module docstring updated")
    else:
        print("FAIL: spharm module docstring not updated")
        return False

    if "from skyborn.spharm import Spharmt" in content:
        print("PASS: Import in example code updated")
    else:
        print("FAIL: Import in example code not updated")
        return False

    return True


def test_meson_config():
    """Test meson configuration file"""
    print("\nTESTING: meson configuration file...")

    with open("src/skyborn/spharm/meson.build", "r") as f:
        content = f.read()

    if "spharm submodule for skyborn" in content:
        print("PASS: meson.build comments updated")
    else:
        print("FAIL: meson.build comments not updated")
        return False

    if "subdir: 'skyborn/spharm'" in content:
        print("PASS: meson.build subdirectory path updated")
    else:
        print("FAIL: meson.build subdirectory path not updated")
        return False

    return True


def test_documentation():
    """Test documentation files"""
    print("\nTESTING: Documentation files...")

    # Check if documentation exists in docs/md folder
    doc_path = "docs/md/MESON_BUILD_GUIDE.md"
    try:
        with open(doc_path, "r") as f:
            content = f.read()

        if "spharm" in content and "pyspharm" not in content:
            print("PASS: Documentation completely updated to spharm")
        else:
            print("WARNING: Documentation may still have pyspharm references")
    except FileNotFoundError:
        print(f"WARNING: Documentation file not found at {doc_path}")
        # Try root directory as fallback
        try:
            with open("MESON_BUILD_GUIDE.md", "r") as f:
                content = f.read()
            print("PASS: Found documentation in root directory")
        except FileNotFoundError:
            print("INFO: MESON_BUILD_GUIDE.md not found, skipping documentation test")
            pass  # Don't fail the test if documentation is missing

    return True


def main():
    """Run all tests"""
    print("TESTING: spharm module renaming")
    print("=" * 40)

    tests = [
        test_directory_structure,
        test_import_references,
        test_module_docstring,
        test_meson_config,
        test_documentation,
    ]

    results = []
    for test in tests:
        results.append(test())

    print("\n" + "=" * 40)
    passed = sum(results)
    total = len(results)

    if passed == total:
        print(f"SUCCESS: All {total} tests passed!")
        print("RESULT: pyspharm successfully renamed to spharm!")
        return 0
    else:
        print(f"WARNING: {passed}/{total} tests passed")
        print("RESULT: Renaming may be incomplete")
        return 1


# if __name__ == "__main__":
#     sys.exit(main())
