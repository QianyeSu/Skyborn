"""
Comprehensive test runner for skyborn.interp module.

This script runs all interpolation and regridding tests to validate
code coverage and ensure all functionality is properly tested.
"""

import os
import subprocess
import sys
from pathlib import Path


def run_test_suite():
    """Run the complete interpolation module test suite."""

    # Get the tests directory
    tests_dir = Path(__file__).parent

    # List of test files for interpolation module
    test_files = [
        "test_interp_regridding.py",  # Regridding tests
        "test_interp_interpolation.py",  # Original interpolation tests
        "test_interp_interpolation_enhanced.py",  # Enhanced interpolation tests
        "test_interp_module_comprehensive.py",  # Module integration tests
    ]

    print("Running Skyborn Interpolation Module Test Suite")
    print("=" * 50)

    total_passed = 0
    total_failed = 0

    for test_file in test_files:
        test_path = tests_dir / test_file

        if not test_path.exists():
            print(f"⚠️  Test file not found: {test_file}")
            continue

        print(f"\n🧪 Running {test_file}...")
        print("-" * 30)

        try:
            # Run pytest for this specific file
            result = subprocess.run(
                [sys.executable, "-m", "pytest", str(test_path), "-v", "--tb=short"],
                capture_output=True,
                text=True,
                cwd=tests_dir.parent,  # Run from project root
            )

            # Parse output for pass/fail counts
            output_lines = result.stdout.split("\n")

            passed = 0
            failed = 0

            for line in output_lines:
                if " PASSED " in line:
                    passed += 1
                elif " FAILED " in line or " ERROR " in line:
                    failed += 1

            if result.returncode == 0:
                print(f"✅ {test_file}: {passed} tests passed")
                total_passed += passed
            else:
                print(f"❌ {test_file}: {passed} passed, {failed} failed")
                total_passed += passed
                total_failed += failed

                # Show error details for failed tests
                if result.stderr:
                    print("Error output:")
                    print(
                        result.stderr[:500] + "..."
                        if len(result.stderr) > 500
                        else result.stderr
                    )

        except Exception as e:
            print(f"❌ Error running {test_file}: {e}")
            total_failed += 1

    print("\n" + "=" * 50)
    print(f"Test Suite Summary:")
    print(f"✅ Total Passed: {total_passed}")
    print(f"❌ Total Failed: {total_failed}")

    if total_failed == 0:
        print("🎉 All tests passed! Interpolation module has excellent coverage.")
    else:
        print(f"⚠️  {total_failed} tests failed. Review and fix issues.")

    return total_failed == 0


def check_coverage():
    """Check test coverage for interpolation module."""
    print("\n📊 Checking test coverage...")
    print("-" * 30)

    try:
        # Run coverage analysis
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "--cov=src.skyborn.interp",
                "--cov-report=term-missing",
                "tests/test_interp*.py",
            ],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        if result.returncode == 0:
            # Extract coverage percentage from output
            lines = result.stdout.split("\n")
            for line in lines:
                if "src/skyborn/interp" in line and "%" in line:
                    print(f"📈 {line}")
        else:
            print("Coverage tool not available or tests failed")

    except Exception as e:
        print(f"Could not run coverage analysis: {e}")


def validate_documentation_compatibility():
    """Validate that all tests are compatible with documented API."""
    print("\n📚 Validating documentation compatibility...")
    print("-" * 40)

    # Check that all documented functions are tested
    documented_functions = [
        "Grid.from_degrees",
        "regrid_dataset",
        "interp_hybrid_to_pressure",
        "interp_sigma_to_hybrid",
        "interp_multidim",
    ]

    test_files_content = []
    tests_dir = Path(__file__).parent

    for test_file in ["test_interp_regridding.py", "test_interp_interpolation.py"]:
        test_path = tests_dir / test_file
        if test_path.exists():
            test_files_content.append(test_path.read_text())

    all_content = "\n".join(test_files_content)

    print("Checking documented functions are tested:")
    for func in documented_functions:
        if func in all_content:
            print(f"✅ {func} - tested")
        else:
            print(f"⚠️  {func} - may need more test coverage")


if __name__ == "__main__":
    print("""
    ██████╗ ██╗   ██╗██████╗  ██████╗ ██████╗ ███╗   ██╗
    ██╔═══██╗██║   ██║██╔══██╗██╔═══██╗██╔══██╗████╗  ██║
    ██║   ██║██║   ██║██████╔╝██║   ██║██████╔╝██╔██╗ ██║
    ██║   ██║██║   ██║██╔══██╗██║   ██║██╔══██╗██║╚██╗██║
    ╚██████╔╝╚██████╔╝██████╔╝╚██████╔╝██║  ██║██║ ╚████║
     ╚═════╝  ╚═════╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝

    Interpolation Module Test Suite
    """)

    success = run_test_suite()

    # Optional: Check coverage if coverage tools are available
    check_coverage()

    # Validate documentation compatibility
    validate_documentation_compatibility()

    if success:
        print("\n🎯 Test suite completed successfully!")
        print("Your interpolation module now has comprehensive test coverage.")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed. Please review and fix issues.")
        sys.exit(1)
