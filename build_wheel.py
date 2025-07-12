#!/usr/bin/env python3
"""
Local wheel building script for skyborn

This script builds wheels locally for testing before pushing to CI.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path


def run_command(cmd, cwd=None):
    """Run a command and handle errors"""
    print(f"Running: {cmd}")
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=cwd, check=True, capture_output=True, text=True
        )
        if result.stdout:
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        if e.stdout:
            print(f"STDOUT: {e.stdout}")
        if e.stderr:
            print(f"STDERR: {e.stderr}")
        return False


def check_dependencies():
    """Check if required tools are available"""
    print("Checking dependencies...")

    # Check Python
    python_version = sys.version_info
    print(
        f"Python version: {python_version.major}.{python_version.minor}.{python_version.micro}"
    )

    if python_version < (3, 9):
        print("ERROR: Python 3.9+ required")
        return False

    # Check gfortran
    try:
        result = subprocess.run(
            ["gfortran", "--version"], capture_output=True, text=True
        )
        print(f"gfortran available: {result.stdout.split()[3]}")
    except FileNotFoundError:
        print("ERROR: gfortran not found. Please install gfortran.")
        return False

    # Check required Python packages
    required_packages = ["build", "wheel", "meson", "ninja", "numpy"]
    missing_packages = []

    for package in required_packages:
        try:
            __import__(package)
            print(f"✓ {package} available")
        except ImportError:
            missing_packages.append(package)
            print(f"✗ {package} missing")

    if missing_packages:
        print(f"Installing missing packages: {', '.join(missing_packages)}")
        cmd = f"python -m pip install {' '.join(missing_packages)}"
        if not run_command(cmd):
            return False

    return True


def clean_build_artifacts():
    """Clean previous build artifacts"""
    print("Cleaning build artifacts...")

    dirs_to_clean = [
        "build",
        "builddir",
        "builddir_test",
        "dist",
        "wheelhouse",
        "*.egg-info",
    ]

    for pattern in dirs_to_clean:
        for path in Path(".").glob(pattern):
            if path.exists():
                print(f"Removing {path}")
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()


def build_wheel():
    """Build wheel using python -m build"""
    print("Building wheel...")

    # Use python -m build which is the modern standard
    cmd = "python -m build --wheel --no-isolation"

    if not run_command(cmd):
        print("ERROR: Wheel build failed")
        return False

    # Check if wheel was created
    dist_dir = Path("dist")
    wheels = list(dist_dir.glob("*.whl"))

    if not wheels:
        print("ERROR: No wheel file found in dist/")
        return False

    for wheel in wheels:
        print(f"✓ Built wheel: {wheel}")

    return True


def test_wheel():
    """Test the built wheel"""
    print("Testing wheel...")

    # Find the wheel file
    dist_dir = Path("dist")
    wheels = list(dist_dir.glob("*.whl"))

    if not wheels:
        print("ERROR: No wheel to test")
        return False

    wheel_path = wheels[0]  # Use the first wheel

    # Create a temporary virtual environment for testing
    test_env = Path("test_env")
    if test_env.exists():
        shutil.rmtree(test_env)

    # Create virtual environment
    if not run_command(f"python -m venv {test_env}"):
        return False

    # Determine pip executable
    if os.name == "nt":  # Windows
        pip_exe = test_env / "Scripts" / "pip.exe"
        python_exe = test_env / "Scripts" / "python.exe"
    else:  # Unix-like
        pip_exe = test_env / "bin" / "pip"
        python_exe = test_env / "bin" / "python"

    # Install dependencies and wheel
    if not run_command(f"{pip_exe} install numpy"):
        return False

    if not run_command(f"{pip_exe} install {wheel_path}"):
        return False

    # Test import
    test_cmd = f'{python_exe} -c "from skyborn.spharm import Spharmt; print(\\"SUCCESS: Wheel test passed\\"); sht = Spharmt(8, 6); print(\\"SUCCESS: Spharmt works!\\")"'

    success = run_command(test_cmd)

    # Cleanup
    shutil.rmtree(test_env)

    return success


def main():
    """Main function"""
    print("🏗️  Local Wheel Builder for Skyborn")
    print("=" * 50)

    # Check current directory
    if not Path("pyproject.toml").exists():
        print("ERROR: Must run from project root (where pyproject.toml is)")
        sys.exit(1)

    # Check dependencies
    if not check_dependencies():
        print("❌ Dependency check failed")
        sys.exit(1)

    # Clean previous builds
    clean_build_artifacts()

    # Build wheel
    if not build_wheel():
        print("❌ Wheel build failed")
        sys.exit(1)

    # Test wheel
    if not test_wheel():
        print("❌ Wheel test failed")
        sys.exit(1)

    print("✅ Wheel build and test completed successfully!")
    print("Wheel files:")
    for wheel in Path("dist").glob("*.whl"):
        print(f"  📦 {wheel}")

    print("\nNext steps:")
    print("1. Test the wheel: pip install dist/*.whl")
    print("2. Upload to PyPI: python -m twine upload dist/*.whl")
    print("3. Or commit and tag for CI build: git tag v0.3.8 && git push --tags")


if __name__ == "__main__":
    main()
