#!/usr/bin/env python3
"""
Installation helper script for skyborn on Windows
This script helps users install the required build dependencies
"""

import subprocess
import sys
import os
import platform


def run_command(cmd, shell=False):
    """Run a command and return success status"""
    try:
        result = subprocess.run(
            cmd, shell=shell, check=True, capture_output=True, text=True
        )
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        return False, e.stderr


def check_compiler(compiler):
    """Check if a compiler is available"""
    success, output = run_command([compiler, "--version"])
    return success


def install_windows_deps():
    """Install Windows build dependencies"""
    print("🔧 Installing Windows build dependencies...")

    # Check if chocolatey is available
    choco_available, _ = run_command(["choco", "--version"])

    if choco_available:
        print("📦 Installing MinGW-w64 via Chocolatey...")
        success, output = run_command(
            ["choco", "install", "mingw", "--force"], shell=True
        )
        if success:
            print("✅ MinGW-w64 installed successfully!")
            print(
                "⚠️  Please restart your terminal and add C:\\tools\\mingw64\\bin to your PATH"
            )
        else:
            print("❌ Failed to install MinGW-w64 via Chocolatey")
            return False
    else:
        print(
            "❌ Chocolatey not found. Please install it first or manually install MinGW-w64"
        )
        print("   Chocolatey: https://chocolatey.org/install")
        print("   MinGW-w64: https://www.mingw-w64.org/downloads/")
        return False

    return True


def install_conda_deps():
    """Install dependencies via conda"""
    print("🐍 Installing conda dependencies...")

    conda_available, _ = run_command(["conda", "--version"])
    if not conda_available:
        print("❌ Conda not found. Please install Anaconda or Miniconda first.")
        return False

    # Install m2w64-toolchain for Windows
    if platform.system() == "Windows":
        print("📦 Installing m2w64-toolchain...")
        success, output = run_command(
            ["conda", "install", "-c", "conda-forge", "m2w64-toolchain", "-y"]
        )
        if not success:
            print("❌ Failed to install m2w64-toolchain")
            return False

    # Install other dependencies
    print("📦 Installing build dependencies...")
    deps = ["meson>=0.64", "ninja", "numpy>=1.23.0"]
    success, output = run_command(
        ["conda", "install", "-c", "conda-forge"] + deps + ["-y"]
    )

    if success:
        print("✅ Dependencies installed successfully!")
        return True
    else:
        print("❌ Failed to install dependencies")
        return False


def main():
    """Main installation helper"""
    print("🚀 Skyborn Installation Helper")
    print("=" * 40)

    system = platform.system()
    print(f"🖥️  Detected OS: {system}")

    # Check if compilers are already available
    has_gcc = check_compiler("gcc")
    has_gfortran = check_compiler("gfortran")

    print(f"🔍 GCC available: {'✅' if has_gcc else '❌'}")
    print(f"🔍 GFortran available: {'✅' if has_gfortran else '❌'}")

    if has_gcc and has_gfortran:
        print("✅ Compilers already available! You can install skyborn directly:")
        print("   pip install skyborn")
        return True

    if system == "Windows":
        print("\n🔧 Setting up Windows build environment...")

        # Try conda first (easier)
        conda_success = install_conda_deps()
        if conda_success:
            print("\n✅ Setup complete! You can now install skyborn:")
            print("   pip install skyborn")
            return True

        # Fall back to chocolatey
        choco_success = install_windows_deps()
        if choco_success:
            print("\n✅ Setup complete! Please restart your terminal and run:")
            print("   pip install skyborn")
            return True

        print("\n❌ Automatic setup failed. Manual installation required:")
        print("   1. Install MinGW-w64: https://www.mingw-w64.org/downloads/")
        print("   2. Add MinGW bin directory to PATH")
        print("   3. Run: pip install skyborn")

    elif system == "Darwin":  # macOS
        print("\n🍎 For macOS, install Xcode Command Line Tools or use Homebrew:")
        print("   brew install gcc")
        print("   pip install skyborn")

    elif system == "Linux":
        print("\n🐧 For Linux, install build tools via your package manager:")
        print("   Ubuntu/Debian: sudo apt install build-essential gfortran")
        print("   RHEL/CentOS: sudo yum install gcc-gfortran")
        print("   pip install skyborn")

    return False


if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n❌ Installation cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        sys.exit(1)
