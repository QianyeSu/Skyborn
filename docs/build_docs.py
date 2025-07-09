#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Documentation builder for Skyborn.

This script builds the English documentation using Sphinx.
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path

# Set proper encoding for Windows
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Configuration
DOCS_DIR = Path(__file__).parent
SOURCE_DIR = DOCS_DIR / "source"
BUILD_DIR = DOCS_DIR / "build" / "html"


def run_command(cmd, cwd=None):
    """Run a shell command and return the result."""
    try:
        # Set proper encoding for Windows
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        if e.stderr:
            print(f"Error output: {e.stderr}")
        return None


def check_dependencies():
    """Check if required packages are installed."""
    required_packages = [
        'sphinx',
        'sphinx-book-theme',
        'myst-nb',
        'sphinx-autodoc-typehints',
        'sphinx-copybutton',
        'toml'
    ]

    missing = []
    for package in required_packages:
        try:
            __import__(package.replace('-', '_'))
        except ImportError:
            missing.append(package)

    if missing:
        print("Missing required packages:")
        for pkg in missing:
            print(f"  - {pkg}")
        print("\nInstall with:")
        print(f"pip install {' '.join(missing)}")
        return False

    return True


def clean_build():
    """Clean the build directory."""
    if BUILD_DIR.parent.exists():
        print("Cleaning build directory...")
        shutil.rmtree(BUILD_DIR.parent)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)


def build_documentation():
    """Build the documentation."""
    print("\nBuilding Skyborn documentation...")
    print(f"Source: {SOURCE_DIR}")
    print(f"Output: {BUILD_DIR}")

    # Create output directory
    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = f"sphinx-build -b html {SOURCE_DIR} {BUILD_DIR}"

    # Set environment variables for proper encoding
    env = os.environ.copy()
    env['PYTHONIOENCODING'] = 'utf-8'
    env['LANG'] = 'en_US.UTF-8'
    env['LC_ALL'] = 'en_US.UTF-8'

    try:
        # Run sphinx-build with proper encoding
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=DOCS_DIR,
            capture_output=True,
            text=True,
            encoding='utf-8',
            env=env,
            check=True
        )
        print("✅ Documentation built successfully!")
        print(f"   Open: {BUILD_DIR}/index.html")
        return True
    except subprocess.CalledProcessError as e:
        print("❌ Failed to build documentation")
        print(f"Command: {cmd}")
        if e.stdout:
            print(f"Output: {e.stdout}")
        if e.stderr:
            print(f"Error: {e.stderr}")
def main():
    """Main build function."""
    parser = argparse.ArgumentParser(description='Build Skyborn documentation')
    parser.add_argument('--clean', action='store_true',
                        help='Clean build directory before building')

    args = parser.parse_args()

    print("🚀 Building Skyborn Documentation")
    print("=" * 40)

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    # Clean build directory if requested
    if args.clean:
        clean_build()

    # Build documentation
    if build_documentation():
        print("\n" + "=" * 40)
        print("🎉 Documentation built successfully!")
        print(f"📖 View documentation: file://{BUILD_DIR.absolute()}/index.html")
    else:
        print("\n" + "=" * 40)
        print("⚠️  Build failed. Check the output above for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
