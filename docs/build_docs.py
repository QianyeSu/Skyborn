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
import re
from pathlib import Path

# Set proper encoding for Windows
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Configuration
DOCS_DIR = Path(__file__).parent
SOURCE_DIR = DOCS_DIR / "source"
BUILD_DIR = DOCS_DIR / "build" / "html"


def get_version_from_init():
    """Get version from skyborn/__init__.py"""
    init_file = DOCS_DIR.parent / "src" / "skyborn" / "__init__.py"

    if init_file.exists():
        with open(init_file, 'r', encoding='utf-8') as f:
            content = f.read()
            # Find the version line
            version_match = re.search(
                r'__version__\s*=\s*["\']([^"\']+)["\']', content)
            if version_match:
                return version_match.group(1)

    return "0.3.7"  # fallback version


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


def setup_entrance_page():
    """Setup entrance page as index.html and rename original index.html to documentation.html."""
    entrance_source = SOURCE_DIR / "entrance.html"
    original_index = BUILD_DIR / "index.html"
    documentation_page = BUILD_DIR / "documentation.html"

    if entrance_source.exists():
        print("Setting up entrance page as main index...")

        # 1. If the original index.html exists, rename it to documentation.html
        if original_index.exists():
            print("Renaming original index.html to documentation.html...")
            shutil.move(str(original_index), str(documentation_page))
            print(f"✅ Original documentation moved to {documentation_page}")

        # 2. Read entrance.html and update version information
        with open(entrance_source, 'r', encoding='utf-8') as f:
            content = f.read()

        # Get current version
        current_version = get_version_from_init()

        # Update version information
        content = re.sub(
            r'Version\s+[\d\.]+',
            f'Version {current_version}',
            content
        )

        # 3. Write the new index.html
        with open(original_index, 'w', encoding='utf-8') as f:
            f.write(content)

        print(
            f"✅ Entrance page set as main index with version {current_version}: {original_index}")

        return True
    else:
        print("⚠️  Entrance page not found, skipping...")
        return False


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

        # Setup entrance page as main index
        setup_entrance_page()

        print(f"   Open: {BUILD_DIR}/index.html (Entry Page)")
        print(f"   Open: {BUILD_DIR}/documentation.html (Documentation)")
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
        print(f"🚪 Entry Page: file://{BUILD_DIR.absolute()}/index.html")
        print(
            f"📖 Documentation: file://{BUILD_DIR.absolute()}/documentation.html")
    else:
        print("\n" + "=" * 40)
        print("⚠️  Build failed. Check the output above for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
