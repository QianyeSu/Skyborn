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

# Set proper encoding for Windows and Read the Docs
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Fix locale issues for Read the Docs
try:
    import locale
    locale.setlocale(locale.LC_ALL, 'C.UTF-8')
except:
    try:
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    except:
        # If both fail, set environment variables
        os.environ['LC_ALL'] = 'C.UTF-8'
        os.environ['LANG'] = 'C.UTF-8'

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
        # Set proper encoding and environment for Read the Docs
        env = os.environ.copy()
        env.update({
            'LC_ALL': 'C.UTF-8',
            'LANG': 'C.UTF-8',
            'PYTHONIOENCODING': 'utf-8'
        })

        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            env=env,
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

        # Check if the entrance page has already been set up (avoid duplicate processing)
        if documentation_page.exists() and original_index.exists():
            # Check if index.html is already the entrance page
            try:
                with open(original_index, 'r', encoding='utf-8') as f:
                    index_content = f.read()
                if 'id="shader"' in index_content:  # entrance.html 的特征
                    print("✅ Entrance page already set up properly, skipping...")
                    return True
            except Exception:
                pass

        # 1. If the original index.html exists and is not entrance page, rename it to documentation.html
        if original_index.exists():
            try:
                with open(original_index, 'r', encoding='utf-8') as f:
                    content = f.read()
                # If index.html contains Sphinx documentation features (and not entrance.html), it needs to be moved
                if 'id="shader"' not in content:
                    print("Renaming original index.html to documentation.html...")
                    if documentation_page.exists():
                        documentation_page.unlink()  # Remove existing documentation.html
                    shutil.move(str(original_index), str(documentation_page))
                    print(
                        f"✅ Original documentation moved to {documentation_page}")
            except Exception as e:
                print(f"Error checking index.html content: {e}")

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

    try:
        # Use Sphinx Python API instead of command line to avoid locale issues
        from sphinx.application import Sphinx
        from sphinx.util.docutils import docutils_namespace

        # Set up Sphinx application
        with docutils_namespace():
            app = Sphinx(
                srcdir=str(SOURCE_DIR),
                confdir=str(SOURCE_DIR),
                outdir=str(BUILD_DIR),
                doctreedir=str(BUILD_DIR.parent / '.doctrees'),
                buildername='html',
                confoverrides={},
                status=sys.stdout,
                warning=sys.stderr,
                freshenv=True,
                warningiserror=False,
                tags=[],
                verbosity=1,
                parallel=1
            )

            # Build all documents
            app.build(None)

        print("✅ Documentation built successfully!")

        # Setup entrance page as main index
        setup_entrance_page()

        print(f"   Open: {BUILD_DIR}/index.html (Entry Page)")
        print(f"   Open: {BUILD_DIR}/documentation.html (Documentation)")
        return True

    except Exception as e:
        print("❌ Failed to build documentation")
        print(f"Error: {e}")
        return False


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
