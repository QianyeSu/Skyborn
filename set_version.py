#!/usr/bin/env python3
"""
Intelligent version synchronization script for skyborn.

Usage:
    python set_version.py [new_version]

This script:
1. If no version specified: checks if _version.py and pyproject.toml are consistent
2. If inconsistent: uses _version.py as source of truth and updates pyproject.toml
3. If version specified: updates both files to the specified version
"""

import sys
import re
from pathlib import Path


def get_version_from_file(file_path, pattern):
    """Extract version from a file using regex pattern"""
    if not file_path.exists():
        return None

    content = file_path.read_text()
    match = re.search(pattern, content)
    return match.group(1) if match else None


def get_current_versions():
    """Get current versions from both files"""
    version_file = Path("src/skyborn/_version.py")
    pyproject_file = Path("pyproject.toml")

    # Get version from _version.py
    version_py = get_version_from_file(
        version_file, r'__version__ = ["\']([^"\']+)["\']'
    )

    # Get version from pyproject.toml
    version_toml = get_version_from_file(
        pyproject_file, r'version = ["\']([^"\']+)["\']'
    )

    return version_py, version_toml


def update_version(new_version):
    """Update version in both _version.py and pyproject.toml"""

    # Update _version.py
    version_file = Path("src/skyborn/_version.py")
    if version_file.exists():
        content = version_file.read_text()
        content = re.sub(
            r'__version__ = ["\'][^"\']+["\']',
            f'__version__ = "{new_version}"',
            content,
        )
        version_file.write_text(content)
        print(f"✓ Updated _version.py to {new_version}")
    else:
        print("✗ _version.py not found!")
        return False

    # Update pyproject.toml
    pyproject_file = Path("pyproject.toml")
    if pyproject_file.exists():
        content = pyproject_file.read_text()

        # Find and replace the version line in [project] section
        lines = content.split("\n")
        in_project_section = False

        for i, line in enumerate(lines):
            if line.strip() == "[project]":
                in_project_section = True
            elif line.startswith("[") and not line.startswith("[project"):
                in_project_section = False
            elif in_project_section and line.startswith("version ="):
                lines[i] = (
                    f'version = "{new_version}"  # Synchronized with src/skyborn/_version.py'
                )
                break

        content = "\n".join(lines)
        pyproject_file.write_text(content)
        print(f"✓ Updated pyproject.toml to {new_version}")
    else:
        print("✗ pyproject.toml not found!")
        return False

    return True


def main():
    # Get current versions
    version_py, version_toml = get_current_versions()

    if version_py is None:
        print("✗ Could not read version from _version.py")
        sys.exit(1)

    if version_toml is None:
        print("✗ Could not read version from pyproject.toml")
        sys.exit(1)

    print(f"Current versions:")
    print(f"  _version.py: {version_py}")
    print(f"  pyproject.toml: {version_toml}")

    # Check if version is provided as argument
    if len(sys.argv) == 2:
        new_version = sys.argv[1]

        # Simple version validation
        if not re.match(r"^\d+\.\d+\.\d+", new_version):
            print("Error: Version should be in format X.Y.Z")
            sys.exit(1)

        print(f"\nSetting version to {new_version}...")

        if update_version(new_version):
            print("Version update completed successfully!")
            print(f"All files now use version {new_version}")
        else:
            print("Version update failed!")
            sys.exit(1)

    elif len(sys.argv) == 1:
        # No version specified, check consistency
        if version_py == version_toml:
            print(f"\n✓ Versions are consistent ({version_py})")
            print("No action needed.")
        else:
            print(f"\n⚠️  Versions are inconsistent!")
            print(f"Using _version.py ({version_py}) as source of truth...")

            if update_version(version_py):
                print("Version synchronization completed successfully!")
                print(f"All files now use version {version_py}")
            else:
                print("Version synchronization failed!")
                sys.exit(1)

    else:
        print("Usage: python set_version.py [new_version]")
        print("Examples:")
        print("  python set_version.py        # Check and sync versions")
        print("  python set_version.py 0.3.8  # Set specific version")
        sys.exit(1)


if __name__ == "__main__":
    main()
