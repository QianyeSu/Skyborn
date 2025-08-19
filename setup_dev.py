#!/usr/bin/env python3
"""
Setup script: Install development environment and pre-commit hooks

Run this script to set up the complete development environment, including pre-commit hooks.
"""

import json
import subprocess
import sys
from pathlib import Path


def run_command(cmd, check=True, shell=False):
    """Run command and print output"""
    print(f"Running: {cmd}")
    try:
        if shell:
            result = subprocess.run(cmd, shell=True, check=check, text=True)
        else:
            result = subprocess.run(cmd.split(), check=check, text=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")
        return False


def check_python_version():
    """Check Python version compatibility"""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 9):
        print("Python 3.9+ is required")
        return False
    print(f"Python {version.major}.{version.minor}.{version.micro}")
    return True


def install_dev_dependencies():
    """Install development dependencies"""
    print("\nInstalling development dependencies...")

    # Install package itself in editable mode
    if not run_command("pip install -e ."):
        return False

    # Install development tools
    dev_packages = [
        "pre-commit",
        "black",
        "isort",
        "flake8",
        "mypy",
        "pytest",
        "pytest-cov",
        "bandit",
        "safety",
    ]

    for package in dev_packages:
        if not run_command(f"pip install {package}"):
            print(f"Failed to install {package}, continuing...")

    return True


def setup_pre_commit():
    """Set up pre-commit hooks"""
    print("\nSetting up pre-commit hooks...")

    # Install pre-commit hooks
    if not run_command("pre-commit install"):
        return False

    # Run once to ensure it works properly
    print("Running pre-commit on all files (first time setup)...")
    run_command("pre-commit run --all-files", check=False)

    return True


def create_vscode_settings():
    """Create VS Code settings for development"""
    print("\nCreating VS Code settings...")

    vscode_dir = Path(".vscode")
    vscode_dir.mkdir(exist_ok=True)

    settings = {
        "python.defaultInterpreterPath": "./venv/bin/python",
        "python.formatting.provider": "black",
        "python.linting.enabled": True,
        "python.linting.flake8Enabled": True,
        "python.linting.mypyEnabled": True,
        "python.sortImports.args": ["--profile", "black"],
        "editor.formatOnSave": True,
        "editor.codeActionsOnSave": {"source.organizeImports": True},
        "[python]": {"editor.rulers": [88], "editor.tabSize": 4},
    }

    with open(vscode_dir / "settings.json", "w") as f:
        json.dump(settings, f, indent=2)

    print("VS Code settings created")


def main():
    """Main function to set up development environment"""
    print("Setting up Skyborn development environment")
    print("=" * 50)

    # Check Python version
    if not check_python_version():
        sys.exit(1)

    # Install dependencies
    if not install_dev_dependencies():
        print("Failed to install dependencies")
        sys.exit(1)

    # Set up pre-commit
    if not setup_pre_commit():
        print("Failed to setup pre-commit")
        sys.exit(1)

    # Create VS Code settings
    create_vscode_settings()

    print("\n" + "=" * 50)
    print("Development environment setup complete!")
    print("\nNext steps:")
    print("1. Start developing code")
    print("2. Pre-commit will run automatically on each commit")
    print("3. Run 'pre-commit run --all-files' to manually check all files")
    print("4. Run 'pytest' to execute tests")
    print("\nUseful commands:")
    print("  pre-commit run --all-files  # Run all checks")
    print("  black src/                  # Format code")
    print("  isort src/                  # Sort imports")
    print("  flake8 src/                 # Code quality check")
    print("  pytest tests/               # Run tests")


if __name__ == "__main__":
    main()
