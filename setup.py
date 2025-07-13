"""
Setup script for Skyborn - Mixed build system with meson for Fortran modules
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from setuptools.command.install import install


class MesonBuildExt(build_ext):
    """Custom build extension to handle meson builds for Fortran modules"""

    def run(self):
        """Run the build process"""
        # Build meson modules first
        self.build_meson_modules()
        # Then run the standard build_ext
        super().run()

    def build_meson_modules(self):
        """Build modules that use meson (like spharm)"""
        meson_modules = [
            {
                "name": "spharm",
                "path": Path("src") / "skyborn" / "spharm",
                "target_dir": Path(self.build_lib) / "skyborn" / "spharm",
            }
        ]

        for module in meson_modules:
            if self.should_build_meson_module(module):
                self.build_single_meson_module(module)

    def should_build_meson_module(self, module):
        """Check if we should build this meson module"""
        meson_build_file = module["path"] / "meson.build"
        return meson_build_file.exists()

    def build_single_meson_module(self, module):
        """Build a single meson module"""
        print(f"Building {module['name']} with meson...")

        module_path = module["path"]
        build_dir = module_path / "build"

        try:
            # Check if meson is available
            subprocess.run(
                ["meson", "--version"], check=True, capture_output=True, text=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"Warning: meson not found. Skipping {module['name']} build.")
            print("To install meson: pip install meson ninja")
            return

        try:
            # Clean and setup meson build
            if build_dir.exists():
                shutil.rmtree(build_dir)

            # Use f2py directly since it works better
            print(f"Using f2py to compile {module['name']}...")

            # Get all .f files in src directory
            src_dir = module_path / "src"
            f_files = list(src_dir.glob("*.f"))
            pyf_file = src_dir / "_spherepack.pyf"

            # Build file list
            build_files = [str(pyf_file)] + [str(f) for f in f_files]

            # Run f2py with explicit output directory
            subprocess.run(
                ["python", "-m", "numpy.f2py", "-c"]
                + build_files
                + ["--build-dir", str(module_path)],
                cwd=str(module_path.parent.parent.parent),  # Run from project root
                check=True,
            )

            print(f"✅ {module['name']} compilation successful!")

        except subprocess.CalledProcessError as e:
            print(f"Warning: Failed to build {module['name']}: {e}")
            print(f"Continuing without {module['name']} extensions...")


class CustomDevelop(develop):
    """Custom develop command that builds meson modules"""

    def run(self):
        # Build meson modules in develop mode
        self.run_command("build_ext")
        super().run()


class CustomInstall(install):
    """Custom install command that ensures meson modules are built"""

    def run(self):
        # Ensure meson modules are built before install
        self.run_command("build_ext")
        super().run()


# Configuration for mixed build
setup_config = {
    "cmdclass": {
        "build_ext": MesonBuildExt,
        "develop": CustomDevelop,
        "install": CustomInstall,
    },
    # Add empty Extension for spharm to trigger build_ext
    "ext_modules": [Extension("skyborn.spharm._spherepack", sources=[])],
}

if __name__ == "__main__":
    setup(**setup_config)
