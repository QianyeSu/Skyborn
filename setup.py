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
import numpy as np

# Check if Cython is available
try:
    from Cython.Build import cythonize

    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

# Force gfortran compiler usage
os.environ["FC"] = os.environ.get("FC", "gfortran")
os.environ["F77"] = os.environ.get("F77", "gfortran")
os.environ["F90"] = os.environ.get("F90", "gfortran")
os.environ["CC"] = os.environ.get("CC", "gcc")


# Check if gfortran is available
def check_gfortran():
    """Check if gfortran is available"""
    try:
        result = subprocess.run(
            ["gfortran", "--version"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            print(
                f"Found gfortran: {result.stdout.split()[4] if len(result.stdout.split()) > 4 else 'unknown version'}"
            )
            return True
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
        pass

    print("Warning: gfortran not found. Fortran extensions may not build correctly.")
    print("Please install gfortran:")
    print("  Linux: sudo apt-get install gfortran")
    print("  macOS: brew install gcc")
    print("  Windows: conda install m2w64-toolchain")
    return False


# Check gfortran availability at setup time
check_gfortran()


def get_gridfill_extensions():
    """Get Cython extensions for gridfill module with cross-platform optimizations"""
    extensions = []

    if HAVE_CYTHON:
        import platform

        # Cross-platform optimization flags based on existing project standards
        # Similar to what we use for Fortran compilation

        # Check compiler type on Windows
        is_msvc = platform.system() == "Windows" and (
            "MSVC" in os.environ.get("CC", "") or "cl.exe" in os.environ.get("CC", "")
        )

        # Check for Apple Silicon (arm64) architecture
        is_macos_arm64 = platform.system() == "Darwin" and platform.machine() == "arm64"

        if is_msvc:
            # MSVC flags for Windows
            extra_compile_args = [
                # Maximum speed optimization (stable, Microsoft recommended)
                "/O2",
                "/Oy",  # Frame pointer omission
                "/GT",  # Support fiber-safe thread-local storage
                # Use SSE2 instructions (widely supported on x86-64)
                "/arch:SSE2",
                # Note: Removed /fp:fast to preserve numerical precision
            ]
        elif is_macos_arm64:
            # Apple Silicon (arm64) optimized flags
            extra_compile_args = [
                "-O3",  # Maximum optimization
                "-march=armv8-a",  # ARM64 architecture
                "-mtune=apple-m1",  # Tune for Apple Silicon
                "-fPIC",  # Position Independent Code
                "-funroll-loops",  # Unroll loops for performance
                "-finline-functions",  # Inline functions
                "-ftree-vectorize",  # Enable vectorization
                "-ffinite-math-only",  # Assume finite math
                "-fno-trapping-math",  # Disable floating-point traps
                "-falign-functions=32",  # Function alignment
            ]
        else:
            # GCC/Clang compatible flags (Linux/x86-64 macOS/MinGW)
            # Using same strategy as Fortran compilation in this project
            extra_compile_args = [
                "-O3",  # Maximum optimization
                # Target x86-64 architecture (portable)
                "-march=x86-64",
                "-mtune=generic",  # Generic tuning (not CPU-specific)
                "-fPIC",  # Position Independent Code
                "-funroll-loops",  # Unroll loops for performance
                "-finline-functions",  # Inline functions
                "-ftree-vectorize",  # Enable vectorization
                # Assume finite math (same as Fortran config)
                "-ffinite-math-only",
                "-fno-trapping-math",  # Disable floating-point traps
                "-falign-functions=32",  # Function alignment
                # Note: Removed -ffast-math to preserve IEEE 754 compliance
            ]

        # Define the Cython extension for gridfill with optimizations
        gridfill_ext = Extension(
            "skyborn.gridfill._gridfill",
            ["src/skyborn/gridfill/_gridfill.pyx"],
            include_dirs=[np.get_include()],
            define_macros=[
                ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
                ("CYTHON_TRACE", "0"),  # Disable tracing for performance
                ("CYTHON_TRACE_NOGIL", "0"),  # Disable nogil tracing
            ],
            extra_compile_args=extra_compile_args,
            language="c",
        )
        extensions.append(gridfill_ext)

        compiler_type = "MSVC" if is_msvc else "GCC/Clang"
        print(f"Found Cython - will build cross-platform optimized gridfill extensions")
        print(f"Using {compiler_type} with flags: {extra_compile_args}")
    else:
        print(
            "Warning: Cython not found - gridfill Cython extensions will not be built"
        )
        print("Install Cython to enable gridfill functionality: pip install Cython")

    return extensions


class MesonBuildExt(build_ext):
    """Custom build extension to handle meson builds for Fortran modules"""

    def run(self):
        """Run the build process"""
        print("DEBUG: MesonBuildExt.run() called")
        # Build meson modules first
        self.build_meson_modules()
        # Then run the standard build_ext
        super().run()

    def build_meson_modules(self):
        """Build modules that use meson (like spharm)"""
        print("DEBUG: build_meson_modules() called")

        # Determine target directory based on --inplace flag
        if self.inplace:
            print("DEBUG: --inplace detected, building to source directory")
            spharm_target = Path("src") / "skyborn" / "spharm"
        else:
            print("DEBUG: Building to build directory")
            spharm_target = Path(self.build_lib) / "skyborn" / "spharm"

        meson_modules = [
            {
                "name": "spharm",
                "path": Path("src") / "skyborn" / "spharm",
                "target_dir": spharm_target,
            }
        ]

        for module in meson_modules:
            print(f"DEBUG: Processing module {module['name']}")
            if self.should_build_meson_module(module):
                print(f"DEBUG: Building module {module['name']} with meson")
                self.build_meson_module(module)
            else:
                print(f"DEBUG: Skipping module {module['name']} - no meson.build found")

    def should_build_meson_module(self, module):
        """Check if we should build this meson module"""
        meson_build_file = module["path"] / "meson.build"
        return meson_build_file.exists()

    def check_meson_available(self):
        """Check if meson and ninja are available"""
        try:
            # Check meson
            result = subprocess.run(
                ["meson", "--version"], capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                return False, "meson not found"

            meson_version = result.stdout.strip()
            print(f"Found meson version: {meson_version}")

            # Check ninja
            result = subprocess.run(
                ["ninja", "--version"], capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                return False, "ninja not found"

            ninja_version = result.stdout.strip()
            print(f"Found ninja version: {ninja_version}")

            return True, None

        except (
            subprocess.TimeoutExpired,
            subprocess.SubprocessError,
            FileNotFoundError,
        ) as e:
            return False, str(e)

    def build_meson_module(self, module):
        """
        Build a meson module using the meson build system.
        """
        print(f"Building {module['name']} with meson build system...")

        # Check if meson and ninja are available
        meson_available, error_msg = self.check_meson_available()
        if not meson_available:
            print(f"ERROR: Meson build tools not available: {error_msg}")
            print("Please install meson and ninja:")
            print("  pip install meson ninja")
            print("  or: conda install meson ninja")
            raise RuntimeError(
                f"Meson build tools required but not available: {error_msg}"
            )

        module_path = module["path"]
        # Use build subdirectory as specified in requirements
        build_dir = module_path / "build"

        try:
            # Clean build directory
            if build_dir.exists():
                print(f"Cleaning existing build directory: {build_dir}")
                shutil.rmtree(build_dir)

            # Setup build directory
            build_dir.mkdir(parents=True, exist_ok=True)

            # Configure meson build
            # Use local 'build' directory and '.' as source when running inside module_path
            # This avoids passing paths that meson will interpret relative to cwd and
            # causing doubled paths like src/skyborn/spharm/src/...
            print(f"Configuring meson build in {build_dir} (cwd={module_path})")
            setup_cmd = [
                "meson",
                "setup",
                "build",  # build directory inside module_path
                ".",  # source is current directory (module_path)
                "--buildtype=release",
                "-Db_lto=true",
            ]

            print(f"Running: {' '.join(setup_cmd)} (cwd={module_path})")

            # Set up environment for conda gfortran across all platforms
            env = os.environ.copy()
            import platform

            conda_prefix = env.get("CONDA_PREFIX", "")
            if conda_prefix:
                system = platform.system()
                current_path = env.get("PATH", "")

                if system == "Windows":
                    # Windows conda environment setup
                    conda_bin = os.path.join(conda_prefix, "bin")
                    conda_library_bin = os.path.join(conda_prefix, "Library", "bin")
                    env["PATH"] = f"{conda_bin};{conda_library_bin};{current_path}"
                    print(
                        f"Enhanced PATH for Windows conda environment: {conda_prefix}"
                    )

                elif system in ["Linux", "Darwin"]:
                    # Linux and macOS conda environment setup
                    conda_bin = os.path.join(conda_prefix, "bin")
                    # On Unix-like systems, use colon separator and prepend to PATH
                    env["PATH"] = f"{conda_bin}:{current_path}"

                    # Add lib directory to LD_LIBRARY_PATH (Linux) or DYLD_LIBRARY_PATH (macOS)
                    conda_lib = os.path.join(conda_prefix, "lib")
                    if system == "Linux":
                        current_lib_path = env.get("LD_LIBRARY_PATH", "")
                        env["LD_LIBRARY_PATH"] = (
                            f"{conda_lib}:{current_lib_path}"
                            if current_lib_path
                            else conda_lib
                        )
                        print(
                            f"Enhanced PATH and LD_LIBRARY_PATH for Linux conda environment: {conda_prefix}"
                        )
                    else:  # macOS
                        current_lib_path = env.get("DYLD_LIBRARY_PATH", "")
                        env["DYLD_LIBRARY_PATH"] = (
                            f"{conda_lib}:{current_lib_path}"
                            if current_lib_path
                            else conda_lib
                        )
                        print(
                            f"Enhanced PATH and DYLD_LIBRARY_PATH for macOS conda environment: {conda_prefix}"
                        )

                else:
                    print(
                        f"Warning: Unknown platform {system}, using basic conda PATH setup"
                    )
                    conda_bin = os.path.join(conda_prefix, "bin")
                    env["PATH"] = f"{conda_bin}:{current_path}"

            subprocess.run(setup_cmd, cwd=str(module_path), check=True, env=env)

            # Build with ninja (run relative to module_path, target 'build')
            print(f"Building with ninja in {build_dir} (cwd={module_path})")
            build_cmd = ["ninja", "-C", "build"]

            print(f"Running: {' '.join(build_cmd)} (cwd={module_path})")
            result = subprocess.run(
                build_cmd,
                cwd=str(module_path),
                check=True,
                capture_output=True,
                text=True,
            )

            if result.stdout:
                print("Build output:", result.stdout)
            if result.stderr:
                print("Build warnings/errors:", result.stderr)

            # File moving is now handled directly in meson.build

            print(f"Meson build for {module['name']} completed successfully!")

            self._built_modules = getattr(self, "_built_modules", set())
            self._built_modules.add(module["name"])

        except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError) as e:
            print(f"ERROR: Meson build failed for {module['name']}: {e}")
            if isinstance(e, subprocess.CalledProcessError):
                print(f"Command failed with exit code: {e.returncode}")
                if hasattr(e, "stdout") and e.stdout:
                    print("Stdout:", e.stdout)
                if hasattr(e, "stderr") and e.stderr:
                    print("Stderr:", e.stderr)
            raise  # Re-raise the exception since we're not using f2py fallback


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
    # Add extensions for both dummy (Windows compatibility) and gridfill
    "ext_modules": [
        Extension("skyborn._dummy", sources=["src/skyborn/_dummy.c"], optional=True)
    ]
    + (
        cythonize(get_gridfill_extensions())
        if HAVE_CYTHON and get_gridfill_extensions()
        else []
    ),
}

if __name__ == "__main__":
    setup(**setup_config)
