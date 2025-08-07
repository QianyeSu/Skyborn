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
        meson_modules = [
            {
                "name": "spharm",
                "path": Path("src") / "skyborn" / "spharm",
                "target_dir": Path(self.build_lib) / "skyborn" / "spharm",
            }
        ]

        for module in meson_modules:
            print(f"DEBUG: Processing module {module['name']}")
            if self.should_build_meson_module(module):
                print(f"DEBUG: Building module {module['name']}")
                self.build_single_meson_module(module)
            else:
                print(f"DEBUG: Skipping module {module['name']} - no meson.build found")

    def should_build_meson_module(self, module):
        """Check if we should build this meson module"""
        meson_build_file = module["path"] / "meson.build"
        return meson_build_file.exists()

    def build_single_meson_module(self, module):
        """
        Build a single meson module using a corrected two-step f2py process.
        """
        print(f"Building {module['name']} with two-step f2py process...")

        module_path = module["path"]
        # CRITICAL: Use an isolated build directory for all generated files.
        # This prevents polluting the source tree and avoids permission errors.
        build_temp = module_path / "build"

        try:
            # Clean and create the isolated build directory
            if build_temp.exists():
                shutil.rmtree(build_temp)
            build_temp.mkdir(parents=True, exist_ok=True)

            src_dir = module_path / "src"
            f_files = list(src_dir.glob("*.f"))
            # f_files = [
            #     src_dir / "lap.f90",
            #     src_dir / "invlap.f90",
            #     src_dir / "sphcom.f",
            #     src_dir / "hrfft.f",
            #     src_dir / "getlegfunc.f",
            #     src_dir / "specintrp.f",
            #     src_dir / "onedtotwod.f",
            #     src_dir / "onedtotwod_vrtdiv.f",
            #     src_dir / "twodtooned.f",
            #     src_dir / "twodtooned_vrtdiv.f",
            #     src_dir / "multsmoothfact.f",
            #     src_dir / "gaqd.f",
            #     src_dir / "shses.f",
            #     src_dir / "shaes.f",
            #     src_dir / "vhaes.f",
            #     src_dir / "vhses.f",
            #     src_dir / "shsgs.f",
            #     src_dir / "shags.f",
            #     src_dir / "vhags.f",
            #     src_dir / "vhsgs.f",
            #     src_dir / "shaec.f",
            #     src_dir / "shagc.f",
            #     src_dir / "shsec.f",
            #     src_dir / "shsgc.f",
            #     src_dir / "vhaec.f",
            #     src_dir / "vhagc.f",
            #     src_dir / "vhsec.f",
            #     src_dir / "vhsgc.f",
            #     src_dir / "ihgeod.f",
            #     src_dir / "alf.f",
            # ]
            pyf_file = src_dir / "_spherepack.pyf"

            # --- STEP 1: Generate wrapper files in the isolated build directory ---
            print("Step 1: Generating C and Fortran wrapper files...")

            generate_cmd = [
                sys.executable,
                "-m",
                "numpy.f2py",
                str(pyf_file),
                "--lower",
                # CRITICAL: Tell f2py the name of the module. This ensures
                # '_spherepackmodule.c' is created with the correct content.
                "-m",
                "_spherepack",
            ]

            print(f"Running generation command: {' '.join(generate_cmd)}")
            # CRITICAL: Run the command from the project root directory.
            subprocess.run(generate_cmd, cwd=str(Path.cwd()), check=True)

            # --- STEP 2: Verify generated files and build the full source list ---
            print("Step 2: Verifying wrappers and preparing source list...")

            # Generated files will be in the current working directory (project root)
            generated_c_wrapper = Path.cwd() / "_spherepackmodule.c"
            generated_f_wrapper = Path.cwd() / "_spherepack-f2pywrappers.f"

            # CRITICAL: Only the C wrapper is mandatory.
            if not generated_c_wrapper.exists():
                raise RuntimeError(
                    f"f2py C wrapper not generated! Looked for: {generated_c_wrapper}"
                )

            print(f"Found mandatory C wrapper: {generated_c_wrapper.name}")

            # Build the final list of source files for the compiler
            compile_sources = [str(f) for f in f_files] + [str(generated_c_wrapper)]

            # Add the optional Fortran wrapper only if it exists
            if generated_f_wrapper.exists():
                print(f"Found optional Fortran wrapper: {generated_f_wrapper.name}")
                compile_sources.append(str(generated_f_wrapper))

            # --- STEP 3: Compile all sources together ---
            print("Step 3: Compiling all sources...")

            # The .pyf file must still be the first argument for f2py to get metadata
            f2py_cmd = (
                [sys.executable, "-m", "numpy.f2py", "-c", str(pyf_file)]
                + compile_sources
                # Use setuptools' global temp dir for output
                + ["--build-dir", str(self.build_temp)]
            )

            # fortran_optim_flags = "-O3 -fPIC -fno-second-underscore -funroll-loops -finline-functions -ftree-vectorize -ffinite-math-only"
            # fortran_optim_flags = "-O3 -fPIC -fno-second-underscore -funroll-loops -finline-functions -ftree-vectorize -ffinite-math-only"
            import platform

            # Check for Apple Silicon (arm64) architecture
            is_macos_arm64 = (
                platform.system() == "Darwin" and platform.machine() == "arm64"
            )

            if is_macos_arm64:
                # Apple Silicon (arm64) optimized flags
                fortran_optim_flags = (
                    "-O3 "
                    "-fPIC "
                    "-fno-second-underscore "
                    "-funroll-loops "
                    "-finline-functions "
                    "-ftree-vectorize "
                    "-ffinite-math-only "
                    "-march=armv8-a "
                    "-mtune=apple-m1 "
                    "-fno-common "
                    "-ftree-loop-im "
                    "-ftree-loop-distribution "
                    "-falign-functions=32"
                )
                c_optim_flags = (
                    "-O3 "
                    "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION "
                    "-march=armv8-a "
                    "-mtune=apple-m1 "
                    "-fPIC "
                    "-fno-trapping-math "
                    "-falign-functions=32"
                )
            else:
                # x86-64 optimized flags (Linux/Windows/Intel macOS)
                fortran_optim_flags = (
                    "-O3 "
                    "-fPIC "
                    "-fno-second-underscore "
                    "-funroll-loops "
                    "-finline-functions "
                    "-ftree-vectorize "
                    "-ffinite-math-only "
                    "-march=x86-64 "
                    "-mtune=generic "
                    "-fno-common "
                    "-ftree-loop-im "
                    "-ftree-loop-distribution "
                    "-falign-functions=32"
                )
                c_optim_flags = (
                    # Highest level optimization, covers most performance enhancements
                    "-O3 "
                    # NumPy API compatibility, very important
                    "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION "
                    # Explicitly target all 64-bit x86 CPUs for optimization, consistent with Fortran, ensures portability
                    "-march=x86-64 "
                    # Tune for generic processors, not specific models, consistent with Fortran
                    "-mtune=generic "
                    # Position Independent Code, required for shared libraries
                    "-fPIC "
                    # Disable floating-point exception traps, prevents crashes due to floating-point errors, enhances robustness
                    "-fno-trapping-math "
                    # Align function entry points, potentially slightly improves instruction cache efficiency
                    "-falign-functions=32"
                )
            # f2py_cmd += [
            #     "--opt=" + fortran_optim_flags,
            #     "--f90flags=" + fortran_optim_flags,
            #     "--cflags=" + c_optim_flags,
            # ]
            print(
                "INFO: Setting compiler flags via environment variables for the Meson backend."
            )
            build_env = os.environ.copy()
            build_env["FFLAGS"] = fortran_optim_flags
            build_env["CFLAGS"] = c_optim_flags
            import platform

            if platform.system() == "Windows":
                f2py_cmd += ["--fcompiler=gnu95"]
            else:
                f2py_cmd += ["--fcompiler=gnu95", "--compiler=unix"]

            print("f2py build command:", " ".join(f2py_cmd))
            print(f"Using FFLAGS: {build_env.get('FFLAGS')}")
            print(f"Using CFLAGS: {build_env.get('CFLAGS')}")
            subprocess.run(
                f2py_cmd,
                # Run from project root so that build dir paths are correct
                cwd=str(Path.cwd()),
                check=True,
                env=build_env,  # <-- This passes the environment to the subprocess
            )

            # --- STEP 4: Clean up generated files and move compiled file ---
            print("Step 4: Cleaning up generated files and moving compiled module...")

            # Clean up generated wrapper files from project root
            if generated_c_wrapper.exists():
                generated_c_wrapper.unlink()
            if generated_f_wrapper.exists():
                generated_f_wrapper.unlink()

            # Find the compiled file (e.g., _spherepack.cpython-312-x86_64-linux-gnu.so)
            # f2py outputs to project root, not to build_temp
            compiled_files = (
                list(Path.cwd().glob("_spherepack*.so"))
                + list(Path.cwd().glob("_spherepack*.pyd"))
                + list(Path.cwd().glob("_spherepack*.dylib"))
            )
            if not compiled_files:
                # Also check build_temp in case f2py behavior changes
                compiled_files = (
                    list(Path(self.build_temp).glob("_spherepack*.so"))
                    + list(Path(self.build_temp).glob("_spherepack*.pyd"))
                    + list(Path(self.build_temp).glob("_spherepack*.dylib"))
                )

            if not compiled_files:
                raise FileNotFoundError(
                    f"Could not find compiled module in project root or {self.build_temp}"
                )

            source_file = compiled_files[0]
            target_dir = module["target_dir"]
            target_dir.mkdir(parents=True, exist_ok=True)

            print(f"Moving {source_file} to {target_dir}")
            shutil.move(str(source_file), str(target_dir))

            print(f"spharm compilation successful!")
            self._built_modules = getattr(self, "_built_modules", set())
            self._built_modules.add(module["name"])

        except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError) as e:
            print(f"ERROR: Failed to build {module['name']}: {e}")
            # Print more info on subprocess errors
            if isinstance(e, subprocess.CalledProcessError):
                print("Stdout:", e.stdout)
                print("Stderr:", e.stderr)
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
