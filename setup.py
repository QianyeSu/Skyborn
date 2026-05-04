"""
Setup script for Skyborn - Mixed build system with meson for Fortran modules
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.install_lib import install_lib

# Check if Cython is available
try:
    from Cython.Build import cythonize

    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

# Check if we're in documentation build mode
DOCS_BUILD_MODE = (
    os.environ.get("SKYBORN_DOCS_BUILD") == "1" or os.environ.get("SKIP_FORTRAN") == "1"
)
CONDA_BUILD_MODE = os.environ.get("CONDA_BUILD") == "1"


def configure_compiler_environment():
    """Configure local compiler defaults without overriding conda-build."""
    if DOCS_BUILD_MODE:
        print("馃摎 Documentation build mode detected - skipping compiler setup")
        return

    if CONDA_BUILD_MODE:
        print("Conda-build environment detected - using externally provided compilers")
        return

    os.environ.setdefault("FC", "gfortran")
    os.environ.setdefault("F77", os.environ["FC"])
    os.environ.setdefault("F90", os.environ["FC"])
    os.environ.setdefault("CC", "gcc")
    print(
        "Using local default compilers: "
        f"FC={os.environ['FC']}, F77={os.environ['F77']}, "
        f"F90={os.environ['F90']}, CC={os.environ['CC']}"
    )


if DOCS_BUILD_MODE:
    print("📚 Documentation build mode detected - skipping Fortran compilation")
else:
    configure_compiler_environment()


# Check if the active Fortran compiler is available
def check_fortran_compiler():
    """Check whether the selected Fortran compiler can be executed."""
    compiler = os.environ.get("FC", "gfortran")
    try:
        result = subprocess.run(
            [compiler, "--version"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            print(
                f"Found Fortran compiler ({compiler}): "
                f"{result.stdout.split()[4] if len(result.stdout.split()) > 4 else 'unknown version'}"
            )
            return True
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
        pass

    print(
        f"Warning: Fortran compiler '{compiler}' not found. "
        "Fortran extensions may not build correctly."
    )
    print("Please install gfortran:")
    print("  Linux: sudo apt-get install gfortran")
    print("  macOS: brew install gcc")
    print("  Windows: conda install m2w64-toolchain")
    print("  Conda-build / conda-forge: rely on the recipe-provided compiler toolchain")
    return False


# Check gfortran availability at setup time (skip in docs mode)
if not DOCS_BUILD_MODE:
    check_fortran_compiler()
else:
    print("📚 Skipping gfortran check in documentation build mode")


# gridfill extensions now handled by meson.build in src/skyborn/gridfill/


class MesonBuildExt(build_ext):
    """Custom build extension to handle meson builds for Fortran modules"""

    def run(self):
        """Run the build process"""
        print("DEBUG: MesonBuildExt.run() called")
        # Build meson modules first
        self.build_meson_modules()
        # Then run the standard build_ext
        super().run()
        self.prune_installed_import_libraries()

    def prune_installed_import_libraries(self):
        """Remove Windows import libraries from the wheel/install staging area."""
        if self.inplace or not getattr(self, "build_lib", None):
            return

        build_lib_path = Path(self.build_lib)
        if not build_lib_path.exists():
            return

        removed = []
        for import_lib in build_lib_path.rglob("*.dll.a"):
            import_lib.unlink()
            removed.append(import_lib)

        if removed:
            print(
                "Removed Windows import libraries from staged package contents: "
                + ", ".join(str(path.relative_to(build_lib_path)) for path in removed)
            )

    def build_meson_modules(self):
        """Build modules that use meson (like spharm)"""
        print("DEBUG: build_meson_modules() called")

        # Skip meson builds in documentation mode
        if DOCS_BUILD_MODE:
            print("📚 Documentation build mode - skipping meson module builds")
            return

        # Determine target directory based on --inplace flag
        if self.inplace:
            print("DEBUG: --inplace detected, building to source directory")
            # spharm_target = Path("src") / "skyborn" / "spharm"
        else:
            print("DEBUG: Building to build directory")
            # spharm_target = Path(self.build_lib) / "skyborn" / "spharm"

        # Auto-discover meson modules based on directory structure
        # Each module should have a meson.build file
        meson_modules = self._discover_meson_modules()

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

            # Configure meson build with custom install directory for wheel builds
            print(f"Configuring meson build in {build_dir} (cwd={module_path})")

            setup_cmd = [
                "meson",
                "setup",
                "build",  # build directory inside module_path
                ".",  # source is current directory (module_path)
                "--buildtype=release",
            ]

            # Conda-forge toolchains already provide aggressive release flags.
            # Leaving Meson LTO enabled has caused two packaging-specific
            # failures:
            # 1. Linux static archives of Fortran helper objects were created
            #    with plain `ar`, which left unresolved symbols behind.
            # 2. macOS Meson configure mis-detected the linker under the
            #    conda-build toolchain environment.
            # Keep local editable/wheel builds on LTO, but turn it off under
            # conda-build for cross-platform package stability.
            if CONDA_BUILD_MODE:
                setup_cmd.append("-Db_lto=false")
            else:
                setup_cmd.append("-Db_lto=true")

            # For wheel builds, configure custom install directory
            if not self.inplace and hasattr(self, "build_lib") and self.build_lib:
                # Tell meson to install to our build directory instead of system
                build_lib_path = Path(self.build_lib).resolve()
                setup_cmd.extend(
                    [
                        f"--python.purelibdir={build_lib_path}",
                        f"--python.platlibdir={build_lib_path}",
                    ]
                )
                print(f"DEBUG: Configuring meson to install to: {build_lib_path}")

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
                    mingw_bin = os.path.join(
                        conda_prefix, "Library", "mingw-w64", "bin"
                    )
                    usr_bin = os.path.join(conda_prefix, "Library", "usr", "bin")
                    path_entries = [conda_bin, conda_library_bin, mingw_bin, usr_bin]
                    env["PATH"] = ";".join(path_entries + [current_path])
                    print(
                        f"Enhanced PATH for Windows conda environment: {conda_prefix}"
                    )

                    if CONDA_BUILD_MODE:
                        # conda-build activation can still point CC/CXX at MSVC on
                        # Windows even when gcc/gfortran packages are installed.
                        # Meson then mixes `cl/link` with MinGW Fortran objects and
                        # fails to link `gfortran.lib`. Keep the Meson subprocess on
                        # one MinGW toolchain so mixed-language extensions link
                        # consistently under conda-forge.
                        toolchain = {
                            "CC": "x86_64-w64-mingw32-gcc.exe",
                            "CXX": "x86_64-w64-mingw32-g++.exe",
                            "FC": "x86_64-w64-mingw32-gfortran.exe",
                            "F77": "x86_64-w64-mingw32-gfortran.exe",
                            "F90": "x86_64-w64-mingw32-gfortran.exe",
                            "AR": "x86_64-w64-mingw32-ar.exe",
                            "RANLIB": "x86_64-w64-mingw32-ranlib.exe",
                        }
                        resolved = {}
                        for key, exe in toolchain.items():
                            path = shutil.which(exe, path=env["PATH"])
                            if path:
                                env[key] = path
                                resolved[key] = path

                        if {"CC", "CXX", "FC"} <= resolved.keys():
                            print(
                                "Overriding Windows conda-build compiler environment for Meson: "
                                + ", ".join(
                                    f"{key}={resolved[key]}"
                                    for key in ("CC", "CXX", "FC")
                                )
                            )
                        else:
                            print(
                                "Warning: could not fully resolve MinGW compilers for Meson on Windows; "
                                "keeping existing compiler environment"
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

                        if CONDA_BUILD_MODE:
                            # Meson on macOS + gfortran may probe for GNU-style tool
                            # names such as `gcc-ar` and plain `ar`. Exporting `AR`
                            # directly causes Meson to treat it as a linker candidate,
                            # so instead provide tiny PATH shims with the names Meson
                            # expects and let its default discovery succeed.
                            import tempfile

                            # Prefer the generic LLVM/binutils programs shipped by
                            # conda-build. The prefixed cross-tool wrappers can
                            # exist but still confuse Meson's `--version`
                            # detection on macOS when a mixed C/Fortran project
                            # probes the static archiver during `project()`.
                            ar_candidates = [
                                os.path.join(conda_bin, "llvm-ar"),
                                os.path.join(conda_bin, "ar"),
                            ]
                            ranlib_candidates = [
                                os.path.join(conda_bin, "llvm-ranlib"),
                                os.path.join(conda_bin, "ranlib"),
                            ]

                            for compiler_key in ("CC", "CXX", "FC", "F77", "F90"):
                                compiler_path = env.get(compiler_key)
                                if not compiler_path:
                                    continue

                                compiler_name = os.path.basename(compiler_path)
                                compiler_dir = os.path.dirname(compiler_path)
                                if compiler_name.endswith(
                                    ("clang", "clang++", "gfortran", "gcc", "g++")
                                ):
                                    prefix = compiler_name.rsplit("-", 1)[0] + "-"
                                    ar_candidates.extend(
                                        [
                                            os.path.join(compiler_dir, prefix + "ar"),
                                            os.path.join(
                                                compiler_dir, prefix + "gcc-ar"
                                            ),
                                        ]
                                    )
                                    ranlib_candidates.extend(
                                        [
                                            os.path.join(
                                                compiler_dir, prefix + "ranlib"
                                            ),
                                            os.path.join(
                                                compiler_dir, prefix + "gcc-ranlib"
                                            ),
                                        ]
                                    )

                            ar_candidates.append("/usr/bin/ar")
                            ranlib_candidates.append("/usr/bin/ranlib")

                            def _pick_first_executable(candidates):
                                seen = set()
                                for candidate in candidates:
                                    if not candidate or candidate in seen:
                                        continue
                                    seen.add(candidate)
                                    if os.path.exists(candidate) and os.access(
                                        candidate, os.X_OK
                                    ):
                                        return candidate
                                return None

                            ar_path = _pick_first_executable(ar_candidates)
                            ranlib_path = _pick_first_executable(ranlib_candidates)

                            if ar_path or ranlib_path:
                                shim_dir = tempfile.mkdtemp(
                                    prefix="skyborn-meson-tools-"
                                )

                                def _write_tool_shim(name, target, version_banner=None):
                                    if not target:
                                        return None
                                    shim_path = os.path.join(shim_dir, name)
                                    with open(
                                        shim_path, "w", encoding="utf-8", newline="\n"
                                    ) as handle:
                                        handle.write("#!/bin/sh\n")
                                        if version_banner:
                                            handle.write(
                                                'if [ "x$1" = "x--version" ]; then\n'
                                            )
                                            handle.write(f'  echo "{version_banner}"\n')
                                            handle.write("  exit 0\n")
                                            handle.write("fi\n")
                                        handle.write(f'exec "{target}" "$@"\n')
                                    os.chmod(shim_path, 0o755)
                                    return shim_path

                                created = []
                                for shim_name in ("ar", "gcc-ar", "gar"):
                                    shim_path = _write_tool_shim(
                                        shim_name,
                                        ar_path,
                                        "GNU ar (Skyborn Meson shim)",
                                    )
                                    if shim_path:
                                        created.append(f"{shim_name}={shim_path}")

                                for shim_name in ("ranlib", "gcc-ranlib"):
                                    shim_path = _write_tool_shim(
                                        shim_name,
                                        ranlib_path,
                                        "GNU ranlib (Skyborn Meson shim)",
                                    )
                                    if shim_path:
                                        created.append(f"{shim_name}={shim_path}")

                                env["PATH"] = f"{shim_dir}:{env['PATH']}"
                                print(
                                    "Prepared macOS Meson tool shims: "
                                    + ", ".join(created)
                                )
                                print(
                                    "Selected macOS archiver tools for Meson: "
                                    f"ar_target={ar_path}, ranlib_target={ranlib_path}"
                                )

                else:
                    print(
                        f"Warning: Unknown platform {system}, using basic conda PATH setup"
                    )
                    conda_bin = os.path.join(conda_prefix, "bin")
                    env["PATH"] = f"{conda_bin}:{current_path}"

                # Meson should detect its own linker / archiver from the
                # compiler driver. On conda-build macOS these environment
                # variables can cause Meson to treat the archiver itself as the
                # linker during configure.
                if CONDA_BUILD_MODE and system == "Darwin":
                    for key in (
                        "AR",
                        "RANLIB",
                        "GCC_AR",
                        "GCC_RANLIB",
                        "LD",
                        "LDSHARED",
                        "CC_LD",
                        "CXX_LD",
                        "FC_LD",
                        "F77_LD",
                        "F90_LD",
                    ):
                        if key in env:
                            env.pop(key, None)

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

            # Install using meson (this will handle the path configuration we set up)
            if not self.inplace and hasattr(self, "build_lib") and self.build_lib:
                print(f"Installing meson build outputs to {self.build_lib}")
                install_cmd = ["meson", "install", "-C", "build", "--only-changed"]
                print(f"Running: {' '.join(install_cmd)} (cwd={module_path})")

                try:
                    install_result = subprocess.run(
                        install_cmd,
                        cwd=str(module_path),
                        check=True,
                        capture_output=True,
                        text=True,
                        env=env,
                    )
                    if install_result.stdout:
                        print("Install output:", install_result.stdout)
                except subprocess.CalledProcessError as e:
                    print(f"ERROR: Meson install failed: {e}")
                    print(
                        "This should not happen with the new setup. Please check meson configuration."
                    )
                    raise

                self.prune_installed_import_libraries()
            else:
                print("Inplace build - extensions handled by meson custom_target")

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

    def _discover_meson_modules(self):
        """
        Auto-discover meson modules by recursively looking for meson.build files
        in skyborn subpackages and their subdirectories.
        """
        modules = []
        skyborn_src = Path("src") / "skyborn"

        def _find_meson_builds(base_path, relative_path=""):
            """Recursively find meson.build files"""
            for subdir in base_path.iterdir():
                if subdir.is_dir():
                    meson_file = subdir / "meson.build"
                    if meson_file.exists():
                        if relative_path:
                            module_name = f"{relative_path}.{subdir.name}"
                        else:
                            module_name = subdir.name
                        print(f"DEBUG: Discovered meson module: {module_name}")
                        modules.append(
                            {
                                "name": module_name,
                                "path": subdir,
                            }
                        )
                    # Continue searching in subdirectories
                    new_relative_path = (
                        f"{relative_path}.{subdir.name}"
                        if relative_path
                        else subdir.name
                    )
                    _find_meson_builds(subdir, new_relative_path)

        # Start recursive search from skyborn root
        _find_meson_builds(skyborn_src)

        return modules


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


class CustomInstallLib(install_lib):
    """Prune transient build artifacts and bytecode from installed package trees."""

    _PRUNE_DIR_NAMES = {
        "__pycache__",
        "build",
        "build_local",
        "meson-private",
        "meson-logs",
        "meson-info",
    }
    _PRUNE_FILE_SUFFIXES = {".pyc", ".pyo"}
    _PRUNE_RELATIVE_DIRS = {
        Path("spharm") / "src",
        Path("gridfill") / "tests",
    }
    _PRUNE_RELATIVE_FILES = {
        Path("spharm") / "fix_f2py_symbols.py",
    }

    def run(self):
        super().run()
        self._prune_install_tree()

    def _prune_install_tree(self):
        install_root = Path(self.install_dir)
        package_root = install_root / "skyborn"
        if not package_root.exists():
            return

        pruned_dirs = []
        for relative_dir in self._PRUNE_RELATIVE_DIRS:
            target = package_root / relative_dir
            if target.exists():
                shutil.rmtree(target, ignore_errors=True)
                pruned_dirs.append(target)

        for path in sorted(
            package_root.rglob("*"), key=lambda item: len(item.parts), reverse=True
        ):
            if path.is_dir() and path.name in self._PRUNE_DIR_NAMES:
                shutil.rmtree(path, ignore_errors=True)
                pruned_dirs.append(path)

        pruned_files = []
        for relative_file in self._PRUNE_RELATIVE_FILES:
            target = package_root / relative_file
            if target.exists():
                target.unlink()
                pruned_files.append(target)

        for path in package_root.rglob("*"):
            if path.is_file() and path.suffix.lower() in self._PRUNE_FILE_SUFFIXES:
                path.unlink()
                pruned_files.append(path)

        if pruned_dirs or pruned_files:
            print(
                "Pruned non-runtime packaging artifacts: "
                f"{len(pruned_dirs)} directories, {len(pruned_files)} files"
            )


# Configuration for mixed build
setup_config = {
    "cmdclass": {
        "build_ext": MesonBuildExt,
        "develop": CustomDevelop,
        "install": CustomInstall,
        "install_lib": CustomInstallLib,
    },
    # Add extensions for dummy (Windows compatibility) only
    # gridfill extensions now handled by meson.build
    "ext_modules": [
        Extension("skyborn._dummy", sources=["src/skyborn/_dummy.c"], optional=True)
    ],
}

if __name__ == "__main__":
    setup(**setup_config)
