# Compiler detection and selection strategy
import os
import shutil
import subprocess
from pathlib import Path


class CompilerDetector:
    """Intelligent compiler detection and selection"""

    COMPILER_PRIORITY = {
        # Ordered by performance and compatibility
        "gfortran": {"priority": 1, "flags": ["-O3", "-fPIC"]},
        "ifort": {"priority": 2, "flags": ["-O3", "-fPIC", "-xHost"]},
        "ifx": {"priority": 3, "flags": ["-O3", "-fPIC"]},
        "flang": {"priority": 4, "flags": ["-O3", "-fPIC"]},
    }

    @classmethod
    def detect_available_compilers(cls):
        """Detect available Fortran compilers on the system"""
        available = {}

        for compiler, config in cls.COMPILER_PRIORITY.items():
            if shutil.which(compiler):
                try:
                    # Test if compiler works
                    result = subprocess.run(
                        [compiler, "--version"],
                        capture_output=True,
                        text=True,
                        timeout=10,
                    )
                    if result.returncode == 0:
                        available[compiler] = {
                            **config,
                            "version": result.stdout.split("\n")[0],
                            "path": shutil.which(compiler),
                        }
                except (subprocess.TimeoutExpired, subprocess.SubprocessError):
                    continue

        return available

    @classmethod
    def select_best_compiler(cls):
        """Select the best available compiler"""
        available = cls.detect_available_compilers()

        if not available:
            raise RuntimeError(
                "No Fortran compiler found. Please install gfortran:\n"
                "  Linux: sudo apt-get install gfortran\n"
                "  macOS: brew install gcc\n"
                "  Windows: conda install m2w64-toolchain"
            )

        # Sort by priority
        best_compiler = min(available.items(), key=lambda x: x[1]["priority"])
        return best_compiler[0], best_compiler[1]

    @classmethod
    def get_optimization_flags(cls, compiler_name, optimization_level="release"):
        """Get compiler optimization flags"""
        base_flags = cls.COMPILER_PRIORITY.get(compiler_name, {}).get("flags", ["-O2"])

        optimization_flags = {
            "debug": ["-g", "-O0"],
            "release": base_flags,
            "performance": base_flags + ["-march=native", "-funroll-loops"],
            "intel_optimized": (
                ["-O3", "-xHost", "-ipo"] if "ifort" in compiler_name else base_flags
            ),
        }

        return optimization_flags.get(optimization_level, base_flags)
