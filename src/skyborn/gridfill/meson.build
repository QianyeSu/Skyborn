# Skyborn gridfill submodule - Cython-based grid filling with cross-platform optimization
project('skyborn-gridfill', 'c',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

# Get Python installation
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Install Python module files
py.install_sources(
    ['__init__.py', 'gridfill.py', 'xarray.py'],
    subdir: 'skyborn/gridfill',
)

# Check for Cython availability
cython_check = run_command(py, '-c', 'import Cython; print("OK")', check: false)
have_cython = cython_check.returncode() == 0

if have_cython
    # Get numpy include directories
    incdir_numpy_cmd = run_command(
        py,
        ['-c', '''
import numpy
import os
# Try multiple potential paths
paths = [
    numpy.get_include(),
    os.path.join(numpy.__path__[0], "core", "include"),
    os.path.join(numpy.__path__[0], "_core", "include")
]
for path in paths:
    if os.path.exists(os.path.join(path, "numpy", "npy_os.h")):
        print(path)
        break
else:
    print(numpy.get_include())  # fallback
'''],
        check: true,
    )
    incdir_numpy = incdir_numpy_cmd.stdout().strip()

    message('NumPy include directory: ' + incdir_numpy)

    # Detect platform and architecture for optimization flags
    host_system = host_machine.system()
    host_cpu = host_machine.cpu_family()

    message('Host system: ' + host_system)
    message('Host CPU family: ' + host_cpu)

    # Cross-platform optimization flags similar to setup.py
    c_args = [
        '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
        '-DCYTHON_TRACE=0',
        '-DCYTHON_TRACE_NOGIL=0',
        '-fPIC',
        # Suppress pointer type warnings that occur on some platforms
        # due to differences in how np.uint32_t is defined
        '-Wno-incompatible-pointer-types',
    ]

    # Use platform-based detection like spharm does
    if host_system == 'windows'
        # Only use MSVC flags if we're actually using MSVC compiler
        cc = meson.get_compiler('c')
        if cc.get_id() == 'msvc'
            c_args += [
                '/O2',     # Maximum optimization
                '/GL',     # Whole program optimization
                '/Oi',     # Enable intrinsic functions
                '/Ot',     # Favor fast code
            ]
        else
            # Windows with GCC/MinGW - use GCC flags
            c_args += [
                '-O3',                    # Maximum optimization
                '-march=x86-64',          # Target x86-64 architecture (portable)
                '-mtune=generic',         # Generic tuning (not CPU-specific)
                '-funroll-loops',         # Unroll loops for performance
                '-finline-functions',     # Inline functions
                '-ftree-vectorize',       # Enable vectorization
                '-ffinite-math-only',     # Assume finite math
                '-fno-trapping-math',     # Disable floating-point traps
                '-falign-functions=32',   # Function alignment
            ]
        endif
    elif host_system == 'darwin' and host_cpu == 'aarch64'
        # Apple Silicon (arm64) optimized flags
        message('Configuring for Apple Silicon (arm64)')
        c_args += [
            '-O3',                    # Maximum optimization
            '-march=armv8-a',         # ARM64 architecture
            '-mtune=apple-m1',        # Tune for Apple Silicon
            '-funroll-loops',         # Unroll loops for performance
            '-finline-functions',     # Inline functions
            '-ftree-vectorize',       # Enable vectorization
            '-ffinite-math-only',     # Assume finite math
            '-fno-trapping-math',     # Disable floating-point traps
            '-falign-functions=32',   # Function alignment
        ]
    else
        # GCC/Clang compatible flags (Linux/x86-64 macOS/MinGW)
        message('Configuring for x86-64 architecture')
        c_args += [
            '-O3',                    # Maximum optimization
            '-march=x86-64',          # Target x86-64 architecture (portable)
            '-mtune=generic',         # Generic tuning (not CPU-specific)
            '-funroll-loops',         # Unroll loops for performance
            '-finline-functions',     # Inline functions
            '-ftree-vectorize',       # Enable vectorization
            '-ffinite-math-only',     # Assume finite math (same as setup.py)
            '-fno-trapping-math',     # Disable floating-point traps
            '-falign-functions=32',   # Function alignment
        ]
    endif

    message('C args: ' + ' '.join(c_args))

    inc_np = include_directories(incdir_numpy)

    # Check for Cython file and build extension
    cython_source_file = '_gridfill.pyx'
    cython_source_check = run_command(
        py,
        ['-c', 'import os; import sys; sys.exit(0 if os.path.exists("' + cython_source_file + '") else 1)'],
        check: false,
    )

    if cython_source_check.returncode() == 0
        message('Found Cython source: ' + cython_source_file)

        # Generate C source from Cython
        cython_c_target = custom_target(
            '_gridfill.c',
            input: cython_source_file,
            output: '_gridfill.c',
            command: [
                py, '-m', 'cython',
                '--3str',  # Python 3 string handling
                '--fast-fail',  # Fail fast on errors
                '--directive', 'language_level=3',  # Python 3 language level
                '--directive', 'boundscheck=False',  # Disable bounds checking for performance
                '--directive', 'wraparound=False',   # Disable negative index wrapping
                '--directive', 'cdivision=True',     # C-style division
                '--directive', 'nonecheck=False',    # Disable None checks
                '--directive', 'profile=False',      # Disable profiling
                '--directive', 'linetrace=False',    # Disable line tracing
                '--output-file', '@OUTPUT@',
                '@INPUT@'
            ],
            build_by_default: true,
            install: false,
        )

        # Build the gridfill extension module
        gridfill_ext = py.extension_module(
            '_gridfill',
            cython_c_target,
            include_directories: inc_np,
            dependencies: [py_dep],
            c_args: c_args,
            install: true,
            install_dir: get_option('python.purelibdir') / 'skyborn' / 'gridfill',
            build_by_default: true,
        )

        # Copy to source directory for --inplace builds
        copy_to_source = custom_target(
            'copy_to_source_for_inplace',
            input: gridfill_ext,
            output: 'inplace_copy_marker.txt',
            command: [
                py,
                '-c',
                '''
import os
import shutil
import glob
import sys
from pathlib import Path

print("=== Gridfill extension source copy check ===")

current_dir = Path.cwd()  # This is the build directory
module_dir = current_dir.parent  # This is the module source directory
print(f"Current dir (build): {current_dir}")
print(f"Module dir (source): {module_dir}")

# Search for compiled extension files
patterns = ["_gridfill*.so", "_gridfill*.pyd", "_gridfill*.dylib"]
built_files = []

# Search directly in the build directory
for pattern in patterns:
    found_files = glob.glob(str(current_dir / pattern), recursive=False)
    built_files.extend(found_files)

# Also search recursively
for pattern in patterns:
    found_files = glob.glob(str(current_dir / "**" / pattern), recursive=True)
    built_files.extend(found_files)

# Remove duplicates
built_files = list(set(built_files))

if built_files:
    print(f"Found built files: {built_files}")

    for built_file_str in built_files:
        built_file = Path(built_file_str)
        source_file = module_dir / built_file.name

        # Copy if source doesn't exist or is older
        should_copy = False
        if not source_file.exists():
            should_copy = True
            print(f"Source {built_file.name} doesn't exist, will copy")
        elif built_file.stat().st_mtime > source_file.stat().st_mtime:
            should_copy = True
            print(f"Source {built_file.name} is older, will copy")

        if should_copy:
            try:
                if source_file.exists():
                    source_file.unlink()
                shutil.copy2(str(built_file), str(source_file))
                print(f"Copied {built_file.name} to source directory successfully")
            except Exception as e:
                print(f"Failed to copy {built_file.name}: {e}")
else:
    print("No compiled extension files found")

# Create marker
with open("''' + '@OUTPUT@' + '''", "w") as f:
    f.write("Gridfill source copy check completed\\n")

print("=== Gridfill extension source copy check completed ===")
'''
            ],
            build_by_default: true,
            console: true
        )

        message('gridfill submodule configured successfully with Cython extensions')
    else
        warning('Cython source file not found: ' + cython_source_file)
        message('Installing Python-only gridfill module')
    endif

else
    warning('Cython not available, skipping gridfill extension generation')
    message('Installing Python-only gridfill module')
endif
