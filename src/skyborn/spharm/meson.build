# Skyborn spharm submodule - Spherical harmonic transforms with cross-platform optimization
project('skyborn-spharm', 'c', 'fortran',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

# Get Python installation
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Find OpenMP dependency
openmp_dep = dependency('openmp', required: false)

# Install Python module files
py.install_sources(
    ['__init__.py', 'spherical_harmonics.py'],
    subdir: 'skyborn/spharm',
)

# Get numpy include directories with robust fallbacks
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

# Try to get f2py include, with multiple fallbacks
incdir_f2py_cmd = run_command(
    py,
    ['-c', '''
import numpy
import os
try:
    import numpy.f2py
    f2py_path = numpy.f2py.get_include()
    if os.path.exists(f2py_path):
        print(f2py_path)
    else:
        print(numpy.get_include())
except:
    print(numpy.get_include())
'''],
    check: true,
)
incdir_f2py = incdir_f2py_cmd.stdout().strip()

message('NumPy include directory: ' + incdir_numpy)
message('F2PY include directory: ' + incdir_f2py)

# Detect platform and architecture for optimization flags
host_system = host_machine.system()
host_cpu = host_machine.cpu_family()

message('Host system: ' + host_system)
message('Host CPU family: ' + host_cpu)

# Cross-platform optimization flags based on setup.py configuration
fortran_args = []
c_args = []

if host_system == 'darwin' and host_cpu == 'aarch64'
    # Apple Silicon (arm64) optimized flags
    message('Configuring for Apple Silicon (arm64)')
    fortran_args = [
        '-O3',                      # Maximum optimization
        '-fPIC',                    # Position independent code
        '-fno-second-underscore',   # Fortran symbol naming convention
        '-funroll-loops',           # Loop unrolling optimization
        '-finline-functions',       # Inline small functions
        '-ftree-vectorize',         # Auto-vectorization for SIMD
        '-ffinite-math-only',       # Assume no NaN/Inf values
        '-ffree-line-length-none',  # Allow arbitrarily long free-form source lines
        '-march=armv8-a',           # ARM64 architecture
        '-mtune=apple-m1',          # Tune for Apple Silicon
        '-fno-common',              # Don't use common blocks
        '-ftree-loop-im',           # Loop invariant motion
        '-ftree-loop-distribution', # Loop distribution
        '-falign-functions=32',     # Function alignment
        '-falign-loops=32',         # Loop alignment optimization
        '-std=legacy',              # Legacy Fortran standard
        '-fopenmp',
    ]
    c_args = [
        # '-O3',                                        # Maximum optimization
        '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', # NumPy API compatibility
        '-march=armv8-a',                             # ARM64 architecture
        '-mtune=apple-m1',                            # Tune for Apple Silicon
        '-fPIC',                                      # Position independent code
        '-fno-trapping-math',                         # Disable floating-point traps
        '-falign-functions=32',                       # Function alignment
        '-fopenmp',
    ]
else
    # x86-64 optimized flags (Linux/Windows/Intel macOS)
    message('Configuring for x86-64 architecture')
    fortran_args = [
        '-O3',                      # Maximum optimization
        '-fPIC',                    # Position independent code
        '-fno-second-underscore',   # Fortran symbol naming convention
        '-funroll-loops',           # Loop unrolling optimization
        '-finline-functions',       # Inline small functions
        '-ftree-vectorize',         # Auto-vectorization for SIMD
        '-ffinite-math-only',       # Assume no NaN/Inf values
        '-ffree-line-length-none',  # Allow arbitrarily long free-form source lines
        '-march=x86-64',            # x86-64 architecture
        '-mtune=generic',           # Generic tuning (portable)
        '-fno-common',              # Don't use common blocks
        '-ftree-loop-im',           # Loop invariant motion
        '-ftree-loop-distribution', # Loop distribution
        '-falign-functions=32',     # Function alignment
        '-falign-loops=32',         # Loop alignment optimization
        '-std=legacy',              # Legacy Fortran standard
        '-fopenmp',
    ]
    c_args = [
        # '-O3',                                        # Maximum optimization
        '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',# NumPy API compatibility
        '-march=x86-64',                              # x86-64 architecture
        '-mtune=generic',                             # Generic tuning (portable)
        '-fPIC',                                      # Position independent code
        '-fno-trapping-math',                         # Disable floating-point traps
        '-falign-functions=32',                       # Function alignment
        '-fopenmp',
    ]
endif

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

# Generate Fortran-Python interface with f2py
try_f2py = run_command(py, '-c', 'import numpy.f2py; print("OK")', check: false)

if try_f2py.returncode() == 0
    inc_np = include_directories(incdir_numpy, incdir_f2py)

    # Verify numpy headers exist
    numpy_header_check = run_command(
        py,
        ['-c', 'import os; import numpy; header_path = os.path.join(numpy.get_include(), "numpy", "npy_os.h"); print("✓" if os.path.exists(header_path) else "✗", header_path)'],
        check: false,
    )
    message('NumPy header check: ' + numpy_header_check.stdout().strip())

    # Get fortranobject.c path with extensive fallbacks
    fortranobject_cmd = run_command(
        py,
        ['-c', '''
import numpy.f2py
import numpy
import os
import sys

# Try multiple potential paths for fortranobject.c
candidates = []

# Method 1: f2py.get_include()
try:
    fo_path = os.path.join(numpy.f2py.get_include(), "fortranobject.c")
    candidates.append(fo_path)
except:
    pass

# Method 2: numpy package paths
numpy_base = numpy.__path__[0]
candidates.extend([
    os.path.join(numpy_base, "f2py", "src", "fortranobject.c"),
    os.path.join(numpy_base, "f2py", "fortranobject.c"),
    os.path.join(numpy_base, "_core", "src", "fortranobject.c"),
])

# Method 3: site-packages paths
try:
    import site
    for site_dir in site.getsitepackages() + [site.getusersitepackages()]:
        if site_dir:
            candidates.extend([
                os.path.join(site_dir, "numpy", "f2py", "src", "fortranobject.c"),
                os.path.join(site_dir, "numpy", "f2py", "fortranobject.c"),
            ])
except:
    pass

# Method 4: sys.path based search
for path in sys.path:
    if path and "numpy" in path:
        candidates.extend([
            os.path.join(path, "numpy", "f2py", "src", "fortranobject.c"),
            os.path.join(path, "numpy", "f2py", "fortranobject.c"),
        ])

# Find first existing file
for candidate in candidates:
    if os.path.exists(candidate):
        print(candidate)
        break
else:
    print("")  # empty if not found
'''],
        check: true,
    )
    fortranobject_c = fortranobject_cmd.stdout().strip()

    # Generate f2py interface files
    spherepack_sources = []
    spherepack_source = custom_target(
        '_spherepackmodule.c',
        input: ['src/_spherepack.pyf'],
        output: ['_spherepackmodule.c', '_spherepack-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    spherepack_sources += spherepack_source

    # Add fortranobject.c if found
    if fortranobject_c != ''
        message('Found fortranobject.c: ' + fortranobject_c)
        fortran_obj_check = run_command(
            py,
            '-c', 'import os; import sys; sys.exit(0 if os.path.exists(r"' + fortranobject_c + '") else 1)',
            check: false,
        )
        if fortran_obj_check.returncode() == 0
            spherepack_sources += files(fortranobject_c)
            message('Added fortranobject.c to build')
        else
            warning('fortranobject.c path found but file does not exist: ' + fortranobject_c)
        endif
    else
        warning('fortranobject.c not found - trying to build without it (may cause issues)')
    endif

    # Build dependency list
    dependencies = [py_dep]
    if openmp_dep.found()
        dependencies += openmp_dep
        message('Using OpenMP for parallel processing')
    else
        warning('OpenMP not found - building without parallel support')
    endif

    # Build the spherepack extension module
    spherepack_ext = py.extension_module(
        '_spherepack',
        [
            # Core SPHEREPACK Fortran sources
            'src/lap.f90',
            'src/invlap.f90',
            'src/sphcom.f90',
            'src/hrfft.f90',
            'src/getlegfunc.f90',
            'src/specintrp.f90',
            'src/onedtotwod.f90',
            'src/onedtotwod_vrtdiv.f90',
            'src/twodtooned.f90',
            'src/twodtooned_vrtdiv.f90',
            'src/multsmoothfact.f90',
            'src/gaqd.f90',
            'src/shses.f90',
            'src/shaes.f90',
            'src/vhaes.f90',
            'src/vhses.f90',
            'src/shsgs.f90',
            'src/shags.f90',
            'src/vhags.f90',
            'src/vhsgs.f90',
            'src/shaec.f90',
            'src/shagc.f90',
            'src/shsec.f90',
            'src/shsgc.f90',
            'src/vhaec.f90',
            'src/vhagc.f90',
            'src/vhsec.f90',
            'src/vhsgc.f90',
            'src/ihgeod.f90',
            'src/alf.f90',
            spherepack_sources,
        ],
        include_directories: inc_np,
        dependencies: dependencies,
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,  # Install the extension module properly
        build_by_default: true,
    )

    # Automatically move compiled extension files to source directory (for --inplace builds)
    # This custom_target executes immediately after spherepack extension compilation completes
    move_extension = custom_target(
        'move_spherepack_extension',
        input: spherepack_ext,  # Depends on compiled extension
        output: 'moved_marker.txt',  # Create marker file to indicate move completion
        command: [
            py,
            '-c',
            '''
import os
import shutil
import glob
import time
import platform
from pathlib import Path

def move_extension_files():
    # Get build directory and source directory
    build_dir = Path("@BUILD_ROOT@")
    source_dir = Path("@SOURCE_ROOT@")

    print(f"Moving compiled extension files from {build_dir} to {source_dir}")

    # Search for compiled extension files
    patterns = ["_spherepack*.so", "_spherepack*.pyd", "_spherepack*.dylib"]
    compiled_files = []

    for pattern in patterns:
        search_pattern = str(build_dir / "**" / pattern)
        found_files = glob.glob(search_pattern, recursive=True)
        compiled_files.extend(found_files)

    if compiled_files:
        for source_file_str in compiled_files:
            source_file = Path(source_file_str)
            destination = source_dir / source_file.name

            print(f"Moving compiled extension: {source_file} -> {destination}")

            # Windows-specific handling for file locking issues
            if platform.system() == "Windows":
                success = False
                max_attempts = 5

                for attempt in range(max_attempts):
                    try:
                        # Handle existing target file first
                        if destination.exists():
                            try:
                                destination.unlink()
                            except OSError:
                                # If file is locked, rename it first
                                backup_name = destination.with_suffix(f".old{int(time.time())}")
                                destination.rename(backup_name)

                        # Copy file (safer than move on Windows)
                        shutil.copy2(str(source_file), str(destination))
                        print(f"Successfully copied {source_file.name}")
                        success = True
                        break

                    except (OSError, shutil.Error) as e:
                        if attempt < max_attempts - 1:
                            print(f"Attempt {attempt + 1} failed: {e}. Retrying in 1 second...")
                            time.sleep(1)
                        else:
                            print(f"Final attempt failed: {e}")

                if not success:
                    print(f"Warning: Unable to copy {source_file.name} after {max_attempts} attempts")
                    print(f"New file is available at: {source_file}")
                    print(f"Please manually copy to: {destination}")

            else:
                # Unix-like systems (Linux/macOS) - simple handling
                try:
                    # Remove existing file if it exists
                    if destination.exists():
                        destination.unlink()

                    # Move file
                    shutil.move(str(source_file), str(destination))
                    print(f"Successfully moved {source_file.name}")

                except (OSError, shutil.Error) as e:
                    print(f"Warning: Unable to move {source_file.name}: {e}")
                    # Try copying instead
                    try:
                        shutil.copy2(str(source_file), str(destination))
                        print(f"Successfully copied {source_file.name}")
                    except (OSError, shutil.Error) as copy_e:
                        print(f"Error: Unable to copy {source_file.name}: {copy_e}")
    else:
        print("Warning: No compiled extension found to move")

    # Create marker file
    with open("@OUTPUT@", "w") as f:
        f.write("Extension moved successfully\\n")

# Execute file moving
move_extension_files()
'''
        ],
        build_by_default: true,  # Execute by default during build
        console: true  # Display output
    )

    message('spharm submodule configured successfully with cross-platform optimizations')
else
    warning('f2py not available, skipping spharm interface generation')
    message('Installing Python-only spharm module (without Fortran extensions)')
endif
