# Skyborn interp submodule - Curvilinear grid interpolation with Fortran acceleration
project(
    'skyborn-interp',
    'c',
    'fortran',
    version: '1.0.0',
    default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

# Get Python installation
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Install Python module files
py.install_sources(
    [
        '__init__.py',
        'interpolation.py',
        'regridding.py',
        'rcm2points.py',
        'rcm2rgrid.py',
        'triple_to_grid.py',
        'errors.py',
        'missing_values.py',
    ],
    subdir: 'skyborn/interp',
)

# Get numpy include directories with robust fallbacks
incdir_numpy_cmd = run_command(
    py,
    [
        '-c', '''
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
''',
    ],
    check: true,
)
incdir_numpy = incdir_numpy_cmd.stdout().strip()

# Try to get f2py include, with multiple fallbacks
incdir_f2py_cmd = run_command(
    py,
    [
        '-c', '''
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
''',
    ],
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

# Cross-platform optimization flags
fortran_args = []
c_args = []

if host_system == 'darwin' and host_cpu == 'aarch64'
    # Apple Silicon (arm64) optimized flags
    message('Configuring for Apple Silicon (arm64)')
    fortran_args = [
        '-O3',
        '-fPIC',
        '-fno-second-underscore',
        '-funroll-loops',
        '-finline-functions',
        '-ftree-vectorize',
        '-ffinite-math-only',
        '-ffree-line-length-none',
        '-march=armv8-a',
        '-mtune=apple-m1',
        '-fno-common',
        '-ftree-loop-im',
        '-ftree-loop-distribution',
        '-falign-functions=32',
        '-falign-loops=32',
        '-std=legacy',
    ]
    c_args = [
        '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
        '-march=armv8-a',
        '-mtune=apple-m1',
        '-fPIC',
        '-fno-trapping-math',
        '-falign-functions=32',
    ]
else
    # x86-64 optimized flags (Linux/Windows/Intel macOS)
    message('Configuring for x86-64 architecture')
    fortran_args = [
        '-O3',
        '-fPIC',
        '-fno-second-underscore',
        '-funroll-loops',
        '-finline-functions',
        '-ftree-vectorize',
        '-ffinite-math-only',
        '-ffree-line-length-none',
        '-march=x86-64',
        '-mtune=generic',
        '-fno-common',
        '-ftree-loop-im',
        '-ftree-loop-distribution',
        '-falign-functions=32',
        '-falign-loops=32',
        '-std=legacy',
    ]
    c_args = [
        '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
        '-march=x86-64',
        '-mtune=generic',
        '-fPIC',
        '-fno-trapping-math',
        '-falign-functions=32',
    ]
endif

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

# Generate Fortran-Python interface with f2py
try_f2py = run_command(py, '-c', 'import numpy.f2py; print("OK")', check: false)

if try_f2py.returncode() == 0
    inc_np = include_directories(incdir_numpy, incdir_f2py)

    # Get fortranobject.c path with extensive fallbacks
    fortranobject_cmd = run_command(
        py,
        [
            '-c', '''
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

# Find first existing file
for candidate in candidates:
    if os.path.exists(candidate):
        print(candidate)
        break
else:
    print("")  # empty if not found
''',
        ],
        check: true,
    )
    fortranobject_c = fortranobject_cmd.stdout().strip()

    # Build RCM2POINTS extension
    rcm2points_sources = []
    rcm2points_source = custom_target(
        'rcm2pointsmodule.c',
        input: ['fortran/rcm2points.pyf'],
        output: ['rcm2pointsmodule.c', 'rcm2points-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    rcm2points_sources += rcm2points_source

    # Build RCM2RGRID extension
    rcm2rgrid_sources = []
    rcm2rgrid_source = custom_target(
        'rcm2rgridmodule.c',
        input: ['fortran/rcm2rgrid.pyf'],
        output: ['rcm2rgridmodule.c', 'rcm2rgrid-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    rcm2rgrid_sources += rcm2rgrid_source

    # Build TRIPLE2GRID extension
    triple2grid_sources = []
    triple2grid_source = custom_target(
        'triple2gridmodule.c',
        input: ['fortran/triple2grid.pyf'],
        output: ['triple2gridmodule.c', 'triple2grid-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    triple2grid_sources += triple2grid_source

    # Build GRID2TRIPLE extension
    grid2triple_sources = []
    grid2triple_source = custom_target(
        'grid2triplemodule.c',
        input: ['fortran/grid2triple.pyf'],
        output: ['grid2triplemodule.c', 'grid2triple-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    grid2triple_sources += grid2triple_source

    # Add fortranobject.c if found
    if fortranobject_c != ''
        message('Found fortranobject.c: ' + fortranobject_c)
        fortran_obj_check = run_command(
            py,
            '-c', 'import os; import sys; sys.exit(0 if os.path.exists(r"'
            + fortranobject_c
            + '") else 1)',
            check: false,
        )
        if fortran_obj_check.returncode() == 0
            rcm2points_sources += files(fortranobject_c)
            rcm2rgrid_sources += files(fortranobject_c)
            triple2grid_sources += files(fortranobject_c)
            grid2triple_sources += files(fortranobject_c)
            message('Added fortranobject.c to build')
        else
            warning(
                'fortranobject.c path found but file does not exist: '
                + fortranobject_c,
            )
        endif
    else
        warning(
            'fortranobject.c not found - trying to build without it (may cause issues)',
        )
    endif

    # Build the rcm2points extension module
    rcm2points_ext = py.extension_module(
        'rcm2points',
        [
            'fortran/rcm2points.f',
            'fortran/rcm2rgrid.f',
            'fortran/linmsg_dp.f',
            'fortran/linint2.f',
            rcm2points_sources,
        ],
        include_directories: inc_np,
        dependencies: [py_dep],
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'interp' / 'fortran',
        build_by_default: true,
    )

    # Build the rcm2rgrid extension module
    rcm2rgrid_ext = py.extension_module(
        'rcm2rgrid',
        [
            'fortran/rcm2rgrid.f',
            'fortran/linmsg_dp.f',
            'fortran/linint2.f',
            rcm2rgrid_sources,
        ],
        include_directories: inc_np,
        dependencies: [py_dep],
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'interp' / 'fortran',
        build_by_default: true,
    )

    # Build the triple2grid extension module
    triple2grid_ext = py.extension_module(
        'triple2grid',
        [
            'fortran/triple2grid.f',
            triple2grid_sources,
        ],
        include_directories: inc_np,
        dependencies: [py_dep],
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'interp' / 'fortran',
        build_by_default: true,
    )

    # Build the grid2triple extension module
    grid2triple_ext = py.extension_module(
        'grid2triple',
        [
            'fortran/grid2triple.f',
            grid2triple_sources,
        ],
        include_directories: inc_np,
        dependencies: [py_dep],
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'interp' / 'fortran',
        build_by_default: true,
    )

    # Copy to source directory for inplace builds
    custom_target(
        'copy_interp_to_source',
        input: [rcm2points_ext, rcm2rgrid_ext, triple2grid_ext, grid2triple_ext],
        output: 'interp_copy_marker.txt',
        command: [
            py,
            '-c',
            '''import os
import shutil
import glob
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent
fortran_dir = module_dir / 'fortran'
fortran_dir.mkdir(exist_ok=True)

patterns = ['*.so', '*.pyd', '*.dylib']
for pattern in patterns:
    for built_file in current_dir.rglob(pattern):
        target = fortran_dir / built_file.name
        if not target.exists() or built_file.stat().st_mtime > target.stat().st_mtime:
            shutil.copy2(str(built_file), str(target))
            print(f'Copied {built_file.name} to source')

import sys
with open(sys.argv[1], 'w') as f:
    f.write('Interp extensions copied' + chr(10))
''',
            '@OUTPUT@',
        ],
        build_by_default: true,
        console: true,
    )

    # Create __init__.py for fortran submodule
    init_file = custom_target(
        'fortran_init',
        output: '__init__.py',
        command: [
            py,
            '-c',
            '''import sys
with open(sys.argv[1], 'w') as f:
    f.write('# Fortran interpolation extensions' + chr(10))
    # Import Fortran-exported symbols from their extension modules
    f.write('from .rcm2points import drcm2points' + chr(10))
    f.write('from .rcm2rgrid import drcm2rgrid, drgrid2rcm' + chr(10))
    f.write('from .triple2grid import triple2grid1' + chr(10))
    f.write('from .grid2triple import grid2triple' + chr(10))
    f.write(chr(10))
    f.write('__all__ = [' + chr(39) + 'drcm2points' + chr(39) + ', ' + chr(39) + 'drcm2rgrid' + chr(39) + ', ' + chr(39) + 'drgrid2rcm' + chr(39) + ', ' + chr(39) + 'triple2grid1' + chr(39) + ', ' + chr(39) + 'grid2triple' + chr(39) + ']' + chr(10))
''',
            '@OUTPUT@',
        ],
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'interp' / 'fortran',
    )

    message('interp submodule configured successfully with Fortran acceleration')
else
    warning('f2py not available, skipping interp Fortran extensions')
endif
