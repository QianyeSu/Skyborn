# Skyborn troposphere submodule - Tropopause calculation with Fortran extensions
project('skyborn-troposphere', 'c', 'fortran',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

# Get Python installation
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Find OpenMP dependency
openmp_dep = dependency('openmp', required: false)

# Install Python module files
py.install_sources([
    '__init__.py',
    'tropopause.py',
    'tropopause_xarray.py'
], subdir: 'skyborn/calc/troposphere')

# Get numpy include directories
incdir_numpy_cmd = run_command(
    py,
    ['-c', '''
import numpy
import os
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
    print(numpy.get_include())
'''],
    check: true,
)
incdir_numpy = incdir_numpy_cmd.stdout().strip()

# Try to get f2py include
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

# Detect platform for optimization flags
host_system = host_machine.system()
host_cpu = host_machine.cpu_family()

message('Host system: ' + host_system)
message('Host CPU family: ' + host_cpu)

# Optimization flags
fortran_args = [
    '-O3',
    '-fPIC',
    '-fno-second-underscore',
    '-funroll-loops',
    '-finline-functions',
    '-ftree-vectorize',
    '-ffree-line-length-none',
    '-fno-common',
    '-std=legacy',
    '-fopenmp',
]

c_args = [
    '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
    '-fPIC',
    '-fno-trapping-math',
    '-fopenmp',
]

if host_system == 'darwin' and host_cpu == 'aarch64'
    # Apple Silicon optimizations
    fortran_args += ['-march=armv8-a', '-mtune=apple-m1']
    c_args += ['-march=armv8-a', '-mtune=apple-m1']
else
    # x86-64 optimizations
    fortran_args += ['-march=x86-64', '-mtune=generic']
    c_args += ['-march=x86-64', '-mtune=generic']
endif

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

# Generate Fortran-Python interface with f2py
try_f2py = run_command(py, '-c', 'import numpy.f2py; print("OK")', check: false)

if try_f2py.returncode() == 0
    inc_np = include_directories(incdir_numpy, incdir_f2py)

    # Get fortranobject.c path
    fortranobject_cmd = run_command(
        py,
        ['-c', '''
import numpy.f2py
import numpy
import os
import sys

candidates = []

try:
    fo_path = os.path.join(numpy.f2py.get_include(), "fortranobject.c")
    candidates.append(fo_path)
except:
    pass

numpy_base = numpy.__path__[0]
candidates.extend([
    os.path.join(numpy_base, "f2py", "src", "fortranobject.c"),
    os.path.join(numpy_base, "f2py", "fortranobject.c"),
    os.path.join(numpy_base, "_core", "src", "fortranobject.c"),
])

for candidate in candidates:
    if os.path.exists(candidate):
        print(candidate)
        break
else:
    print("")
'''],
        check: true,
    )
    fortranobject_c = fortranobject_cmd.stdout().strip()

    # Generate f2py interface files for tropopause calculation
    tropopause_sources = []
    tropopause_source = custom_target(
        'tropopause_heightmodule.c',
        input: ['fortran/tropopause_height.pyf'],
        output: ['tropopause_heightmodule.c', 'tropopause_height-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    tropopause_sources += tropopause_source

    # Add fortranobject.c if found
    if fortranobject_c != ''
        message('Found fortranobject.c: ' + fortranobject_c)
        fortran_obj_check = run_command(
            py,
            '-c', 'import os; import sys; sys.exit(0 if os.path.exists(r"' + fortranobject_c + '") else 1)',
            check: false,
        )
        if fortran_obj_check.returncode() == 0
            tropopause_sources += files(fortranobject_c)
            message('Added fortranobject.c to build')
        else
            warning('fortranobject.c path found but file does not exist: ' + fortranobject_c)
        endif
    else
        warning('fortranobject.c not found - trying to build without it')
    endif

    # Build dependency list
    dependencies = [py_dep]
    if openmp_dep.found()
        dependencies += openmp_dep
        message('Using OpenMP for parallel processing')
    else
        warning('OpenMP not found - building without parallel support')
    endif

    # Build the tropopause extension module
    tropopause_ext = py.extension_module(
        'tropopause_height',
        [
            # Core tropopause Fortran sources
            'fortran/stattrop_dp.f',
            'fortran/stndAtmUS76_dp.f',
            'fortran/tropopause_height.f90',
            tropopause_sources,
        ],
        include_directories: inc_np,
        dependencies: dependencies,
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'troposphere',
        build_by_default: true,
    )

    # Copy to source directory for development
    copy_to_source = custom_target(
        'copy_tropopause_to_source',
        input: tropopause_ext,
        output: 'tropopause_copy_marker.txt',
        command: [
            py,
            '-c',
            '''
import os
import shutil
import glob
from pathlib import Path

current_dir = Path.cwd()
build_dir = current_dir / "build"

patterns = ["tropopause_height*.so", "tropopause_height*.pyd", "tropopause_height*.dylib"]
built_files = []

for pattern in patterns:
    found_files = glob.glob(str(build_dir / "**" / pattern), recursive=True)
    built_files.extend(found_files)

for built_file_str in built_files:
    built_file = Path(built_file_str)
    source_file = current_dir / built_file.name

    try:
        if source_file.exists():
            source_file.unlink()
        shutil.copy2(str(built_file), str(source_file))
        print(f"✓ Copied {built_file.name} to source directory")
    except Exception as e:
        print(f"✗ Failed to copy {built_file.name}: {e}")

with open("''' + '@OUTPUT@' + '''", "w") as f:
    f.write("Tropopause extension copy completed\\n")
'''
        ],
        build_by_default: true,
        console: true
    )

    message('troposphere submodule configured successfully with tropopause extensions')
else
    warning('f2py not available, skipping tropopause interface generation')
    message('Installing Python-only troposphere module (without tropopause Fortran extensions)')
endif
