# Skyborn GPI submodule - Genesis Potential Index and Tropical Cyclone Potential Intensity with Fortran extensions
project('skyborn-gpi', 'c', 'fortran',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

# Get Python installation
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Install Python module files
py.install_sources([
    '__init__.py',
    'core.py',
    'xarray.py',
], subdir: 'skyborn/calc/GPI')

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

openmp_dependencies = []
openmp_fortran_args = []
openmp_c_args = []
openmp_link_args = []
fortran_compiler = meson.get_compiler('fortran')
c_compiler = meson.get_compiler('c')
is_gcc_toolchain = fortran_compiler.get_id() == 'gcc' and c_compiler.get_id() == 'gcc'
if is_gcc_toolchain
    if host_system == 'windows'
        message('Enabling OpenMP flags for GCC-based Windows GPI build')
    endif
    openmp_dep = dependency('openmp', required: false)
    if openmp_dep.found()
        openmp_dependencies += openmp_dep
        message('Using OpenMP dependency for parallel processing')
    else
        warning('OpenMP dependency not found - relying on compiler/linker flags only')
    endif
    openmp_fortran_args += ['-fopenmp']
    openmp_c_args += ['-fopenmp']
    openmp_link_args += ['-fopenmp']
elif host_system != 'windows'
    openmp_dep = dependency('openmp', required: false)
    if openmp_dep.found()
        openmp_dependencies += openmp_dep
        message('Using OpenMP dependency for parallel processing')
    else
        warning('OpenMP dependency not found')
    endif
endif

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
]

c_args = [
    '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
    '-fPIC',
    '-fno-trapping-math',
]

fortran_args += openmp_fortran_args
c_args += openmp_c_args

if host_system == 'darwin' and host_cpu == 'aarch64'
    # Apple Silicon optimizations
    fortran_args += ['-march=armv8-a', '-mtune=apple-m1']
    c_args += ['-march=armv8-a', '-mtune=apple-m1']
elif host_system == 'windows'
    # Windows-specific optimizations
    fortran_args += ['-march=x86-64', '-mtune=generic']
    c_args += ['-march=x86-64', '-mtune=generic']
else
    # Linux x86-64 optimizations
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

    # Generate f2py interface files for tropical cyclone potential intensity
    gpi_sources = []
    gpi_source = custom_target(
        'tropical_cyclone_potential_intensitymodule.c',
        input: ['fortran/tropical_cyclone_potential_intensity.pyf'],
        output: ['tropical_cyclone_potential_intensitymodule.c', 'tropical_cyclone_potential_intensity-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    gpi_sources += gpi_source

    # Add fortranobject.c if found
    if fortranobject_c != ''
        message('Found fortranobject.c: ' + fortranobject_c)
        fortran_obj_check = run_command(
            py,
            '-c', 'import os; import sys; sys.exit(0 if os.path.exists(r"' + fortranobject_c + '") else 1)',
            check: false,
        )
        if fortran_obj_check.returncode() == 0
            gpi_sources += files(fortranobject_c)
            message('Added fortranobject.c to build')
        else
            warning('fortranobject.c path found but file does not exist: ' + fortranobject_c)
        endif
    else
        warning('fortranobject.c not found - trying to build without it')
    endif

    # Build dependency list
    dependencies = [py_dep] + openmp_dependencies

    # Build the GPI extension module
    # Use subdir parameter to install directly into the Python package structure
    gpi_ext = py.extension_module(
        'tropical_cyclone_potential_intensity',
        [
            # Core GPI Fortran sources
            'fortran/tropical_cyclone_potential_intensity.f90',
            gpi_sources,
        ],
        include_directories: inc_np,
        dependencies: dependencies,
        fortran_args: fortran_args,
        c_args: c_args,
        link_args: openmp_link_args,
        install: true,   # Install to custom location if provided
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'GPI',  # Use setuptools path if available
        build_by_default: true,
    )

    # Always create a copy target for source directory (for --inplace builds)
    # This target will run after the main extension is built, but only copy if needed
    copy_to_source = custom_target(
        'copy_to_source_for_inplace',
        input: gpi_ext,
        output: 'inplace_copy_marker.txt',
        command: [
            py,
            '-c',
            '''
import os
import shutil
import sys
from pathlib import Path

print("=== Extension source copy check ===")

current_dir = Path.cwd()  # This is the build directory itself
module_dir = current_dir.parent  # This is the module source directory
print(f"Current dir (build): {current_dir}")
print(f"Module dir (source): {module_dir}")

built_file = Path(sys.argv[1]).resolve()
source_file = module_dir / built_file.name

print(f"Built extension: {built_file}")

if not built_file.exists():
    raise FileNotFoundError(f"Expected built extension does not exist: {built_file}")

should_copy = False
if not source_file.exists():
    should_copy = True
    print(f"Source {source_file.name} doesn't exist, will copy")
elif built_file.stat().st_mtime > source_file.stat().st_mtime:
    should_copy = True
    print(f"Source {source_file.name} is older, will copy")
else:
    print(f"Source {source_file.name} already up to date")

if should_copy:
    if source_file.exists():
        try:
            source_file.unlink()
        except PermissionError:
            backup = source_file.with_suffix(source_file.suffix + '.locked_old')
            if backup.exists():
                backup.unlink()
            source_file.rename(backup)
    shutil.copy2(str(built_file), str(source_file))
    print(f"Copied {built_file.name} to source directory successfully")

# Create marker
with open("''' + '@OUTPUT@' + '''", "w") as f:
    f.write("Source copy check completed\\n")

print("=== Extension source copy check completed ===")
'''
            ,
            '@INPUT@',
        ],
        build_by_default: true,  # Always run this check
        console: true
    )

    message('GPI submodule configured successfully with tropical cyclone potential intensity extensions')
else
    warning('f2py not available, skipping GPI interface generation')
    message('Installing Python-only GPI module (without Fortran extensions)')
endif
