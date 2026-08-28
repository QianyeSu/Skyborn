# Skyborn calc.filter submodule - Lanczos filtering with Fortran acceleration
project('skyborn-filter', 'c', 'fortran',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()
build_py = find_program('build_python', required: false)
if not build_py.found()
    build_py = py
endif

# Install Python sources
py.install_sources(
    ['__init__.py', 'lanczos.py'],
    subdir: 'skyborn/calc/filter',
)

# Get NumPy include directory
incdir_numpy_cmd = run_command(
    build_py,
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

# Get Python include directories
python_include_dirs_cmd = run_command(
    py,
    ['-c', '''
import os
import sys
import sysconfig
from pathlib import Path

paths = []

def add(path):
    if not path:
        return
    normalized = os.path.abspath(path)
    if os.path.exists(os.path.join(normalized, "Python.h")) and normalized not in paths:
        paths.append(normalized)

version = f"python{sys.version_info[0]}.{sys.version_info[1]}"

for key in ("INCLUDEPY", "CONFINCLUDEPY"):
    add(sysconfig.get_config_var(key))

for key in ("include", "platinclude"):
    add(sysconfig.get_path(key))

roots = [
    sys.prefix,
    sys.exec_prefix,
    sys.base_prefix,
    sys.base_exec_prefix,
    str(Path(sys.executable).resolve().parent.parent),
]

for root in roots:
    add(os.path.join(root, "include", version))
    add(os.path.join(root, "include"))
    add(os.path.join(root, "Include"))
    add(os.path.join(root, "Headers"))

for path in paths:
    print(path)
'''],
    check: true,
)
python_include_dirs_text = python_include_dirs_cmd.stdout().strip()
python_include_dirs = []
if python_include_dirs_text != ''
    python_include_dirs = python_include_dirs_text.split('\n')
endif

message('NumPy include directory: ' + incdir_numpy)
message('Python include directories: ' + python_include_dirs_text)

host_system = host_machine.system()
host_cpu = host_machine.cpu_family()

message('Host system: ' + host_system)
message('Host CPU family: ' + host_cpu)

package_install_dir = get_option('python.purelibdir') / 'skyborn' / 'calc' / 'filter'

# Configure OpenMP
openmp_dependencies = []
openmp_fortran_args = []
openmp_c_args = []
openmp_link_args = []
fortran_compiler = meson.get_compiler('fortran')
c_compiler = meson.get_compiler('c')
is_gcc_toolchain = fortran_compiler.get_id() == 'gcc' and c_compiler.get_id() == 'gcc'

if is_gcc_toolchain
    openmp_dep = dependency('openmp', required: false)
    if openmp_dep.found()
        openmp_dependencies += openmp_dep
        message('Using OpenMP dependency for parallel processing')
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
        warning('OpenMP not found - building without parallel support')
    endif
else
    message('Skipping OpenMP flags on Windows for filter build')
endif

# Configure compiler flags
fortran_args = []
c_args = []

if host_system == 'darwin' and host_cpu == 'aarch64'
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
    ]
    c_args = [
        '-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION',
        '-march=armv8-a',
        '-mtune=apple-m1',
        '-fPIC',
        '-fno-trapping-math',
        '-falign-functions=32',
    ]
else
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
    ]
    c_args = [
        '-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION',
        '-march=x86-64',
        '-mtune=generic',
        '-fPIC',
        '-fno-trapping-math',
        '-falign-functions=32',
    ]
endif

fortran_args += openmp_fortran_args
c_args += openmp_c_args

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

if python_include_dirs.length() == 1
    inc_np = include_directories(
        incdir_numpy,
        python_include_dirs[0],
    )
else
    inc_np = include_directories(incdir_numpy)
endif

dependencies = [py_dep] + openmp_dependencies

# Build extension module
lanczos_ext = py.extension_module(
    '_lanczos_core',
    ['_lanczos_core_module.c', '_lanczos_core.f90'],
    include_directories: inc_np,
    dependencies: dependencies,
    fortran_args: fortran_args,
    c_args: c_args,
    link_args: openmp_link_args,
    install: true,
    install_dir: package_install_dir,
    build_by_default: true,
)

# Copy to source directory for development
copy_to_source = custom_target(
    'copy_lanczos_to_source',
    input: lanczos_ext,
    output: 'lanczos_copy_marker.txt',
    command: [
        py,
        '-c',
        '''
import os
import shutil
import glob
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent
print(f"Current dir: {current_dir}")
print(f"Module dir: {module_dir}")

patterns = [
    "_lanczos_core*.so",
    "_lanczos_core*.pyd",
    "_lanczos_core*.dylib",
]

built_files = []
for pattern in patterns:
    found = glob.glob(str(current_dir / pattern), recursive=False)
    built_files.extend(found)
    found = glob.glob(str(current_dir / "**" / pattern), recursive=True)
    built_files.extend(found)

built_files = list(set(built_files))

if built_files:
    print(f"Found built files: {built_files}")
    for built_file_str in built_files:
        built_file = Path(built_file_str)
        source_file = module_dir / built_file.name

        should_copy = False
        if not source_file.exists():
            should_copy = True
        elif built_file.stat().st_mtime > source_file.stat().st_mtime:
            should_copy = True

        if should_copy:
            try:
                if source_file.exists():
                    source_file.unlink()
                shutil.copy2(str(built_file), str(source_file))
                print(f"Copied {built_file.name} to source directory")
            except Exception as e:
                print(f"Failed to copy {built_file.name}: {e}")
else:
    print("No compiled extension files found")

with open("''' + '@OUTPUT@' + '''", "w") as f:
    f.write("Copy completed\\n")
'''
    ],
    build_by_default: true,
    console: true
)

message('filter submodule configured successfully')
