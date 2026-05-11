# Skyborn GPI submodule - Genesis Potential Index and Tropical Cyclone Potential Intensity
project(
    'skyborn-gpi',
    ['c', 'fortran'],
    version: '1.0.0',
    default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()
build_py = find_program('build_python', required: false)
if not build_py.found()
    build_py = py
endif

py.install_sources(
    [
        '__init__.py',
        'core.py',
        'xarray.py',
    ],
    subdir: 'skyborn/calc/GPI',
)

incdir_numpy_cmd = run_command(
    build_py,
    [
        '-c',
        '''
import numpy
import os
paths = [
    numpy.get_include(),
    os.path.join(numpy.__path__[0], "core", "include"),
    os.path.join(numpy.__path__[0], "_core", "include"),
]
for path in paths:
    if os.path.exists(os.path.join(path, "numpy", "npy_os.h")):
        print(path)
        break
else:
    print(numpy.get_include())
''',
    ],
    check: true,
)
incdir_numpy = incdir_numpy_cmd.stdout().strip()

message('NumPy include directory: ' + incdir_numpy)

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
    endif
endif

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
    '-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION',
    '-fPIC',
    '-fno-trapping-math',
]

fortran_args += openmp_fortran_args
c_args += openmp_c_args

if host_system == 'darwin' and host_cpu == 'aarch64'
    fortran_args += ['-march=armv8-a', '-mtune=apple-m1']
    c_args += ['-march=armv8-a', '-mtune=apple-m1']
elif host_system == 'windows'
    fortran_args += ['-march=x86-64', '-mtune=generic']
    c_args += ['-march=x86-64', '-mtune=generic']
else
    fortran_args += ['-march=x86-64', '-mtune=generic']
    c_args += ['-march=x86-64', '-mtune=generic']
endif

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

gpi_ext = py.extension_module(
    'tropical_cyclone_potential_intensity',
    [
        'tropical_cyclone_potential_intensity_module.c',
        'fortran/tropical_cyclone_potential_intensity.f90',
        'fortran/tropical_cyclone_potential_intensity_bindings.f90',
    ],
    include_directories: include_directories(incdir_numpy),
    dependencies: [py_dep] + openmp_dependencies,
    fortran_args: fortran_args,
    c_args: c_args,
    link_args: openmp_link_args,
    install: true,
    install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'GPI',
    build_by_default: true,
)

custom_target(
    'copy_gpi_to_source',
    input: gpi_ext,
    output: 'inplace_copy_marker.txt',
    command: [
        py,
        '-c',
        '''import shutil
import sys
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent

for pattern in (
    'tropical_cyclone_potential_intensity*.so',
    'tropical_cyclone_potential_intensity*.pyd',
    'tropical_cyclone_potential_intensity*.dylib',
):
    for built_file in current_dir.rglob(pattern):
        target = module_dir / built_file.name
        if not target.exists() or built_file.stat().st_mtime > target.stat().st_mtime:
            shutil.copy2(str(built_file), str(target))
            print(f'Copied {built_file.name} to source')

with open(Path(sys.argv[1]), 'w', encoding='utf-8') as marker:
    marker.write('gpi copied\\n')
''',
        '@OUTPUT@',
    ],
    build_by_default: true,
    console: true,
)

message('GPI submodule configured successfully with direct Fortran/C extension')
