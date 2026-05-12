# Skyborn DCAPE submodule - direct Fortran/C extension
project(
  'skyborn-dcape',
  ['c', 'fortran'],
  version: '1.0.0',
  default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()
build_py = find_program('build_python', required: false)
if not build_py.found()
  build_py = py
endif

py.install_sources(
  [
    '__init__.py',
    'core.py',
  ],
  subdir: 'skyborn/calc/dcape',
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
inc_np = include_directories(incdir_numpy_cmd.stdout().strip())

host_system = host_machine.system()
host_cpu = host_machine.cpu_family()
cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')

message('NumPy include directory: ' + incdir_numpy_cmd.stdout().strip())
message('Host system: ' + host_system)
message('Host CPU family: ' + host_cpu)

c_args = ['-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION']
fortran_args = [
  '-O3',
  '-fPIC',
  '-fno-second-underscore',
  '-funroll-loops',
  '-finline-functions',
  '-ftree-vectorize',
  '-ffree-line-length-none',
  '-fno-common',
]

if host_system == 'windows' and cc.get_id() == 'msvc'
  c_args += ['/O2']
else
  c_args += ['-O3', '-fPIC', '-fno-trapping-math']
endif

if fc.get_id() == 'gcc'
  fortran_args += ['-fimplicit-none']
endif

if host_system == 'darwin' and host_cpu == 'aarch64'
  fortran_args += ['-march=armv8-a', '-mtune=apple-m1']
  c_args += ['-march=armv8-a', '-mtune=apple-m1']
else
  fortran_args += ['-march=x86-64', '-mtune=generic']
  c_args += ['-march=x86-64', '-mtune=generic']
endif

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

dcape_ext = py.extension_module(
  'dcape_core',
  [
    'dcape_core_module.c',
    'fortran/dcape.f90',
    'fortran/dcape_bindings.f90',
  ],
  include_directories: inc_np,
  dependencies: [py_dep],
  c_args: c_args,
  fortran_args: fortran_args,
  install: true,
  install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'dcape',
  build_by_default: true,
)

custom_target(
  'copy_to_source_for_inplace',
  input: dcape_ext,
  output: 'inplace_copy_marker.txt',
  command: [
    py,
    '-c',
    '''
import glob
import shutil
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent

patterns = ["dcape_core*.so", "dcape_core*.pyd", "dcape_core*.dylib"]
built_files = []
for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / pattern), recursive=False))
for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / "**" / pattern), recursive=True))

for built_file_str in sorted(set(built_files)):
    built_file = Path(built_file_str)
    target = module_dir / built_file.name
    if not target.exists() or built_file.stat().st_mtime > target.stat().st_mtime:
        shutil.copy2(str(built_file), str(target))
        print(f"Copied {built_file.name} to source")

with open(Path("@OUTPUT@"), "w", encoding="utf-8") as marker:
    marker.write("dcape copied\\n")
''',
  ],
  build_by_default: true,
  console: true,
)

message('dcape submodule configured successfully with direct Fortran/C extension')
