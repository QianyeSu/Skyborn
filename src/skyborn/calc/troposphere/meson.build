# Skyborn troposphere submodule - Tropopause calculation with direct Fortran/C extension
project(
  'skyborn-troposphere',
  ['c', 'fortran'],
  version: '1.0.0',
  default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

py.install_sources(
  [
    '__init__.py',
    'core.py',
    'xarray.py',
  ],
  subdir: 'skyborn/calc/troposphere',
)

incdir_numpy_cmd = run_command(
  py,
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

openmp_dependencies = []
openmp_fortran_args = []
openmp_c_args = []
if host_system == 'windows'
  message('Skipping OpenMP flags on Windows for troposphere build')
else
  openmp_dep = dependency('openmp', required: false)
  if openmp_dep.found()
    openmp_dependencies += openmp_dep
    message('Using OpenMP dependency for parallel processing')
  else
    warning('OpenMP dependency not found; relying on compiler flags only')
  endif
  openmp_fortran_args += ['-fopenmp']
  openmp_c_args += ['-fopenmp']
endif

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

fortran_args += openmp_fortran_args
c_args += openmp_c_args

if host_system == 'darwin' and host_cpu == 'aarch64'
  fortran_args += ['-march=armv8-a', '-mtune=apple-m1']
  c_args += ['-march=armv8-a', '-mtune=apple-m1']
else
  fortran_args += ['-march=x86-64', '-mtune=generic']
  c_args += ['-march=x86-64', '-mtune=generic']
endif

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

tropopause_ext = py.extension_module(
  'tropopause_height',
  [
    'tropopause_height_module.c',
    'fortran/stattrop_dp.f',
    'fortran/stndAtmUS76_dp.f',
    'fortran/tropopause_height.f90',
    'fortran/tropopause_height_bindings.f90',
  ],
  include_directories: inc_np,
  dependencies: [py_dep] + openmp_dependencies,
  c_args: c_args,
  fortran_args: fortran_args,
  install: true,
  install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'troposphere',
  build_by_default: true,
)

copy_to_source = custom_target(
  'copy_to_source_for_inplace',
  input: tropopause_ext,
  output: 'inplace_copy_marker.txt',
  command: [
    py,
    '-c',
    '''
import glob
import shutil
from pathlib import Path

print("=== Extension source copy check ===")

current_dir = Path.cwd()
module_dir = current_dir.parent
print(f"Current dir (build): {current_dir}")
print(f"Module dir (source): {module_dir}")

patterns = ["tropopause_height*.so", "tropopause_height*.pyd", "tropopause_height*.dylib"]
built_files = []

for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / pattern), recursive=False))
for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / "**" / pattern), recursive=True))

built_files = list(set(built_files))

if built_files:
    print(f"Found built files: {built_files}")
    for built_file_str in built_files:
        built_file = Path(built_file_str)
        source_file = module_dir / built_file.name
        should_copy = False

        if not source_file.exists():
            should_copy = True
            print(f"Source {built_file.name} does not exist, will copy")
        elif built_file.stat().st_mtime > source_file.stat().st_mtime:
            should_copy = True
            print(f"Source {built_file.name} is older, will copy")

        if should_copy:
            try:
                if source_file.exists():
                    source_file.unlink()
                shutil.copy2(str(built_file), str(source_file))
                print(f"Copied {built_file.name} to source directory successfully")
            except Exception as exc:
                print(f"Failed to copy {built_file.name}: {exc}")
else:
    print("No compiled extension files found")
    print("Available files in current directory:")
    if current_dir.exists():
        for f in current_dir.rglob("*"):
            if f.is_file():
                print(f"  {f}")
    else:
        print("  Current directory does not exist")

with open("''' + '@OUTPUT@' + '''", "w", encoding="utf-8") as handle:
    handle.write("Source copy check completed\\n")

print("=== Extension source copy check completed ===")
''',
  ],
  build_by_default: true,
  console: true,
)

message('troposphere submodule configured successfully with direct Fortran/C extension')
