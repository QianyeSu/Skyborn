# Skyborn private Mann-Kendall kernels
project(
    'skyborn-mann-kendall-core',
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
    ],
    subdir: 'skyborn/calc/mann_kendall',
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

c_args = [
    '-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION',
    '-fPIC',
    '-fno-trapping-math',
]

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

mk_ext = py.extension_module(
    'mann_kendall_core',
    [
        'mann_kendall_core_module.c',
        'fortran/mann_kendall_core.f90',
        'fortran/mann_kendall_core_bindings.f90',
    ],
    include_directories: include_directories(incdir_numpy),
    dependencies: [py_dep],
    fortran_args: fortran_args,
    c_args: c_args,
    install: true,
    install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'mann_kendall',
    build_by_default: true,
)

copy_to_source = custom_target(
    'copy_to_source_for_inplace',
    input: mk_ext,
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

patterns = ["mann_kendall_core*.so", "mann_kendall_core*.pyd", "mann_kendall_core*.dylib"]
built_files = []

for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / pattern), recursive=False))
    built_files.extend(glob.glob(str(current_dir / "**" / pattern), recursive=True))

built_files = list(set(built_files))

if built_files:
    print(f"Found built files: {built_files}")
    for built_file_str in built_files:
        built_file = Path(built_file_str)
        source_file = module_dir / built_file.name

        if not source_file.exists():
            should_copy = True
            print(f"Source {built_file.name} doesn't exist, will copy")
        elif built_file.stat().st_mtime > source_file.stat().st_mtime:
            should_copy = True
            print(f"Source {built_file.name} is older, will copy")
        else:
            should_copy = False

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

with open("''' + '@OUTPUT@' + '''", "w", encoding="utf-8") as f:
    f.write("Source copy check completed\\n")

print("=== Extension source copy check completed ===")
''',
    ],
    build_by_default: true,
    console: true,
)

message('mann_kendall_core submodule configured successfully with direct Fortran/C extension')
