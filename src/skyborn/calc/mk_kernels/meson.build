# Skyborn private Mann-Kendall kernels - optional Fortran backend
project('skyborn-mk-kernels', 'c', 'fortran',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

py.install_sources([
    '__init__.py',
], subdir: 'skyborn/calc/mk_kernels')

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
    '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
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

try_f2py = run_command(py, '-c', 'import numpy.f2py; print("OK")', check: false)

if try_f2py.returncode() == 0
    inc_np = include_directories(incdir_numpy, incdir_f2py)

    fortranobject_cmd = run_command(
        py,
        ['-c', '''
import numpy.f2py
import numpy
import os

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

    mk_sources = []
    mk_source = custom_target(
        'mk_kernels_backendmodule.c',
        input: ['fortran/mk_kernels_backend.pyf'],
        output: ['mk_kernels_backendmodule.c', 'mk_kernels_backend-f2pywrappers.f'],
        command: [
            py,
            '-c',
            '''
import builtins
import sys
from numpy.f2py import f2py2e
from numpy.f2py import crackfortran

def _patched_openhook(filename, mode):
    return builtins.open(filename, mode, encoding="utf-8")

crackfortran.openhook = _patched_openhook
f2py2e.run_main(["--lower", sys.argv[1], "--build-dir", sys.argv[2]])
''',
            '@INPUT@',
            '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    mk_sources += mk_source

    if fortranobject_c != ''
        message('Found fortranobject.c: ' + fortranobject_c)
        fortran_obj_check = run_command(
            py,
            '-c', 'import os; import sys; sys.exit(0 if os.path.exists(r"' + fortranobject_c + '") else 1)',
            check: false,
        )
        if fortran_obj_check.returncode() == 0
            mk_sources += files(fortranobject_c)
            message('Added fortranobject.c to build')
        else
            warning('fortranobject.c path found but file does not exist: ' + fortranobject_c)
        endif
    else
        warning('fortranobject.c not found - trying to build without it')
    endif

    mk_ext = py.extension_module(
        'mk_kernels_backend',
        [
            'fortran/mk_kernels_backend.f90',
            mk_sources,
        ],
        include_directories: inc_np,
        dependencies: [py_dep],
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'mk_kernels',
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

patterns = ["mk_kernels_backend*.so", "mk_kernels_backend*.pyd", "mk_kernels_backend*.dylib"]
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

with open("''' + '@OUTPUT@' + '''", "w") as f:
    f.write("Source copy check completed\\n")

print("=== Extension source copy check completed ===")
'''
        ],
        build_by_default: true,
        console: true,
    )

    message('mk_kernels submodule configured successfully with Mann-Kendall extensions')
else
    warning('f2py not available, skipping mk_kernels interface generation')
    message('Installing Python-only mk_kernels module (without Fortran extensions)')
endif
