# Skyborn growth-rate compiled normal-mode kernels
project(
    'skyborn-growth-rate',
    'c',
    'fortran',
    version: '1.0.0',
    default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

py.install_sources(
    ['__init__.py'],
    subdir: 'skyborn/calc/growth_rate',
)

incdir_numpy_cmd = run_command(
    py,
    [
        '-c', '''
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
except Exception:
    print(numpy.get_include())
''',
    ],
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
    '-ffinite-math-only',
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
        [
            '-c', '''
import numpy.f2py
import numpy
import os

candidates = []
try:
    candidates.append(os.path.join(numpy.f2py.get_include(), "fortranobject.c"))
except Exception:
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
''',
        ],
        check: true,
    )
    fortranobject_c = fortranobject_cmd.stdout().strip()

    growth_rate_sources = []
    growth_rate_source = custom_target(
        'growth_rate_kernelsmodule.c',
        input: ['fortran/growth_rate_kernels.pyf'],
        output: ['growth_rate_kernelsmodule.c', 'growth_rate_kernels-f2pywrappers.f'],
        command: [
            py,
            '-m', 'numpy.f2py',
            '--lower', '@INPUT@',
            '--build-dir', '@OUTDIR@',
        ],
        build_by_default: true,
        install: false,
    )
    growth_rate_sources += growth_rate_source

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
            growth_rate_sources += files(fortranobject_c)
            message('Added fortranobject.c to build')
        else
            warning(
                'fortranobject.c path found but file does not exist: '
                + fortranobject_c,
            )
        endif
    else
        warning('fortranobject.c not found - trying to build without it')
    endif

    growth_rate_ext = py.extension_module(
        'growth_rate_kernels',
        [
            'fortran/growth_rate_kernels.f90',
            'fortran/qzhes.f90',
            'fortran/qzit.f90',
            'fortran/qzval.f90',
            'fortran/epslon.f90',
            growth_rate_sources,
        ],
        include_directories: inc_np,
        dependencies: [py_dep],
        fortran_args: fortran_args,
        c_args: c_args,
        install: true,
        install_dir: get_option('python.purelibdir') / 'skyborn' / 'calc' / 'growth_rate',
        build_by_default: true,
    )

    custom_target(
        'copy_growth_rate_to_source',
        input: growth_rate_ext,
        output: 'inplace_copy_marker.txt',
        command: [
            py,
            '-c',
            '''import shutil
import sys
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent
target_dir = module_dir

for pattern in ('growth_rate_kernels*.so', 'growth_rate_kernels*.pyd', 'growth_rate_kernels*.dylib'):
    for built_file in current_dir.rglob(pattern):
        target = target_dir / built_file.name
        if not target.exists() or built_file.stat().st_mtime > target.stat().st_mtime:
            shutil.copy2(str(built_file), str(target))
            print(f'Copied {built_file.name} to source')

with open(Path(sys.argv[1]), 'w', encoding='utf-8') as marker:
    marker.write('growth_rate copied\\n')
''',
            '@OUTPUT@',
        ],
        build_by_default: true,
        console: true,
    )

    message('growth_rate submodule configured successfully')
else
    warning('f2py not available, skipping growth-rate core interface generation')
endif
