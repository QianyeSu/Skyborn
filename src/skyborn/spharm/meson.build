# Skyborn spharm submodule - Spherical harmonic transforms with cross-platform optimization
project('skyborn-spharm', 'c', 'fortran',
  version : '1.0.0',
  default_options : ['warning_level=2', 'buildtype=release', 'b_lto=true'])

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()
build_py = find_program('build_python', required: false)
if not build_py.found()
    build_py = py
endif

py.install_sources(
    ['__init__.py', 'spherical_harmonics.py', 'reduced_gaussian.py'],
    subdir: 'skyborn/spharm',
)

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

package_install_dir = get_option('python.purelibdir') / 'skyborn' / 'spharm'

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
    message('Skipping OpenMP flags on Windows for spharm build')
endif

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
        '-std=legacy',
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
        '-std=legacy',
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

numpy_header_check = run_command(
    build_py,
    ['-c', 'import os; import numpy; header_path = os.path.join(numpy.get_include(), "numpy", "npy_os.h"); print("OK" if os.path.exists(header_path) else "MISSING", header_path)'],
    check: false,
)
message('NumPy header check: ' + numpy_header_check.stdout().strip())

dependencies = [py_dep] + openmp_dependencies

common_fortran_sources = [
    'src/lap.f90',
    'src/invlap.f90',
    'src/sphcom.f90',
    'src/hrfft.f90',
    'src/multsmoothfact.f90',
    'src/gaqd.f90',
    'src/shags.f90',
]

common_fortran_lib = static_library(
    '_spherepack_common',
    common_fortran_sources,
    dependencies: dependencies,
    fortran_args: fortran_args,
    link_args: openmp_link_args,
    build_by_default: true,
)

direct_fortran_sources = [
    'src/spherepack_bindings.f90',
    'src/onedtotwod.f90',
    'src/twodtooned.f90',
    'src/shaes.f90',
    'src/vhaes.f90',
    'src/vhags.f90',
    'src/shses.f90',
    'src/vhses.f90',
    'src/shsgs.f90',
    'src/vhsgs.f90',
    'src/shaec.f90',
    'src/vhaec.f90',
    'src/shagc.f90',
    'src/vhagc.f90',
    'src/shsec.f90',
    'src/vhsec.f90',
    'src/shsgc.f90',
    'src/vhsgc.f90',
    'src/ihgeod.f90',
    'src/getlegfunc.f90',
    'src/specintrp.f90',
    'src/alf.f90',
]

reduced_direct_fortran_sources = [
    'src/reduced_bindings.f90',
    'src/reduced_gaussian_scalar.f90',
]

spherepack_direct_ext = py.extension_module(
    '_spherepack',
    ['_spherepack_module.c', '_spherepack_reduced_module.c'] + direct_fortran_sources + reduced_direct_fortran_sources,
    include_directories: inc_np,
    dependencies: dependencies,
    link_with: common_fortran_lib,
    fortran_args: fortran_args,
    c_args: c_args,
    link_args: openmp_link_args,
    install: true,
    install_dir: package_install_dir,
    build_by_default: true,
)

copy_targets = [spherepack_direct_ext]

copy_to_source = custom_target(
    'copy_to_source_for_inplace',
    input: copy_targets,
    output: 'inplace_copy_marker.txt',
    command: [
        py,
        '-c',
        '''
import os
import shutil
import glob
import sys
from pathlib import Path

print("=== Extension source copy check ===")

current_dir = Path.cwd()
module_dir = current_dir.parent
print(f"Current dir (build): {current_dir}")
print(f"Module dir (source): {module_dir}")

patterns = [
    "_spherepack*.so",
    "_spherepack*.pyd",
    "_spherepack*.dylib",
    "_spherepack_reduced*.so",
    "_spherepack_reduced*.pyd",
    "_spherepack_reduced*.dylib",
    "lib_spherepack_common*.so",
    "lib_spherepack_common*.dll",
    "lib_spherepack_common*.dylib",
    "_spherepack_common*.so",
    "_spherepack_common*.dll",
    "_spherepack_common*.dylib",
]
built_files = []

for pattern in patterns:
    found_files = glob.glob(str(current_dir / pattern), recursive=False)
    built_files.extend(found_files)

for pattern in patterns:
    found_files = glob.glob(str(current_dir / "**" / pattern), recursive=True)
    built_files.extend(found_files)

built_files = list(set(built_files))

if built_files:
    print(f"Found built files: {built_files}")

    stale_patterns = [
        "_spherepack.cp*-win_amd64.pyd",
        "_spherepack_direct.cp*-win_amd64.pyd",
        "_spherepack_reduced.cp*-win_amd64.pyd",
        "_spherepack*.so",
        "_spherepack*.dylib",
    ]
    keep_names = set()
    for stale_pattern in stale_patterns:
        for stale_file in module_dir.glob(stale_pattern):
            if stale_file.name in keep_names:
                continue
            try:
                stale_file.unlink()
                print(f"Removed stale source artifact {stale_file.name}")
            except Exception as e:
                print(f"Failed to remove stale source artifact {stale_file.name}: {e}")

    for stale_common in (
        list(module_dir.glob("_spherepack_common*.dll"))
        + list(module_dir.glob("_spherepack_common*.so"))
        + list(module_dir.glob("_spherepack_common*.dylib"))
        + list(module_dir.glob("lib_spherepack_common*.dll"))
        + list(module_dir.glob("lib_spherepack_common*.so"))
        + list(module_dir.glob("lib_spherepack_common*.dylib"))
    ):
        try:
            stale_common.unlink()
            print(f"Removed stale shared common artifact {stale_common.name}")
        except Exception as e:
            print(f"Failed to remove stale shared common artifact {stale_common.name}: {e}")

    for built_file_str in built_files:
        built_file = Path(built_file_str)
        if built_file.name.startswith("_spherepack_reduced"):
            continue
        source_file = module_dir / built_file.name

        should_copy = False
        if not source_file.exists():
            should_copy = True
            print(f"Source {built_file.name} doesn't exist, will copy")
        elif built_file.stat().st_mtime > source_file.stat().st_mtime:
            should_copy = True
            print(f"Source {built_file.name} is older, will copy")

        if should_copy:
            try:
                if source_file.exists():
                    source_file.unlink()
                shutil.copy2(str(built_file), str(source_file))
                print(f"Copied {built_file.name} to source directory successfully")
            except Exception as e:
                print(f"Failed to copy {built_file.name}: {e}")
else:
    print("No compiled extension files found")
    print("Available files in current directory:")
    if current_dir.exists():
        for f in current_dir.rglob("*"):
            if f.is_file():
                print(f"  {f}")
    else:
        print("  Current directory does not exist")

with open("''' + '@OUTPUT@' + '''", "w") as f:
    f.write("Source copy check completed\\n")

print("=== Extension source copy check completed ===")
'''
    ],
    build_by_default: true,
    console: true
)

message('spharm submodule configured successfully with cross-platform optimizations')
