project(
    'skyborn-plot-nclcurly-native',
    'c',
    version: '1.0.0',
    default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

py.install_sources(
    ['__init__.py'],
    subdir: 'skyborn/plot/nclcurly_native',
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
incdir_numpy = incdir_numpy_cmd.stdout().strip()
inc_np = include_directories(incdir_numpy)

host_system = host_machine.system()
cc = meson.get_compiler('c')

c_args = [
    '-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION',
]

if host_system == 'windows' and cc.get_id() == 'msvc'
    c_args += ['/O2', '/GL', '/Oi', '/Ot']
else
    c_args += ['-O3', '-fno-math-errno', '-fPIC']
endif

nclcurly_ext = py.extension_module(
    '_ncl_curly_core',
    '_ncl_curly_core.c',
    include_directories: inc_np,
    dependencies: [py_dep],
    c_args: c_args,
    install: true,
    install_dir: get_option('python.purelibdir') / 'skyborn' / 'plot' / 'nclcurly_native',
    build_by_default: true,
)

copy_to_source = custom_target(
    'copy_to_source_for_inplace',
    input: nclcurly_ext,
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
patterns = ["_ncl_curly_core*.so", "_ncl_curly_core*.pyd", "_ncl_curly_core*.dylib"]
built_files = []

for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / pattern), recursive=False))
for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / "**" / pattern), recursive=True))

for built_file_str in sorted(set(built_files)):
    built_file = Path(built_file_str)
    target = module_dir / built_file.name
    should_copy = (not target.exists()) or (built_file.stat().st_mtime > target.stat().st_mtime)
    if should_copy:
        if target.exists():
            target.unlink()
        shutil.copy2(str(built_file), str(target))

with open("''' + '@OUTPUT@' + '''", "w", encoding="utf-8") as marker:
    marker.write("nclcurly native inplace copy completed\\n")
''',
    ],
    build_by_default: true,
    console: true,
)
