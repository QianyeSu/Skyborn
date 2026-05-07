project(
  'skyborn-ectrans-backend',
  ['c', 'fortran'],
  version: '0.1.0',
  default_options: ['warning_level=2', 'buildtype=release', 'b_lto=false'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy_cmd = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
)
inc_np = include_directories(incdir_numpy_cmd.stdout().strip())

host_system = host_machine.system()
cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')

c_args = ['-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION']
fortran_args = []

if host_system == 'windows' and cc.get_id() == 'msvc'
  c_args += ['/O2']
else
  c_args += ['-O3']
endif

if fc.get_id() == 'gcc'
  fortran_args += ['-O3']
endif

ectrans_ext = py.extension_module(
  '_ectrans_backend',
  [
    '_ectrans_backend_module.c',
    'ectrans_gridmeta.f90',
    'ectrans_scalar.f90',
    'ectrans_vector.f90',
    'PARKIND1.F90',
    'YOMHOOK.F90',
    'TPM_CONSTANTS.F90',
    'TPM_DIM.F90',
    'TPM_DISTR.F90',
    'TPM_FIELDS.F90',
    'TPM_GEOMETRY.F90',
    'GSTATS.F90',
    'SUWAVEDI_MOD.F90',
    'PREPSNM_MOD.F90',
    'PRFI1B_MOD.F90',
    'LDFOU2_MOD.F90',
    'UPDSPB_MOD.F90',
    'VDTUV_MOD.F90',
    'UVTVD_MOD.F90',
    'VD2UV_MOD.F90',
    'VD2UV_CTL_MOD.F90',
    'runtime_state_mod.f90',
    'ldfou2_uv_scaling.f90',
    'vordiv_to_uv.f90',
    'uv_to_vordiv.f90',
    '../src/alf.f90',
    '../src/gaqd.f90',
    '../src/getlegfunc.f90',
    '../src/hrfft.f90',
    '../src/specintrp.f90',
    '../src/sphcom.f90',
    '../src/vhsgs.f90',
  ],
  include_directories: inc_np,
  dependencies: [py_dep],
  c_args: c_args,
  fortran_args: fortran_args,
  install: true,
  install_dir: get_option('python.purelibdir') / 'skyborn' / 'spharm',
  build_by_default: true,
)

copy_to_source = custom_target(
  'copy_ectrans_extension_to_source',
  input: ectrans_ext,
  output: 'inplace_copy_marker.txt',
  command: [
    py,
    '-c',
    '''
import glob
import shutil
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent.parent
patterns = ["_ectrans_backend*.so", "_ectrans_backend*.pyd", "_ectrans_backend*.dylib"]
built_files = []
for pattern in patterns:
    built_files.extend(glob.glob(str(current_dir / pattern), recursive=False))
    built_files.extend(glob.glob(str(current_dir / "**" / pattern), recursive=True))
built_files = list(set(built_files))
for built_file_str in built_files:
    built_file = Path(built_file_str)
    source_file = module_dir / built_file.name
    if (not source_file.exists()) or built_file.stat().st_mtime > source_file.stat().st_mtime:
        if source_file.exists():
            source_file.unlink()
        shutil.copy2(str(built_file), str(source_file))
with open("''' + '@OUTPUT@' + '''", "w", encoding="utf-8") as handle:
    handle.write("ectrans source copy completed\\n")
''',
  ],
  build_by_default: true,
  console: true,
)
