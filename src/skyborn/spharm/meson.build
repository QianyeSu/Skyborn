# spharm submodule for skyborn - Spherical harmonic transforms
# Simplified configuration for better compatibility

# Get Python installation
py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Install Python module files
py.install_sources([
    '__init__.py',
    'spharm.py',
  ],
  subdir: 'skyborn/spharm'
)

# Get numpy include directories
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

# Generate Fortran-Python interface
spherepack_source = custom_target('_spherepackmodule.c',
    input: ['src/_spherepack.pyf'],
    output: ['_spherepackmodule.c', '_spherepack-f2pywrappers.f'],
    command: [py, '-m', 'numpy.f2py', '--lower', '@INPUT@', '--build-dir', '@OUTDIR@'],
    build_by_default: true
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Simple fortran compilation flags
fortran_args = ['-O3', '-fPIC', '-fno-second-underscore']
c_args = ['-O3']

# Build the spherepack extension module
py.extension_module('_spherepack',
    [
        # Core SPHEREPACK sources
        'src/sphcom.f90',
        'src/getlegfunc.f90',
        'src/specintrp.f90',
        'src/onedtotwod.f90',
        'src/onedtotwod_vrtdiv.f90',
        'src/twodtooned.f90',
        'src/twodtooned_vrtdiv.f90',
        'src/multsmoothfact.f90',
        'src/lap.f90',
        'src/invlap.f90',
        'src/gaqd.f90',
        'src/shses.f90',
        'src/shaes.f90',
        'src/vhaes.f90',
        'src/vhses.f90',
        'src/shsgs.f90',
        'src/shags.f90',
        'src/vhags.f90',
        'src/vhsgs.f90',
        'src/hrfft.f90',
        'src/shaec.f90',
        'src/shagc.f90',
        'src/shsec.f90',
        'src/shsgc.f90',
        'src/vhaec.f90',
        'src/vhagc.f90',
        'src/vhsec.f90',
        'src/vhsgc.f90',
        'src/ihgeod.f90',
        'src/alf.f90',
        spherepack_source
    ],
    include_directories: inc_np,
    dependencies: py_dep,
    fortran_args: fortran_args,
    c_args: c_args,
    install: true,
    subdir: 'skyborn/spharm'
)

message('spharm submodule configured')
