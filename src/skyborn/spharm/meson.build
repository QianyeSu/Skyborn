# spharm submodule for skyborn - OPTIMIZED VERSION
# Spherical harmonic transforms for atmospheric/oceanic modeling
# Performance optimizations: Modern Fortran, OpenMP, vectorization

# Install Python module files - updated structure
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

# Fortran compiler optimization flags
fortran_args = []
c_args = []

# Get optimization mode from environment (portable vs native)
optimization_mode = get_env('OPTIMIZATION_MODE', 'native')

# GCC/gfortran optimizations
if meson.get_compiler('fortran').get_id() == 'gcc'
    if optimization_mode == 'portable'
        # Portable optimization for wheels
        fortran_args += [
            '-O3',                    # Maximum optimization
            '-march=x86-64',          # Portable x86-64 baseline
            '-mtune=generic',         # Generic tuning
            '-funroll-loops',         # Aggressive loop unrolling
            '-ffast-math',            # Fast math optimizations
            '-ftree-vectorize',       # Enable vectorization
            '-fopenmp',               # OpenMP support
            '-fPIC',                  # Position independent code
            '-fno-second-underscore'  # F2PY compatibility
        ]
        c_args += [
            '-O3',
            '-march=x86-64',
            '-mtune=generic',
            '-fopenmp',
            '-ffast-math'
        ]
    else
        # Native optimization for source builds
        fortran_args += [
            '-O3',                    # Maximum optimization
            '-march=native',          # Optimize for current CPU
            '-mtune=native',         # Tune for current CPU
            '-funroll-loops',        # Aggressive loop unrolling
            '-ffast-math',           # Fast math optimizations
            '-ftree-vectorize',      # Enable vectorization
            '-fopenmp',              # OpenMP support
            '-fPIC',                 # Position independent code
            '-fno-second-underscore' # F2PY compatibility
        ]
        c_args += [
            '-O3',
            '-march=native',
            '-mtune=native',
            '-fopenmp',
            '-ffast-math'
        ]
    endif
endif

# Intel Fortran optimizations (if available)
if meson.get_compiler('fortran').get_id() == 'intel'
    if optimization_mode == 'portable'
        # Portable optimization for wheels
        fortran_args += [
            '-O3',
            '-march=x86-64',         # Portable x86-64 baseline
            '-ipo',                  # Interprocedural optimization
            '-fp-model', 'fast=2',   # Fast floating point
            '-fma',                  # Fused multiply-add
            '-qopenmp',              # OpenMP support
            '-fPIC'
        ]
        c_args += [
            '-O3',
            '-march=x86-64',
            '-ipo',
            '-fp-model', 'fast=2',
            '-qopenmp'
        ]
    else
        # Native optimization for source builds
        fortran_args += [
            '-O3',
            '-xHost',               # Optimize for current CPU
            '-ipo',                 # Interprocedural optimization
            '-fp-model', 'fast=2',  # Fast floating point
            '-fma',                 # Fused multiply-add
            '-qopenmp',             # OpenMP support
            '-fPIC'
        ]
        c_args += [
            '-O3',
            '-xHost',
            '-ipo',
            '-fp-model', 'fast=2',
            '-qopenmp'
        ]
    endif
endif

# Build the optimized spherepack extension module
py.extension_module('_spherepack',
    [
        # Use optimized .f90 files with proper dependency order
        'src/sphcom.f90',               # First: common module that others depend on
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
        # SPHEREPACK core library sources - optimized .f90 files
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
    incdir_f2py / 'fortranobject.c',
    include_directories: inc_np,
    dependencies: [py_dep, numpy_dep],
    fortran_args: fortran_args,
    c_args: c_args,
    install: true,
    subdir: 'skyborn/spharm'
)

message('spharm submodule configured with FULLY OPTIMIZED Fortran extensions')
message('Optimization level: O3 with native CPU optimizations')
message('Features: OpenMP parallelization, vectorization, fast math')
message('ALL FORTRAN FILES OPTIMIZED: Modern Fortran (.f90), parallel computing, vectorization')
message('Optimized functions: vorticity/divergence, gradients, Legendre, FFT, Gaussian quadrature, spherical harmonics')
