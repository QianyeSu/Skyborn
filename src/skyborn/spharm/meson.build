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

# Get numpy include directories with robust fallbacks
incdir_numpy_cmd = run_command(py,
  ['-c', '''
import numpy
import os
# Try multiple potential paths
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
    print(numpy.get_include())  # fallback
'''],
  check : true
)
incdir_numpy = incdir_numpy_cmd.stdout().strip()

# Try to get f2py include, with multiple fallbacks
incdir_f2py_cmd = run_command(py,
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
    check : true
)
incdir_f2py = incdir_f2py_cmd.stdout().strip()

# Debug output for CI troubleshooting
message('NumPy include directory: ' + incdir_numpy)
message('F2PY include directory: ' + incdir_f2py)

# Generate Fortran-Python interface with better error handling
spherepack_source = custom_target('_spherepackmodule.c',
    input: ['src/_spherepack.pyf'],
    output: ['_spherepackmodule.c', '_spherepack-f2pywrappers.f'],
    command: [py, '-m', 'numpy.f2py', '--lower', '@INPUT@', '--build-dir', '@OUTDIR@'],
    build_by_default: true,
    install: false
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Verify numpy headers exist (for CI debugging)
numpy_header_check = run_command(py,
  ['-c', 'import os; import numpy; header_path = os.path.join(numpy.get_include(), "numpy", "npy_os.h"); print("✓" if os.path.exists(header_path) else "✗", header_path)'],
  check: false
)
message('NumPy header check: ' + numpy_header_check.stdout().strip())

# Robust fortran compilation flags for different environments
fortran_args = ['-O3', '-fPIC', '-fno-second-underscore']
c_args = ['-O3', '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']

# Additional debugging for CI
message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

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
