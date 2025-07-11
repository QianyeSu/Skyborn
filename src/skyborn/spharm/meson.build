# spharm submodule for skyborn - Spherical harmonic transforms
# Simplified configuration for better compatibility

# Get Python installation - use the same as parent
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Install Python module files
py.install_sources([
    '__init__.py',
    'spherical_harmonics.py',  # Renamed from spharm.py to avoid naming conflicts
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
spherepack_sources = []
try_f2py = run_command(py, '-c', 'import numpy.f2py; print("OK")', check: false)

if try_f2py.returncode() == 0
  # Define variables for successful f2py case
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

  # Get fortranobject.c path with extensive fallbacks for CI environments
  fortranobject_cmd = run_command(py,
    ['-c', '''
import numpy.f2py
import numpy
import os
import sys

# Try multiple potential paths for fortranobject.c
candidates = []

# Method 1: f2py.get_include()
try:
    fo_path = os.path.join(numpy.f2py.get_include(), "fortranobject.c")
    candidates.append(fo_path)
except:
    pass

# Method 2: numpy package paths
numpy_base = numpy.__path__[0]
candidates.extend([
    os.path.join(numpy_base, "f2py", "src", "fortranobject.c"),
    os.path.join(numpy_base, "f2py", "fortranobject.c"),
    os.path.join(numpy_base, "_core", "src", "fortranobject.c"),
])

# Method 3: site-packages paths (for pip installs)
try:
    import site
    for site_dir in site.getsitepackages() + [site.getusersitepackages()]:
        if site_dir:
            candidates.extend([
                os.path.join(site_dir, "numpy", "f2py", "src", "fortranobject.c"),
                os.path.join(site_dir, "numpy", "f2py", "fortranobject.c"),
            ])
except:
    pass

# Method 4: sys.path based search
for path in sys.path:
    if path and "numpy" in path:
        candidates.extend([
            os.path.join(path, "numpy", "f2py", "src", "fortranobject.c"),
            os.path.join(path, "numpy", "f2py", "fortranobject.c"),
        ])

# Find first existing file
for candidate in candidates:
    if os.path.exists(candidate):
        print(candidate)
        break
else:
    print("")  # empty if not found
'''],
    check: true
  )
  fortranobject_c = fortranobject_cmd.stdout().strip()

  # Try to generate the interface using the core working functions only
  spherepack_source = custom_target('_spherepackmodule.c',
      input: ['src/_spherepack_core.pyf'],
      output: ['_spherepackmodule.c', '_spherepack-f2pywrappers2.f90'],
      command: [py, '-m', 'numpy.f2py', '--lower', '@INPUT@', '--build-dir', '@OUTDIR@'],
      build_by_default: true,
      install: false
  )
  spherepack_sources += spherepack_source

  # Add fortranobject.c if found, otherwise try to build without it
  if fortranobject_c != ''
    message('Found fortranobject.c: ' + fortranobject_c)
    # Use Python to check file existence (works on all platforms)
    fortran_obj_check = run_command(py, '-c', 'import os; import sys; sys.exit(0 if os.path.exists(r"' + fortranobject_c + '") else 1)', check: false)
    if fortran_obj_check.returncode() == 0
      spherepack_sources += files(fortranobject_c)
      message('Added fortranobject.c to build')
    else
      warning('fortranobject.c path found but file does not exist: ' + fortranobject_c)
    endif
  else
    warning('fortranobject.c not found - trying to build without it (may cause issues)')
  endif

  message('f2py interface generation: SUCCESS')

  # Build the spherepack extension module with modernized .f90 files
  py.extension_module('_spherepack',
      [
          # Core SPHEREPACK sources - using .f90 modernized files
          'src/sphcom.f90',
          'src/hrfft.f90',
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
          spherepack_sources
      ],
      include_directories: inc_np,
      dependencies: py_dep,
      fortran_args: fortran_args,
      c_args: c_args,
      install: true,
      subdir: 'skyborn/spharm'
  )

  message('spharm submodule configured successfully')
else
  warning('f2py not available, skipping spharm interface generation')
  # Still install Python files even if Fortran extension fails
  message('Installing Python-only spharm module (without Fortran extensions)')
endif
