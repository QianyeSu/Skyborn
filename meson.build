project('skyborn', 'c',
  version: '0.3.8',
  license: 'BSD-3-Clause',
  meson_version: '>= 0.64.0',
  default_options: [
    'warning_level=2',
    'buildtype=release',
    'optimization=3',
  ]
)

# Check for Fortran compiler availability with better detection
fortran_available = false
fortran_compiler = find_program('gfortran', required: false)
if not fortran_compiler.found()
  # Try alternative Fortran compilers
  fortran_compiler = find_program('ifort', required: false)
  if not fortran_compiler.found()
    fortran_compiler = find_program('ifx', required: false)
  endif
endif

if fortran_compiler.found()
  # Try to add Fortran language support
  add_languages('fortran', required: false, native: false)
  fortran_available = true
  message('Fortran compiler found: building with Fortran extensions')
else
  message('Fortran compiler not found: building without Fortran extensions')
  message('Users can still install skyborn, but spharm module will not be available')
endif

# Import Python functionality - explicitly use current python
py_mod = import('python')
py = py_mod.find_installation('python', pure: false)
py_dep = py.dependency()

# Check for required Python packages using run_command
numpy_dep = declare_dependency(
  dependencies: py_dep
)

# Add subdirectories
subdir('src')

# Install Python package files
py.install_sources(
  [
    'src/skyborn/__init__.py',
    'src/skyborn/causality.py',
    'src/skyborn/gradients.py',
  ],
  subdir: 'skyborn'
)

# Install submodules
subdir('src/skyborn/calc')
subdir('src/skyborn/conversion')
subdir('src/skyborn/interp')
subdir('src/skyborn/plot')
subdir('src/skyborn/ROF')

# Try to build spharm with better error handling
spharm_enabled = get_option('spharm')
if not spharm_enabled.disabled() and fortran_available
  # Check if we can build spharm
  can_build_spharm = true

  # Test numpy/f2py availability
  try_spharm = run_command(py,
    ['-c', 'import numpy; import numpy.f2py; print("OK")'],
    check: false
  )

  if try_spharm.returncode() != 0
    can_build_spharm = false
    message('WARNING: numpy/f2py not available for spharm')
  endif

  # Test numpy headers
  numpy_check = run_command(py,
    ['-c', 'import numpy; import os; print("OK" if os.path.exists(os.path.join(numpy.get_include(), "numpy", "npy_os.h")) else "FAIL")'],
    check: false
  )

  if numpy_check.stdout().strip() != 'OK'
    can_build_spharm = false
    message('WARNING: numpy headers not found for spharm')
  endif

  if can_build_spharm
    message('Building spharm module with numpy/f2py support')
    subdir('src/skyborn/spharm')
  else
    message('Skipping spharm module due to missing dependencies')
  endif
elif not fortran_available
  message('Skipping spharm module: Fortran compiler not available')
endif

# Summary
message('Skyborn configuration complete!')
message('Version: ' + meson.project_version())
message('Build type: ' + get_option('buildtype'))
