project('skyborn', 'c', 'fortran',
  version: '0.3.8',
  license: 'BSD-3-Clause',
  meson_version: '>= 0.64.0',
  default_options: [
    'warning_level=2',
    'buildtype=release',
    'optimization=3',
  ]
)

# Import Python functionality
py_mod = import('python')
py = py_mod.find_installation(pure: false)
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
if not spharm_enabled.disabled()
  try_spharm = run_command(py,
    ['-c', 'import numpy; import numpy.f2py; print("OK")'],
    check: false
  )
  if try_spharm.returncode() == 0
    message('Building spharm module with numpy/f2py support')
    subdir('src/skyborn/spharm')
  else
    if spharm_enabled.enabled()
      warning('Failed to detect numpy/f2py for spharm module, but continuing...')
    else
      message('Skipping spharm module due to missing dependencies')
    endif
  endif
endif

# Summary
message('Skyborn configuration complete!')
message('Version: ' + meson.project_version())
message('Build type: ' + get_option('buildtype'))
