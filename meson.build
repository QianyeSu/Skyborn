project('skyborn', 'c', 'fortran',
  version: '0.3.7',
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
subdir('src/skyborn/spharm')

# Summary
message('Skyborn configuration complete!')
message('Version: ' + meson.project_version())
message('Build type: ' + get_option('buildtype'))
