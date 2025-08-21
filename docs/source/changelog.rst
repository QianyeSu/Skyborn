Changelog
=========

Version 0.3.12.post1 (Current)
-------------------------------

**🔧 Critical Bug Fixes**

* **Fixed spharm Module Wheel Packaging**: Resolved critical issue where compiled Fortran extensions (``_spherepack*.pyd`` files) were missing from wheel distributions built via GitHub CI

  - **Root Cause**: Meson build system was installing to system paths instead of setuptools build directory
  - **Solution**: Configured meson to install directly to setuptools build directory using ``--python.purelibdir`` and ``--python.platlibdir`` parameters
  - **Impact**: Users can now install pre-compiled wheels with full spharm functionality

* **Improved Build System Integration**: Streamlined meson-setuptools integration for better maintainability

  - **Enhanced setup.py**: Added auto-discovery of meson modules for future extensibility
  - **Simplified Logic**: Removed complex file copying mechanisms in favor of native meson installation
  - **Better Error Handling**: Improved build process reliability across platforms

* **Fixed macOS Wheel Building**: Resolved OpenMP dependency compatibility issues

  - **Issue**: ``libgomp.1.dylib`` required minimum macOS 14.0 target version
  - **Solution**: Set ``MACOSX_DEPLOYMENT_TARGET=14.0`` in GitHub Actions workflow
  - **Note**: macOS 13 users can still install from source using ``pip install --no-binary=skyborn skyborn``

**🛠️ Technical Improvements**

* **Enhanced Meson Configuration**:

  - Changed from ``install: false`` to ``install: true`` with proper ``install_dir`` configuration
  - Maintained smart copying logic for ``--inplace`` builds
  - Enhanced cross-platform compatibility

* **Streamlined GitHub Actions**:

  - Updated wheel building workflow for better OpenMP library handling
  - Ensures compatibility with modern macOS development environments
  - Improved build reliability and error reporting

**✅ Validation**

* Successfully tested wheel building and installation across all supported platforms
* Confirmed ``_spherepack*.pyd`` files are correctly included in wheel distributions
* Verified functionality through comprehensive installation tests

Version 0.3.12
-------------------------------

**🚀 Major Performance Enhancements**

* **Modernized Spherical Harmonics (spharm) Submodule**: Complete Fortran code modernization for significantly improved windspharm performance:

  - **Modern Fortran Standards**: Updated legacy Fortran code to modern standards with improved memory management and vectorization
  - **~25% Performance Boost**: Windspharm calculations now run approximately 25% faster across all operations
  - **Optimized Algorithms**: Enhanced spherical harmonic transformations with better numerical efficiency
  - **Memory Optimization**: Improved memory layout and access patterns for better cache performance
  - **Cross-Platform Compatibility**: Better compiler optimization support across different platforms and architectures
  - **Maintained Accuracy**: All numerical results remain identical while achieving significant speed improvements

* **Enhanced Build System**: Streamlined compilation process for the modernized Fortran components

* **Python 3.13 Support**: Added full compatibility with Python 3.13:

  - **Wheel Distribution**: Pre-compiled wheels now available for Python 3.13 across all supported platforms
  - **Build System Compatibility**: Updated build configuration to support Python 3.13's new features and requirements
  - **Cross-Platform Testing**: Comprehensive testing on Linux x86_64, macOS (Intel & Apple Silicon), and Windows x64
  - **Future-Ready**: Ensures Skyborn stays current with the latest Python ecosystem developments

**🔧 Technical Improvements**

* **Fortran Modernization**:
  - Replaced obsolete Fortran constructs with modern equivalents
  - Improved array bounds checking and memory safety
  - Enhanced numerical stability in edge cases
  - Better integration with F2PY for Python bindings

* **Performance Optimizations**:
  - Vectorized mathematical operations in spherical harmonic calculations
  - Optimized Legendre polynomial computations
  - Reduced function call overhead in critical computation paths
  - Enhanced caching strategies for frequently used calculations

* **Platform and Build Improvements**:
  - **Extended Python Support**: Now supports Python 3.9, 3.10, 3.11, 3.12, and 3.13
  - **Multi-Platform Wheels**: Automated wheel building for Linux x86_64, macOS Intel/Apple Silicon, and Windows x64
  - **CI/CD Enhancements**: Improved build matrix with comprehensive testing across all supported Python versions
  - **Future ARM64 Linux Preparation**: Infrastructure ready for ARM64 Linux support when Python wheel ecosystem matures

**📊 Performance Benchmarks**

Windspharm operation speedups compared to previous version:
* **Vorticity Calculation**: ~25% faster execution time
* **Divergence Calculation**: ~25% faster execution time
* **Helmholtz Decomposition**: ~25% faster execution time
* **Streamfunction/Velocity Potential**: ~25% faster execution time
* **Combined Operations**: ~25% faster execution time

Version 0.3.11
-------------------------------

**🚀 Major Performance Improvements**

* **Optimized Mann-Kendall Trend Analysis**: Completely rewritten for significantly improved performance:

  - **Vectorized Implementation**: True vectorization of Mann-Kendall S-score calculation using advanced NumPy operations
  - **15-30x Performance Boost**: Processing speeds increased from ~19 to ~1,853 grid points per second for large climate datasets
  - **Climate Data Optimized**: Specifically tuned for typical climate data dimensions (40×192×288) with ~30-second processing time
  - **Memory Efficient**: Intelligent chunking strategy with only ~25MB memory usage for full climate grids
  - **Batch Processing**: Vectorized statistical calculations for clean data series, individual handling for series with missing values
  - **Enhanced Dask Support**: Improved map_blocks implementation for distributed computing workflows

* **Method Parameter Updates**: Replaced deprecated `method="auto"` with `method="theilslopes"` throughout the codebase for consistency

**🔧 Technical Improvements**

* **Simplified Import Structure**: Removed conditional/backup import logic in favor of direct scipy.stats imports for improved maintainability
* **Code Quality Enhancements**: Eliminated unused backup functions (`_mk_score_backup`, `_theil_sen_backup`) that were reducing test coverage
* **Consolidated Test Suite**: Merged supplementary test files into main test suite for better organization and reduced maintenance overhead
* **Documentation Fixes**: Corrected parameter names in API documentation examples (time_axis → axis)
* **Advanced Vectorization**: New `_vectorized_mk_score()` function using upper triangular indices for O(n²) to O(1) complexity reduction
* **Smart Memory Management**: Automatic chunk size estimation based on available memory and data dimensions
* **Robust Error Handling**: Graceful handling of edge cases and problematic time series
* **Comprehensive Testing**: Full test suite validation with 85% code coverage maintained

**🎨 UI/UX Improvements**

* **Dark Mode Compatibility**: Fixed notification color gradients for better visibility in dark themes:

  - Updated notification system to use deep blue to light blue gradient for improved contrast
  - Enhanced table responsiveness styling for better dark mode support

* **Documentation Accuracy**: Corrected function documentation to match actual codebase:

  - Fixed plot module function listings to reflect actual available functions
  - Removed non-existent functions from documentation (plot_field, plot_vector_field, plot_streamlines, plot_contour)
  - Added proper documentation for actual functions (add_equal_axes, createFigure, curved_quiver, add_curved_quiverkey)
  - Updated windspharm interface references for accurate Sphinx linking
  - Standardized "XArray" to "Xarray" throughout documentation

**📊 Performance Benchmarks**

For typical climate data analysis scenarios:

* **Small datasets** (50×20×30): 6.3x speedup (251 → 1,578 points/sec)
* **Medium datasets** (100×30×40): 14.8x speedup (74 → 1,093 points/sec)
* **Large datasets** (200×40×50): 31.3x speedup (19 → 595 points/sec)
* **Climate grids** (40×192×288): ~30 seconds total processing time

Version 0.3.10
-------------------------------

**🚀 New Features**

* **Advanced GridFill Module**: Major expansion of grid filling capabilities for atmospheric data interpolation:

  - **New XArray Interface**: Modern `skyborn.gridfill.xarray` module with automatic coordinate detection
  - **Comprehensive Tutorial**: Interactive Jupyter notebook demonstrating wind field gap filling techniques
  - **Multiple Interpolation Methods**: Basic Poisson, high-precision, zonal initialization, and relaxation parameter tuning
  - **Physical Validation**: Component-wise vs direct speed filling comparison for vector wind fields
  - **Quality Assessment**: Grid coverage validation and interpolation accuracy metrics

* **Rossby Wave Source Analysis**: Added comprehensive Rossby wave source calculation capabilities to the windspharm module:

  - New ``rossbywavesource()`` method in both standard and xarray interfaces
  - Implements the Sardeshmukh & Hoskins (1988) formulation: S = -ζₐ∇·v - v_χ·∇ζₐ
  - Support for custom truncation levels and Earth's angular velocity parameters
  - CF-compliant metadata for xarray output with proper units and standard names


**🔧 Improvements**

* **Test File Consolidation**: Merged duplicate gridfill test files for better maintainability
* **Better Grid Handling**: Improved spherical harmonic truncation validation for different grid sizes
* **Documentation Updates**: Enhanced gallery with new Rossby wave source visualization examples

**📚 Documentation**

* **New GridFill Tutorial**: Complete interactive demonstration including:

  - Advanced data interpolation techniques with real atmospheric wind data
  - Missing data simulation and quality assessment methodologies
  - Component-wise vs direct approach comparison for vector fields
  - Publication-quality visualizations with integer colorbar formatting
  - Performance analysis and best practices for atmospheric applications

* **New Tutorial Notebooks**: Added comprehensive examples for:

  - Rossby wave source analysis and visualization
  - Grid filling techniques with atmospheric data
  - Longitude coordinate system transformations

* **Enhanced Gallery**: Updated with new visualization examples including:

  - ``windspharm_rossby_wave_source_truncations.png`` showing truncation effects
  - ``gridfill_missing_data_overview.png`` demonstrating gap filling scenarios
  - ``gridfill_component_vs_direct_comparison.png`` showing physical constraint preservation
  - Improved figure captions and mathematical formulations
  - Better integration of notebook examples

**🧪 Testing**

* **Expanded Test Coverage**: Added comprehensive tests for new Rossby wave source functionality
* **Grid Size Validation**: Enhanced parameter validation for different grid resolutions
* **Cross-interface Testing**: Verified consistency between standard and xarray interfaces

**Technical Notes**

* All existing functionality remains backward compatible
* Enhanced error handling for grid size limitations in spherical harmonic calculations
* Improved memory efficiency for large-scale atmospheric analysis

Version 0.3.9
------------------------

**New Features**

* **Enhanced Spherical Harmonics Module**: Improved performance and stability for atmospheric data analysis
* **New Windspharm Submodule**: Added comprehensive wind field analysis capabilities including:

  - Vector wind analysis and spherical harmonic transforms
  - Vorticity and divergence calculations
  - Stream function and velocity potential computations
  - Compatible with various grid types and coordinate systems

* **Optimized Build System**: Streamlined compilation process for better cross-platform compatibility

**🔧 Improvements**

* **Better Error Handling**: Enhanced error messages and debugging information
* **Performance Optimizations**: Faster execution for large-scale atmospheric calculations
* **Code Quality**: Improved type hints and documentation coverage

**🐛 Bug Fixes**

* **Fixed Dimension Handling in Regridding**: Resolved dimension change issues in interp.regridding.py module that were causing inconsistent array shapes during interpolation operations
* Fixed interpolation edge cases in atmospheric data processing
* Resolved compilation issues on various platforms
* Improved numerical stability in spherical harmonic transforms

**📚 Documentation**

* **Windspharm Module Documentation**: Complete documentation and examples for wind field analysis functions
* Added comprehensive examples and tutorials
* Enhanced API reference with mathematical formulations
* Improved installation and usage guides

**🔧 Technical Details**

* **Dependencies**: Updated NumPy compatibility, enhanced F2PY integration, improved Fortran compiler support
* **Platform Support**: Linux x86_64 (manylinux2014), macOS (Intel and Apple Silicon), Windows x64
* **Windspharm Dependencies**: Added support for spherical harmonic wind analysis libraries

Version 0.3.8
--------------

**🔧 Bug Fixes**

* **fix**: remove obsolete Fortran wrapper file spherepack-f2pywrappers.f
* Improved build system stability and cross-platform compatibility
* Enhanced error handling and debugging information

**📚 Documentation**

* Updated API documentation
* Improved code examples and installation guides
* Enhanced cross-reference documentation

Version 0.3.7
--------------

**✨ New Features**

* **Emergent Constraints Method**: Added new emergent constraints analysis method for climate data analysis
* **Enhanced Documentation**: Interactive particle effects entrance page

**🔧 Improvements**

* Optimized documentation structure and user interface
* Updated interactive documentation entry page with particle effects
* Improved cross-platform compatibility
* Enhanced code quality and test coverage

**📚 Documentation**

* New particle effects documentation entrance page
* Updated API documentation
* Improved code examples and usage guides
* Enhanced Sphinx Book Theme with blue color scheme

**🐛 Bug Fixes**

* Fixed minor issues and improved code quality
* Resolved documentation build issues
* Enhanced error handling

Version 0.3.6
--------------

* Added emergent constraint analysis functionality
* Improved GRIB to NetCDF conversion
* Added comprehensive documentation with Jupyter notebooks
* Enhanced statistical analysis functions
