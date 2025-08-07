Changelog
=========

Version 0.3.10 (Current)
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
