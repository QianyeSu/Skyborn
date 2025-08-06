Changelog
=========

Version 0.3.9 (Current)
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
