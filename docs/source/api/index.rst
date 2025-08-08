API Reference
=============

.. toctree::
   :maxdepth: 2

   calculations
   mann_kendall
   conversion
   plotting
   interpolation
   gridfill
   gradients
   causality
   spharm
   windspharm

This section provides detailed API documentation for all Skyborn modules.

Calculation Functions
---------------------

Statistical analysis and emergent constraint methods.

Mann-Kendall Trend Analysis
---------------------------

High-performance trend detection for climate data with vectorized implementations optimized for large multidimensional arrays.

**Version 0.3.11 Highlights:**
- 15-30x performance improvements through vectorization
- Climate data optimized (40×192×288 processed in ~30 seconds)
- Memory efficient chunking (~25MB usage)
- Enhanced Dask support for distributed computing

Conversion Functions
--------------------

Data format conversion utilities.

Plotting Functions
------------------

Visualization and plotting utilities.

Interpolation Functions
-----------------------

Data interpolation and regridding.

GridFill Functions
------------------

Advanced atmospheric data interpolation using Poisson equation solvers. Provides gap-filling and smoothing capabilities for irregular atmospheric datasets with physically-based methods.

**Version 0.3.10 Highlights:**
- Enhanced xarray interface with automatic coordinate detection
- Comprehensive tutorial and examples
- Improved performance and convergence handling
- Publication-ready visualization tools

Gradient Functions
------------------

Spatial and temporal gradient calculations.

Causality Functions
-------------------

Causality analysis methods.

Spherical Harmonics Functions
-----------------------------

Spherical harmonic transforms and spectral analysis for atmospheric and oceanic data.

Vector Wind Analysis Functions
------------------------------

Spherical harmonic vector wind analysis including vorticity, divergence, streamfunction, and velocity potential calculations.

.. note::
   The ROF (Regularized Optimal Fingerprinting) module; currently under development and testing.
   Documentation will be available in future releases.

ROF Functions
-------------

Reduced Order Form methods for climate model analysis.
