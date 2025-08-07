Skyborn Modules
===============

This section provides an overview of all Skyborn modules and their functionality.

Core Modules
------------

The Skyborn package contains the following main modules:

* **calc**: Statistical calculations and emergent constraint methods
* **conversion**: Data format conversion utilities (GRIB to NetCDF)
* **gridfill**: Advanced data interpolation using Poisson equation solvers
* **interp**: Interpolation and regridding functions
* **plot**: Specialized plotting and visualization tools
* **gradients**: Spatial and temporal gradient calculations
* **causality**: Causal analysis methods (Granger, Liang)
* **spharm**: Spherical harmonic transforms and analysis (New in v0.3.8+)
* **windspharm**: Vector wind field analysis on the sphere (New in v0.3.8+)
* **ROF**: Regularized Optimal Fingerprinting (ROF) method; currently under development and testing.

Module Overview
---------------

**calc**
   Statistical calculations and emergent constraint methods for climate data analysis

**conversion**
   Data format conversion utilities (GRIB to NetCDF) with enhanced error handling

**gridfill** *(New in v0.3.10)*
   Advanced data interpolation using Poisson equation solvers for filling missing values in gridded data.
   Provides multiple interfaces (standard, iris, xarray) and algorithms for atmospheric and oceanic applications.

**interp**
   Interpolation and regridding functions with improved dimension handling

**plot**
   Specialized atmospheric data visualization with curved quiver plots

**gradients**
   Spatial and temporal gradient calculations for atmospheric fields

**causality**
   Granger and Liang causality analysis for climate variable relationships

**spharm** *(New in v0.3.8+)*
   Spherical harmonic transforms and analysis for global atmospheric data.
   Provides fast spectral transforms, filtering, and grid conversions.

**windspharm** *(New in v0.3.8+)*
   Comprehensive vector wind field analysis on the sphere including:

   - Vorticity and divergence calculations
   - Stream function and velocity potential computations
   - Helmholtz decomposition (rotational/irrotational components)
   - Spectral truncation and filtering
   - Multiple interface support (xarray, standard, tools)

**ROF**
   Regularized Optimal Fingerprinting (ROF) method; currently under development and testing.

Getting Started
---------------

To use any module, import Skyborn and access the module:

.. code-block:: python

   import skyborn as skb

   # Use calculation functions
   pdf = skb.calc.gaussian_pdf(mu=0, sigma=1, x=x_values)

   # Convert data formats
   skb.conversion.grib_to_netcdf('input.grib', 'output.nc')

   # Fill missing data with GridFill
   from skyborn.gridfill.xarray import fill
   filled_data = fill(data_with_gaps, eps=1e-4)

   # Create visualizations
   fig, ax = skb.plot.createFigure((10, 6), 1, 1)

   # NEW: Use spherical harmonic analysis
   from skyborn.spharm import Spharmt
   spharm = Spharmt(nlon=144, nlat=72, gridtype='gaussian')

   # NEW: Analyze wind fields
   from skyborn.windspharm.xarray import VectorWind
   vw = VectorWind(u_wind, v_wind)
   vorticity = vw.vorticity()
   divergence = vw.divergence()

**Version Highlights**

* **v0.3.7**: Enhanced documentation with interactive particle effects entrance page
* **v0.3.8**: Added comprehensive windspharm module and improved spherical harmonics support
* **v0.3.9**: Current version with optimized build system and enhanced cross-platform compatibility
* **v0.3.10**: Major GridFill module expansion with advanced interpolation capabilities and comprehensive tutorial

For detailed API documentation, see :doc:`../api/index`.

**Quick Links**

* :doc:`../notebooks/gridfill_tutorial` - Complete tutorial for advanced data interpolation
* :doc:`../notebooks/windspharm_tutorial` - Complete tutorial for wind field analysis
* :doc:`../api/gridfill` - GridFill data interpolation documentation
* :doc:`../api/spharm` - Spherical harmonic transforms documentation
* :doc:`../api/windspharm` - Vector wind analysis documentation
