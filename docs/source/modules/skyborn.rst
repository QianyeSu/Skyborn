Skyborn Modules
===============

This section provides an overview of all Skyborn modules and their functionality.

Core Modules
------------

The Skyborn package contains the following main modules:

* **calc**: Statistical calculations and emergent constraint methods
* **conversion**: Data format conversion utilities (GRIB to NetCDF)
* **interp**: Interpolation and regridding functions
* **plot**: Specialized plotting and visualization tools
* **gradients**: Spatial and temporal gradient calculations
* **causality**: Causal analysis methods (Granger, Liang)
* **ROF**: Regularized Optimal Fingerprinting (ROF) method; currently under development and testing.

Module Overview
---------------

**calc**
   Statistical calculations and emergent constraint methods

**conversion**
   Data format conversion utilities (GRIB to NetCDF)

**interp**
   Interpolation and regridding functions

**plot**
   Specialized atmospheric data visualization

**gradients**
   Spatial and temporal gradient calculations

**causality**
   Granger and Liang causality analysis

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

   # Create visualizations
   fig, ax = skb.plot.createFigure((10, 6), 1, 1)

For detailed API documentation, see :doc:`../api/index`.
