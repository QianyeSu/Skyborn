Interpolation
=======================

The Skyborn interpolation module provides data interpolation and regridding capabilities for atmospheric and climate data.

Regridding Classes
------------------

.. autoclass:: skyborn.interp.Grid
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: skyborn.interp.Regridder
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: skyborn.interp.NearestRegridder
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: skyborn.interp.BilinearRegridder
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: skyborn.interp.ConservativeRegridder
   :members:
   :undoc-members:
   :show-inheritance:


Regridding Functions
--------------------

.. autofunction:: skyborn.interp.regrid_dataset

.. autofunction:: skyborn.interp.nearest_neighbor_indices

Interpolation Functions
-----------------------

.. autofunction:: skyborn.interp.interp_hybrid_to_pressure

.. autofunction:: skyborn.interp.interp_sigma_to_hybrid

.. autofunction:: skyborn.interp.interp_multidim


Curvilinear and Unstructured Interpolation
------------------------------------------

.. autofunction:: skyborn.interp.rcm2points

.. autofunction:: skyborn.interp.rcm2rgrid

.. autofunction:: skyborn.interp.rgrid2rcm

.. autofunction:: skyborn.interp.grid_to_triple

.. autofunction:: skyborn.interp.triple_to_grid


Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import xarray as xr

   # Load data
   data = xr.open_dataset('input_data.nc')

   # Create regridder
   regridder = skb.interp.ConservativeRegridder(
       source_grid=data,
       target_grid=target_grid
   )

   # Regrid data
   regridded_data = regridder.regrid_array(data['temperature'])

   # Convert gridded data to triples (x, y, value) and back to a rectilinear grid
   from skyborn.interp import grid_to_triple, triple_to_grid
   triples = grid_to_triple(data['temperature'], data['lon'], data['lat'])
   # Define target 1D coordinates
   x_out = xr.DataArray(np.linspace(float(data['lon'].min()), float(data['lon'].max()), 100), dims=('x',))
   y_out = xr.DataArray(np.linspace(float(data['lat'].min()), float(data['lat'].max()), 80), dims=('y',))
   gridded = triple_to_grid(triples[2], triples[0], triples[1], x_out, y_out, method=1, domain=1.0)
