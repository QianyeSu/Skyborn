Interpolation
=======================

The Skyborn interpolation module provides data interpolation and regridding capabilities for atmospheric and climate data.

Regridding Classes
------------------

.. autoclass:: skyborn.interp.Grid
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. autoclass:: skyborn.interp.Regridder
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. autoclass:: skyborn.interp.NearestRegridder
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. autoclass:: skyborn.interp.BilinearRegridder
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. autoclass:: skyborn.interp.ConservativeRegridder
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

Regridding Functions
--------------------

.. autofunction:: skyborn.interp.regrid_dataset
   :no-index:

.. autofunction:: skyborn.interp.nearest_neighbor_indices
   :no-index:

Interpolation Functions
-----------------------

.. autofunction:: skyborn.interp.interp_hybrid_to_pressure
   :no-index:

.. autofunction:: skyborn.interp.interp_sigma_to_hybrid
   :no-index:

.. autofunction:: skyborn.interp.interp_multidim
   :no-index:

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
