Gradient Functions
==================

Spatial and temporal gradient calculation utilities for atmospheric and climate data.

Spatial Gradients
-----------------

.. autofunction:: skyborn.gradients.calculate_gradient
   :no-index:

.. autofunction:: skyborn.gradients.calculate_meridional_gradient
   :no-index:

.. autofunction:: skyborn.gradients.calculate_zonal_gradient
   :no-index:

.. autofunction:: skyborn.gradients.calculate_vertical_gradient
   :no-index:

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import xarray as xr

   # Load atmospheric data
   data = xr.open_dataset('temperature_data.nc')
   temp = data['temperature']

   # Calculate meridional temperature gradient
   meridional_grad = skb.calculate_meridional_gradient(temp)

   # Calculate zonal temperature gradient
   zonal_grad = skb.calculate_zonal_gradient(temp)

   # Calculate vertical gradient (if pressure levels available)
   vertical_grad = skb.calculate_vertical_gradient(temp, data['pressure'])
       data['temperature'],
       data['latitude']
   )

   # Calculate temporal gradient
   temp_grad = skb.temporal_gradient(
       data['temperature'],
       data['time']
   )
