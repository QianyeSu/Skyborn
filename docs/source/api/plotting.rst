Plotting
==================

The Skyborn plotting module provides specialized visualization functions for atmospheric and climate data.

Plotting Utilities
------------------

.. autofunction:: skyborn.plot.add_equal_axes

.. autofunction:: skyborn.plot.createFigure

Curly Vector Plots
------------------

.. autofunction:: skyborn.plot.curly_vector

The same ``skyborn.plot.curly_vector`` entry point supports both call styles:

- ``curly_vector(ax, x, y, u, v, ...)`` for direct NumPy / Matplotlib input
- ``curly_vector(ds, x="lon", y="lat", u="u", v="v", ax=ax, ...)`` for xarray dataset input

.. autofunction:: skyborn.plot.curly_vector_key

.. autoclass:: skyborn.plot.CurlyVectorPlotSet
   :members:

.. autoclass:: skyborn.plot.CurlyVectorKey
   :members:

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import matplotlib.pyplot as plt
   import numpy as np
   import xarray as xr

   # Create figure with equal axes
   fig, ax = plt.subplots()
   ax_new = skb.plot.add_equal_axes(ax, 'right', 0.1, 0.2)

   # Low-level array API
   lon = np.linspace(-180, 180, 50)
   lat = np.linspace(-90, 90, 25)
   u = np.random.random((25, 50))
   v = np.random.random((25, 50))

   skb.plot.curly_vector(ax, lon, lat, u, v)

   # Dataset wrapper API
   ds = xr.Dataset(
       {"u": (("lat", "lon"), u), "v": (("lat", "lon"), v)},
       coords={"lon": lon, "lat": lat},
   )
   skb.plot.curly_vector(ds, x="lon", y="lat", u="u", v="v", ax=ax)

Visualization Examples
----------------------

The plotting module is designed for creating publication-quality atmospheric visualizations:

.. code-block:: python

   import skyborn as skb
   import matplotlib.pyplot as plt

   # Create figure
   fig = skb.plot.createFigure(figsize=(12, 8), dpi=300)
   ax = fig.add_subplot(111)

   # Plot wind vectors from arrays
   curly = skb.plot.curly_vector(ax, lon, lat, u_wind, v_wind)

For complete visualization examples, see :doc:`../gallery`.
