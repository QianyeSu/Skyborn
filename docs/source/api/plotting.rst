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

The recommended module-oriented import path for the vector feature family is
``skyborn.plot.vector``:

.. code-block:: python

   from skyborn.plot.vector import curly_vector, curly_vector_key

The package-level names ``skyborn.plot.curly_vector`` and
``skyborn.plot.curly_vector_key`` remain available and forward to the same
public behavior.

The same ``skyborn.plot.curly_vector`` entry point supports both call styles:

- ``curly_vector(ax, x, y, u, v, ...)`` for direct NumPy / Matplotlib input
- ``curly_vector(ds, x="lon", y="lat", u="u", v="v", ax=ax, ...)`` for xarray dataset input

.. autofunction:: skyborn.plot.curly_vector_key

Scatter Stippling
-----------------

.. autofunction:: skyborn.plot.scatter

The same ``skyborn.plot.scatter`` entry point supports the standard
Matplotlib-style call forms:

- ``scatter(ax, x, y, ...)``
- ``scatter(x, y, ..., ax=ax)``
- ``scatter(x, y, ...)``

For gridded significance stippling, pass 1D grid axes together with a 2D
``where`` or ``mask`` field. Skyborn expands the grid candidates and then
applies NCL-style display-space thinning so dense lat-lon maps and
latitude-pressure sections do not need manual ``[::step]`` tuning. In most
cases, ``density`` is the only spacing control you need; ``distance`` is
available when you want to override the retained-point spacing explicitly.
For masked grids with inferable cell geometry, the default placement now fills
the selected cells with interior candidates before thinning so month-lat,
profile, and curvilinear-grid plots can place dots between coordinate centers.
Use ``placement="points"`` when you explicitly want the older node-based
behavior.

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import matplotlib.pyplot as plt
   import numpy as np
   import xarray as xr
   from skyborn.plot.vector import curly_vector

   # Create figure with equal axes
   fig, ax = plt.subplots()
   ax_new = skb.plot.add_equal_axes(ax, 'right', 0.1, 0.2)

   # Low-level array API
   lon = np.linspace(-180, 180, 50)
   lat = np.linspace(-90, 90, 25)
   u = np.random.random((25, 50))
   v = np.random.random((25, 50))

   curly_vector(ax, lon, lat, u, v)

   # Dataset wrapper API
   ds = xr.Dataset(
       {"u": (("lat", "lon"), u), "v": (("lat", "lon"), v)},
       coords={"lon": lon, "lat": lat},
   )
   curly_vector(ds, x="lon", y="lat", u="u", v="v", ax=ax)

   # Display-space-thinned stippling on a gridded significance mask
   p = xr.DataArray(
       np.random.random((25, 50)),
       dims=("lat", "lon"),
       coords={"lon": lon, "lat": lat},
   )
   skb.plot.scatter(ax, lon, lat, where=p < 0.05, density=2, s=4, c="0.3")

   # Vertical-profile stippling uses the same API
   level = np.linspace(1000, 100, 12)
   profile_sig = np.random.random((12, 25)) < 0.08
   fig2, ax2 = plt.subplots()
   skb.plot.scatter(ax2, lat, level, where=profile_sig, density=2, s=4, c="0.2")
   ax2.invert_yaxis()

Visualization Examples
----------------------

The plotting module is designed for creating publication-quality atmospheric visualizations:

.. code-block:: python

   import skyborn as skb
   import matplotlib.pyplot as plt
   from skyborn.plot.vector import curly_vector

   # Create figure
   fig = skb.plot.createFigure(figsize=(12, 8), dpi=300)
   ax = fig.add_subplot(111)

   # Plot wind vectors from arrays
   curly = curly_vector(ax, lon, lat, u_wind, v_wind)

   # Add significance stippling from a gridded mask
   stipple = skb.plot.scatter(ax, lon, lat, where=p_values < 0.05, density=2)

For complete visualization examples, see :doc:`../gallery`.
