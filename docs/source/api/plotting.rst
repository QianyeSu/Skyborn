Plotting
==================

The Skyborn plotting module provides specialized visualization functions for atmospheric and climate data.

Plotting Utilities
------------------

.. autofunction:: skyborn.plot.add_equal_axes
   :no-index:

.. autofunction:: skyborn.plot.createFigure
   :no-index:

Curved Quiver Plots
-------------------

.. autofunction:: skyborn.plot.curved_quiver
   :no-index:

Modular Plotting
----------------

.. automodule:: skyborn.plot.modplot
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import matplotlib.pyplot as plt
   import numpy as np

   # Create figure with equal axes
   fig, ax = plt.subplots()
   ax_new = skb.plot.add_equal_axes(ax, 'right', 0.1, 0.2)

   # Create curved quiver plot
   lon = np.linspace(-180, 180, 50)
   lat = np.linspace(-90, 90, 25)
   u = np.random.random((25, 50))
   v = np.random.random((25, 50))

   skb.plot.curved_quiver(ax, lon, lat, u, v)

Visualization Examples
----------------------

The plotting module is designed for creating publication-quality atmospheric visualizations:

.. code-block:: python

   import skyborn as skb
   import matplotlib.pyplot as plt

   # Create figure
   fig, ax = skb.plot.createFigure((12, 8), 1, 1)

   # Plot wind vectors
   quiver = skb.plot.curved_quiver(ax, lon, lat, u_wind, v_wind)

For complete visualization examples, see :doc:`../gallery`.
