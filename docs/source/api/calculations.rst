Calculation
=====================

The Skyborn calculation module provides statistical and mathematical functions for climate data analysis.

Statistical Functions
---------------------

.. autofunction:: skyborn.calc.linear_regression
   :no-index:

.. autofunction:: skyborn.calc.pearson_correlation
   :no-index:

.. autofunction:: skyborn.calc.spearman_correlation
   :no-index:

.. autofunction:: skyborn.calc.kendall_correlation
   :no-index:

Emergent Constraint Functions
-----------------------------

.. autofunction:: skyborn.calc.gaussian_pdf
   :no-index:

.. autofunction:: skyborn.calc.emergent_constraint_posterior
   :no-index:

.. autofunction:: skyborn.calc.emergent_constraint_prior
   :no-index:

Legacy Functions (for backward compatibility)
---------------------------------------------

.. autofunction:: skyborn.calc.calc_GAUSSIAN_PDF
   :no-index:

.. autofunction:: skyborn.calc.calc_PDF_EC
   :no-index:

.. autofunction:: skyborn.calc.find_std_from_PDF
   :no-index:

.. autofunction:: skyborn.calc.calc_PDF_EC_PRIOR
   :no-index:

Utility Functions
-----------------

.. autofunction:: skyborn.calc.calculate_potential_temperature
   :no-index:

.. autofunction:: skyborn.calc.convert_longitude_range
   :no-index:

Mann-Kendall Trend Analysis
---------------------------

.. autofunction:: skyborn.calc.mann_kendall_test
   :no-index:

.. autofunction:: skyborn.calc.mann_kendall_multidim
   :no-index:

.. autofunction:: skyborn.calc.mann_kendall_xarray
   :no-index:

.. autofunction:: skyborn.calc.trend_analysis
   :no-index:

.. autofunction:: skyborn.calc.mk_test
   :no-index:

.. autofunction:: skyborn.calc.mk_multidim
   :no-index:

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import numpy as np

   # Statistical analysis
   x = np.random.randn(100)
   y = 2 * x + np.random.randn(100) * 0.5

   slope, intercept, r_value = skb.linear_regression(x, y)
   correlation = skb.pearson_correlation(x, y)

   # Emergent constraint analysis
   x_values = np.linspace(-3, 3, 100)
   pdf = skb.gaussian_pdf(mu=0, sigma=1, x=x_values)

   # Apply emergent constraint
   posterior_mean, posterior_std = skb.emergent_constraint_posterior(
       prior_mean=3.0, prior_std=1.5,
       obs_mean=0.5, obs_std=0.2,
       relationship_slope=2.0, relationship_intercept=0.1
   )

Emergent Constraints
====================

The emergent constraint module implements methods for reducing uncertainty in climate projections
by leveraging relationships between observable present-day quantities and uncertain future projections.

For a complete interactive example, see :doc:`../notebooks/ecs_emergent_constraints_analysis`.
