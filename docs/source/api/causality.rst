Causality Functions
===================

Causality analysis methods for climate data.

Granger Causality
-----------------

.. autofunction:: skyborn.causality.granger_causality
   :no-index:

Liang Causality
---------------

.. autofunction:: skyborn.causality.liang_causality
   :no-index:

.. autofunction:: skyborn.causality.liang

Significance Testing
--------------------

.. autofunction:: skyborn.causality.signif_isopersist

.. autofunction:: skyborn.causality.signif_isospec

Utility Functions
-----------------

.. autofunction:: skyborn.causality.ar1_fit_evenly

.. autofunction:: skyborn.causality.phaseran

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import numpy as np

   # Generate sample time series
   x = np.random.randn(1000)
   y = np.random.randn(1000)

   # Test Granger causality
   gc_result = skb.granger_causality(x, y, maxlag=10)

   # Calculate Liang causality
   liang_result = skb.liang_causality(x, y)

   print(f"Granger causality results: {gc_result}")
   print(f"Liang causality: {liang_result}")

Example Usage
-------------

.. code-block:: python

   import skyborn as skb
   import numpy as np

   # Generate sample time series
   x = np.random.randn(1000)
   y = np.random.randn(1000)

   # Test Granger causality
   gc_result = skb.granger_causality_test(x, y, max_lag=10)

   # Calculate Liang causality
   liang_result = skb.liang_causality(x, y)

   print(f"Granger causality p-value: {gc_result['p_value']}")
   print(f"Liang causality: {liang_result}")
