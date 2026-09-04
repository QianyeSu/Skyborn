Causality
==================

.. automodule:: skyborn.calc.causality
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: granger_causality, liang_causality

Overview
--------

The ``skyborn.calc.causality`` module provides methods for analyzing causal relationships between time series in atmospheric and climate data. This module implements both Granger causality and Liang-Kleeman information flow methods with comprehensive significance testing. ``skyborn.causality`` remains available as a compatibility import.

Key Features
------------

- **Granger Causality**: Statistical test for causality based on prediction improvement
- **Liang Information Flow**: Physically-based causality measure using information theory
- **Significance Testing**: Multiple methods for statistical significance assessment
- **AR(1) Modeling**: Autoregressive model fitting for red noise generation
- **Phase Randomization**: Surrogate data generation for null hypothesis testing

Implementation and Build
------------------------

The maintained Python implementation is in
``src/skyborn/calc/causality/core.py``. The public package exports are defined
in ``src/skyborn/calc/causality/__init__.py``. The historical
``skyborn.causality`` module remains a thin compatibility entry point.

Liang information flow and the batch surrogate calculations use a direct
C/Fortran extension built by ``src/skyborn/calc/causality/meson.build``:

.. code-block:: text

   src/skyborn/calc/causality/
   |-- __init__.py
   |-- core.py
   |-- causality_core_module.c
   |-- fortran/
   |   |-- causality_core.f90
   |   `-- causality_core_bindings.f90
   `-- meson.build

There is deliberately no Python ``bindings.py`` file. The C extension is the
thin NumPy boundary, and the Fortran bindings are compiled as part of the
Meson target.

Surrogate generation is also batched at the Python/NumPy boundary:

* ``phaseran`` uses batched real FFTs and preserves the historical random
  phase stream.
* ``sm_ar1_sim`` preserves the historical random stream and can use a
  Fortran/OpenMP AR(1) filtering kernel for larger surrogate ensembles.
* Liang significance calculations pass the generated surrogate matrices to
  the Fortran/OpenMP batch kernel.

The :class:`~skyborn.calc.causality.Surrogate` facade provides a common
``generate`` entry point for ``isospec`` phase-randomized and ``isopersist``
AR(1) null models. The AR(1) native path accelerates filtering after NumPy has
generated the innovations; this keeps seeded output compatible with the
historical statsmodels implementation.

On Windows, build and test this module from the ``skyborn_dev`` Anaconda
environment:

.. code-block:: powershell

   $env:CC = "F:\Anaconda3\envs\skyborn_dev\Library\bin\gcc.exe"
   $env:FC = "F:\Anaconda3\envs\skyborn_dev\Library\bin\gfortran.exe"
   $env:PATH = "F:\Anaconda3\envs\skyborn_dev\Library\bin;$env:PATH"
   $env:PYTHONPATH = "D:\Skyborn\src"
   F:\Anaconda3\envs\skyborn_dev\python.exe -m mesonbuild.mesonmain compile -C src\skyborn\calc\causality\build_gfortran
   F:\Anaconda3\envs\skyborn_dev\python.exe -m pytest -q -o addopts= tests\test_causality.py

Methods Available
-----------------

Causality Analysis
~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.causality.granger_causality

.. autofunction:: skyborn.calc.causality.liang_causality

Significance Testing
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.causality.signif_isopersist

.. autofunction:: skyborn.calc.causality.signif_isospec

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.causality.ar1_fit_evenly

.. autoclass:: skyborn.calc.causality.Surrogate
   :members:

.. autofunction:: skyborn.calc.causality.phaseran

Theoretical Background
----------------------

Granger Causality
~~~~~~~~~~~~~~~~~

Granger causality tests whether past values of one time series help predict another time series beyond what can be predicted from the target series alone. The null hypothesis states that the second series does NOT Granger-cause the first.

**Mathematical Foundation:**

For time series X and Y, Y Granger-causes X if:

.. math::

   \sigma^2(X_t | X_{t-1}, X_{t-2}, \ldots) > \sigma^2(X_t | X_{t-1}, X_{t-2}, \ldots, Y_{t-1}, Y_{t-2}, \ldots)

where :math:`\sigma^2` denotes the prediction error variance.

Liang-Kleeman Information Flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Liang method quantifies information flow between time series using rigorous information theory principles. It measures the rate of information transfer from one series to another.

**Key Metrics:**

- **T21**: Information flow from series 2 to series 1
- **tau21**: Normalized information flow (relative to total information)
- **Z**: Total information in the system

.. math::

   T_{2 \rightarrow 1} = \frac{C_{12}}{C_{11}} \cdot \frac{-C_{21}\frac{dC_{11}}{dt} + C_{11}\frac{dC_{21}}{dt}}{|C|}

where :math:`C_{ij}` are covariance matrix elements.

Usage Examples
--------------

Basic Granger Causality Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import skyborn.calc.causality as scaus
   import numpy as np

   # Generate sample atmospheric time series
   np.random.seed(42)
   n_samples = 1000

   # Temperature-like series
   temp = np.cumsum(np.random.randn(n_samples)) * 0.1

   # Pressure-like series with some dependence on temperature
   pressure = np.zeros(n_samples)
   for i in range(1, n_samples):
       pressure[i] = 0.7 * pressure[i-1] + 0.3 * temp[i-1] + np.random.randn()

   # Test if temperature Granger-causes pressure
   gc_result = scaus.granger_causality(pressure, temp, maxlag=5)

   # Extract p-values for different lags
   for lag in gc_result:
       f_stat = gc_result[lag][0]['ssr_ftest'][0]
       p_value = gc_result[lag][0]['ssr_ftest'][1]
       print(f"Lag {lag}: F-statistic = {f_stat:.3f}, p-value = {p_value:.3f}")

Liang Information Flow Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import skyborn.calc.causality as scaus
   import numpy as np

   # Generate coupled time series
   np.random.seed(123)
   n = 500

   # Series 1: autonomous dynamics
   x1 = np.zeros(n)
   for i in range(1, n):
       x1[i] = 0.8 * x1[i-1] + np.random.randn()

   # Series 2: driven by series 1
   x2 = np.zeros(n)
   for i in range(1, n):
       x2[i] = 0.5 * x2[i-1] + 0.4 * x1[i-1] + np.random.randn()

   # Calculate Liang causality with significance testing
   result = scaus.liang_causality(x2, x1, signif_test='isospec', nsim=1000)

   print(f"Information flow (T21): {result['T21']:.4f}")
   print(f"Normalized flow (tau21): {result['tau21']:.4f}")
   print(f"Total information (Z): {result['Z']:.4f}")

   # Check significance
   sig_level = 0.05
   lower_bound = result['T21_noise'][1]  # 2.5th percentile
   upper_bound = result['T21_noise'][-2]  # 97.5th percentile

   if result['T21'] > upper_bound or result['T21'] < lower_bound:
       print(f"Causality is significant at {(1-sig_level)*100}% level")
   else:
       print("Causality is not significant")

Atmospheric Science Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import skyborn.calc.causality as scaus
   import numpy as np
   import matplotlib.pyplot as plt

   # Simulate ENSO-like and temperature-like indices
   def generate_enso_temp_data(n_years=50):
       n_months = n_years * 12
       t = np.arange(n_months)

       # ENSO-like oscillation (irregular ~3-7 year cycle)
       enso_base = np.sin(2 * np.pi * t / 42) + 0.5 * np.sin(2 * np.pi * t / 84)
       enso_noise = np.random.randn(n_months) * 0.5
       enso = enso_base + enso_noise

       # Temperature anomaly influenced by ENSO with lag
       temp = np.zeros(n_months)
       for i in range(3, n_months):
           temp[i] = 0.6 * temp[i-1] + 0.3 * enso[i-3] + np.random.randn() * 0.3

       return enso, temp

   # Generate data
   enso_index, temp_anomaly = generate_enso_temp_data(40)

   # Test causality: Does ENSO cause temperature changes?
   liang_result = scaus.liang_causality(
       temp_anomaly, enso_index,
       signif_test='isopersist',
       nsim=2000
   )

   print("ENSO → Temperature Analysis:")
   print(f"Information Flow: {liang_result['T21']:.4f}")
   print(f"Normalized Flow: {liang_result['tau21']:.4f}")

   # Compare with Granger causality
   granger_result = scaus.granger_causality(temp_anomaly, enso_index, maxlag=6)

   print("\\nGranger Causality Results:")
   for lag in [1, 3, 6]:
       if lag in granger_result:
           p_val = granger_result[lag][0]['ssr_ftest'][1]
           print(f"Lag {lag}: p-value = {p_val:.4f}")

Significance Testing Methods
----------------------------

Isopersistent Testing (``signif_isopersist``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tests significance using AR(1) surrogates with the same persistence (autocorrelation) as the original data. This method preserves the red noise characteristics while removing any causal relationships.

**When to use:**
- When your data shows significant autocorrelation
- For testing against red noise null hypothesis
- When computational efficiency is important

Isospectral Testing (``signif_isospec``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tests significance using phase-randomized surrogates that preserve the power spectrum of the original data. This method maintains spectral properties while destroying phase relationships.

**When to use:**
- When preserving spectral characteristics is important
- For more conservative significance testing
- When dealing with complex periodic behaviors

Interpretation Guidelines
-------------------------

Granger Causality Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **p-value < 0.05**: Reject null hypothesis; evidence for Granger causality
- **p-value ≥ 0.05**: Fail to reject null hypothesis; no evidence for causality
- Consider multiple lag values to capture different timescale relationships
- Be aware that Granger causality tests statistical precedence, not physical causation

Liang Information Flow Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **T21 > 0**: Positive information flow from series 2 to series 1
- **T21 < 0**: Negative information flow (information destruction)
- **|tau21| → 1**: Strong relative causality
- **|tau21| → 0**: Weak relative causality
- Compare against significance bounds from surrogate testing

Best Practices
--------------

Data Preparation
~~~~~~~~~~~~~~~~

1. **Stationarity**: Ensure time series are stationary or properly detrended
2. **Length**: Use sufficiently long time series (typically > 100 points)
3. **Sampling**: Ensure consistent and appropriate sampling rates
4. **Missing Data**: Handle gaps appropriately before analysis

Statistical Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Multiple Testing**: Apply correction for multiple hypothesis testing
2. **Lag Selection**: Test multiple lag values for Granger causality
3. **Surrogate Count**: Use adequate number of surrogates (≥ 1000) for significance testing
4. **Cross-Validation**: Validate results on independent data when possible

Physical Interpretation
~~~~~~~~~~~~~~~~~~~~~~~

1. **Mechanism**: Consider physical mechanisms that could explain causal relationships
2. **Timescales**: Match analysis timescales to relevant physical processes
3. **Confounding**: Be aware of potential confounding variables
4. **Bidirectionality**: Test causality in both directions

Related Modules
---------------

- :mod:`skyborn.calc`: Statistical calculations and emergent constraint analysis
- :mod:`skyborn.interp`: Time series interpolation and regridding
- :mod:`skyborn.plot`: Visualization tools for causality results

References
----------

**Granger Causality:**

1. Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438.

2. Granger, C. W. J. (1980). Testing for causality: A personal viewpoint. *Journal of Economic Dynamics and Control*, 2, 329-352.

**Liang-Kleeman Information Flow:**

3. Liang, X. S. (2013). The Liang-Kleeman Information Flow: Theory and Applications. *Entropy*, 15, 327-360.

4. Liang, X. S. (2014). Unraveling the cause-effect relation between time series. *Physical Review E*, 90, 052150.

5. Liang, X. S. (2016). Information flow and causality as rigorous notions ab initio. *Physical Review E*, 94, 052201.

**Surrogate Methods:**

6. Prichard, D., & Theiler, J. (1994). Generating surrogate data for time series with several simultaneously measured variables. *Physical Review Letters*, 73(7), 951-954.
