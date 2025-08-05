.. _notebooks:

Jupyter Notebooks
=================

This section contains interactive Jupyter notebooks demonstrating various Skyborn capabilities.

.. toctree::
   :maxdepth: 2
   :caption: Available Notebooks:

   ecs_emergent_constraints_analysis.ipynb
   windspharm_tutorial.ipynb

Emergent Constraints Analysis
-----------------------------

The :doc:`ecs_emergent_constraints_analysis` notebook provides a comprehensive
demonstration of emergent constraint methods for climate sensitivity analysis.

**Key Features:**

* Interactive data visualization
* Step-by-step constraint methodology
* Real climate model data analysis
* Uncertainty quantification
* Statistical validation

**Topics Covered:**

* Loading and preprocessing climate data
* Calculating inter-model relationships
* Applying observational constraints
* Quantifying uncertainty reduction
* Visualizing results

Vector Wind Analysis
--------------------

The :doc:`windspharm_tutorial` notebook demonstrates spherical harmonic vector wind
analysis using the Skyborn windspharm module.

**Key Features:**

* Spherical harmonic transforms
* Vorticity and divergence calculations
* Helmholtz decomposition
* Streamfunction and velocity potential
* Time series analysis

**Topics Covered:**

* Creating realistic wind data
* Computing dynamical quantities
* Decomposing wind into components
* Visualizing results
* Performance optimization
* Error handling

Getting Started
---------------

To run these notebooks locally:

1. Install Skyborn with development dependencies
2. Launch Jupyter Lab or Notebook
3. Navigate to the notebooks directory
4. Execute cells step by step

.. code-block:: bash

   pip install skyborn[dev]
   jupyter lab docs/source/notebooks/

Note
----

All notebooks include fallback synthetic data generation in case the required
data files are not available, ensuring they can be executed in any environment.
