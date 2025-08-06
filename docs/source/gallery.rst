Gallery
=======

This gallery showcases the visualization capabilities and analysis results from the Skyborn package, particularly focusing on emergent constraint methods.

Interactive Notebooks
----------------------

Emergent Constraints Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a detailed, interactive analysis, see our comprehensive Jupyter notebook:

:doc:`notebooks/ecs_emergent_constraints_analysis`

This notebook demonstrates:

* Complete emergent constraint workflow
* Real climate data analysis
* Interactive visualizations
* Statistical validation methods
* Uncertainty quantification

Emergent Constraints Analysis Dashboard
---------------------------------------

Overview Dashboard
~~~~~~~~~~~~~~~~~~

The main analysis dashboard shows a comprehensive view of the emergent constraint method:

.. image:: images/emergent_constraints_dashboard.png
   :alt: Emergent Constraints Dashboard
   :width: 100%

*Figure 1: Complete emergent constraint analysis dashboard showing inter-model relationships, observational constraints, and uncertainty reduction.*

Key Components
~~~~~~~~~~~~~~

1. **Inter-model Relationship**

   * Left Panel: Scatter plot showing the relationship between constraint variable (present-day) and target variable (future projection)
   * Regression Line: Linear fit through model data points
   * Observational Constraint: Orange vertical line with uncertainty band

2. **Observational PDF**

   * Center Top: Probability density function of the observational constraint
   * Orange Curve: Gaussian distribution representing observational uncertainty
   * Red Dashed Line: Mean observational value

3. **Constraint Effect Comparison**

   * Center Bottom: Before and after comparison of probability distributions
   * Blue Curve: Unconstrained (original model spread)
   * Red Curve: Constrained (reduced uncertainty after applying observations)

4. **Uncertainty Reduction Statistics**

   * Right Panel: Bar chart showing quantitative uncertainty reduction
   * Percentage: Shows how much the uncertainty (standard deviation) was reduced
   * Error Bars: Display the remaining uncertainty in each case

Method Comparison
-----------------

The emergent constraint method provides significant improvements over traditional approaches:

.. image:: images/method_comparison.png
   :alt: Method Comparison
   :width: 80%

*Figure 2: Comparison of traditional vs. emergent constraint methods showing uncertainty reduction.*

ECS Analysis Results
--------------------

.. image:: images/ecs_emergent_constraints_analysis.png
   :alt: ECS Emergent Constraints Analysis
   :width: 100%

*Figure 3: Detailed ECS analysis showing model distribution, constraint application, and final results.*

Getting Started
---------------

To run these analyses yourself:

1. **Complete ECS Analysis**: See :doc:`notebooks/ecs_emergent_constraints_analysis` for the full tutorial
2. **Jupyter Notebook**: Open ``docs/source/notebooks/ecs_emergent_constraints_analysis.ipynb``
3. **Simple Demo**: Try ``examples/emergent_constraints_demo.ipynb`` for a quick start

Example Code
~~~~~~~~~~~~

.. code-block:: python

   import skyborn as skb
   import numpy as np

   # Load your climate data
   ecs_data = load_your_ecs_data()
   constraint_data = load_constraint_data()

   # Apply emergent constraint
   pdf = skb.gaussian_pdf(obs_mean, obs_std, x_grid)
   correlation = skb.pearson_correlation(constraint_data, ecs_data)

   # Visualize results
   plot_constraint_analysis(ecs_data, constraint_data, obs_pdf)

Technical Details
-----------------

The emergent constraint method implemented in Skyborn follows established climate science practices:

* **Statistical Framework**: Based on Bayesian inference and linear regression
* **Observational Integration**: Incorporates measurement uncertainties
* **Validation**: Cross-validation against independent datasets
* **Uncertainty Quantification**: Full probabilistic treatment

Spherical Harmonic Wind Analysis (Windspharm)
=============================================

Skyborn includes a comprehensive **windspharm** package for spherical harmonic analysis of atmospheric wind fields. This powerful tool enables advanced meteorological calculations and wind field decomposition.

Overview
--------

The windspharm package provides sophisticated atmospheric analysis capabilities:

.. image:: images/windspharm_vorticity_divergence.png
   :alt: Windspharm Vorticity and Divergence Analysis
   :width: 100%

*Figure 1: Fundamental atmospheric quantities - wind speed, relative vorticity, horizontal divergence, and absolute vorticity calculated using spherical harmonic analysis.*

.. note::
   **Key Features:**

   * **Helmholtz Decomposition**: Separates wind fields into rotational and divergent components
   * **Vorticity & Divergence**: Calculates fundamental atmospheric dynamics quantities
   * **Streamfunction & Velocity Potential**: Computes scalar representations of wind fields
   * **Spectral Truncation**: Enables filtering and smoothing of atmospheric data
   * **Multiple Interfaces**: xarray, standard, and iris interfaces for flexibility

Core Calculations
-----------------

**1. Fundamental Quantities**

The package calculates essential atmospheric dynamics fields:

* **Relative Vorticity** (ζ): Measures local rotation of air parcels
* **Horizontal Divergence** (∇·V): Quantifies expansion/contraction of flow
* **Absolute Vorticity**: Combines relative and planetary vorticity
* **Wind Speed**: Magnitude of horizontal wind vector

**2. Helmholtz Decomposition**

Separates any wind field into two fundamental components:

.. image:: images/windspharm_helmholtz.png
   :alt: Windspharm Helmholtz Decomposition
   :width: 100%

*Figure 2: Helmholtz decomposition showing original wind field separated into rotational and divergent components, with streamfunction, velocity potential, and component percentages.*

* **Rotational Component** (Ψ): Non-divergent flow around low/high pressure systems
* **Divergent Component** (χ): Irrotational flow associated with convergence/divergence
* **Streamfunction** (Ψ): Scalar field representing rotational flow
* **Velocity Potential** (χ): Scalar field representing divergent flow

**3. Advanced Analysis**

* **Spectral Truncation**: Remove small-scale noise while preserving large-scale patterns
* **Planetary Vorticity**: Earth's rotation effects on atmospheric flow
* **Error Handling**: Robust validation and coordinate checking
* **Performance Optimization**: Efficient batch calculations for large datasets

Streamfunction and Velocity Potential Analysis
-----------------------------------------------

.. image:: images/windspharm_streamfunction_potential.png
   :alt: Windspharm Streamfunction and Velocity Potential
   :width: 100%

*Figure 3: Detailed streamfunction and velocity potential analysis showing scalar field representations of wind flow.*

.. image:: images/windspharm_sfvp_analysis.png
   :alt: Windspharm SFVP Comparison Analysis
   :width: 100%

*Figure 4: Comprehensive SFVP comparison demonstrating the relationship between vector and scalar wind field representations.*

Component Analysis and Comparison
----------------------------------

.. image:: images/windspharm_component_comparison.png
   :alt: Windspharm Component Comparison
   :width: 100%

*Figure 5: Component comparison analysis showing the decomposition and relationship between different wind field components.*

Advanced Features
-----------------

**Gradient Analysis**

.. image:: images/windspharm_gradients.png
   :alt: Windspharm Gradient Analysis
   :width: 100%

*Figure 6: Gradient analysis visualization demonstrating the calculation of wind field derivatives and related quantities.*

**Spectral Truncation Effects**

.. image:: images/windspharm_truncation_comparison.png
   :alt: Windspharm Truncation Comparison
   :width: 100%

*Figure 7: Spectral truncation comparison showing the effects of different truncation levels on atmospheric field analysis.*

Mathematical Foundation
-----------------------

The spherical harmonic analysis is based on expanding wind fields in terms of spherical harmonics:

.. math::

   u(\lambda, \theta) = \sum_{n=0}^{N} \sum_{m=-n}^{n} u_n^m Y_n^m(\lambda, \theta)

Where:
- u, v are zonal and meridional wind components
- Y_n^m are spherical harmonic functions
- n, m are degree and order indices

Interactive Tutorial
--------------------

**Comprehensive Tutorial**: :doc:`notebooks/windspharm_tutorial`

The complete windspharm tutorial demonstrates:

1. **Data Loading**: Working with NetCDF atmospheric data
2. **Basic Calculations**: Vorticity, divergence, and wind speed
3. **Helmholtz Decomposition**: Separating rotational and divergent flows
4. **Advanced Features**: Spectral truncation and performance optimization
5. **Visualization**: Creating publication-quality atmospheric plots
6. **Best Practices**: Memory management and error handling

Example Applications
--------------------

**Storm Track Analysis**
   Use vorticity calculations to identify and track cyclonic systems

**Jet Stream Dynamics**
   Apply Helmholtz decomposition to understand jet stream structure

**Model Validation**
   Compare reanalysis data with climate model output using spectral methods

**Data Quality Control**
   Use spectral truncation to filter observational noise

Getting Started
---------------

.. code-block:: python

   from skyborn.windspharm.xarray import VectorWind
   import xarray as xr

   # Load your wind data
   ds = xr.open_dataset('Era5_Windfield_Data.nc')

   # Create VectorWind object
   vw = VectorWind(ds.u, ds.v)

   # Calculate fundamental quantities
   vorticity = vw.vorticity()
   divergence = vw.divergence()

   # Perform Helmholtz decomposition
   uchi, vchi, upsi, vpsi = vw.helmholtz()

   # Get streamfunction and velocity potential
   streamfunction = vw.streamfunction()
   velocity_potential = vw.velocitypotential()

Technical Notes
---------------

**Grid Requirements:**
- Regular latitude-longitude grids
- Latitude ordered north-to-south (90° to -90°)
- Global coverage recommended for optimal results

**Performance:**
- Use batch calculations (e.g., ``vrtdiv()``) for better performance
- Consider memory vs. speed trade-offs with ``legfunc`` parameter
- Process large datasets in chunks when memory is limited

**Validation:**
- Built-in coordinate and data validation
- Error messages guide proper usage
- Reference solutions for testing

References
----------

**Emergent Constraints:**
* **Methodology**: Cox, P. M., et al. (2013). Nature, 494(7437), 341-344
* **Implementation**: Based on https://github.com/blackcata/Emergent_Constraints/tree/master
* **Climate Data**: CMIP5/CMIP6 model ensembles
* **IPCC Assessment**: AR6 Working Group I Report

**Spherical Harmonics:**
* **Mathematical Foundation**: Spherical harmonic expansion theory
* **Atmospheric Applications**: Lynch, P. (2006). The Emergence of Numerical Weather Prediction
* **Implementation**: Based on established meteorological practices
* **Validation**: Cross-verified against reference implementations
