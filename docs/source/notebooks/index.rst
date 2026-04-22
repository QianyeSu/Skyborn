Tutorials
=========

This page collects the tutorial notebooks shipped with Skyborn. Each card opens
the full notebook page with figures, code, and dataset notes, while the global
sidebar stays compact with a single ``Tutorials`` entry.

.. toctree::
   :hidden:
   :titlesonly:

   mann_kendall_tutorial
   gridfill_tutorial
   windspharm_tutorial
   ecs_emergent_constraints_analysis

.. raw:: html

   <div class="skyborn-tutorial-grid">
     <a class="skyborn-tutorial-card" href="mann_kendall_tutorial.html">
       <span class="skyborn-tutorial-card__media">
         <img class="skyborn-tutorial-card__cover" src="../images/precipitation_trends_comparison_1979_2014.png" alt="Global precipitation trend comparison figure used for the Mann-Kendall tutorial" loading="lazy">
       </span>
       <span class="skyborn-tutorial-card__body">
         <span class="skyborn-tutorial-card__eyebrow">Statistics</span>
         <h2>Mann-Kendall Tutorial</h2>
         <p>Real GPCP precipitation trends, four-method speed comparison, and multiple Mann-Kendall test families on global data.</p>
         <span class="skyborn-tutorial-card__meta">Real-data notebook</span>
       </span>
     </a>
     <a class="skyborn-tutorial-card" href="gridfill_tutorial.html">
       <span class="skyborn-tutorial-card__media">
         <img class="skyborn-tutorial-card__cover" src="../images/gridfill_comprehensive_comparison.png" alt="GridFill comparison figure used for the GridFill tutorial" loading="lazy">
       </span>
       <span class="skyborn-tutorial-card__body">
         <span class="skyborn-tutorial-card__eyebrow">Interpolation</span>
         <h2>GridFill Tutorial</h2>
         <p>Gap-filling workflows for atmospheric grids, including validation, method comparison, and publication-style figures.</p>
         <span class="skyborn-tutorial-card__meta">Workflow notebook</span>
       </span>
     </a>
     <a class="skyborn-tutorial-card" href="windspharm_tutorial.html">
       <span class="skyborn-tutorial-card__media">
         <img class="skyborn-tutorial-card__cover" src="../images/windspharm_sfvp_analysis.png" alt="Streamfunction and velocity potential figure used for the Windspharm tutorial" loading="lazy">
       </span>
       <span class="skyborn-tutorial-card__body">
         <span class="skyborn-tutorial-card__eyebrow">Dynamics</span>
         <h2>Windspharm Tutorial</h2>
         <p>Spherical-harmonic wind analysis with real data, vorticity, divergence, Helmholtz decomposition, and related diagnostics.</p>
         <span class="skyborn-tutorial-card__meta">Diagnostics notebook</span>
       </span>
     </a>
     <a class="skyborn-tutorial-card" href="ecs_emergent_constraints_analysis.html">
       <span class="skyborn-tutorial-card__media">
         <img class="skyborn-tutorial-card__cover" src="../images/ecs_emergent_constraints_analysis.png" alt="Emergent constraints dashboard used for the emergent constraints tutorial" loading="lazy">
       </span>
       <span class="skyborn-tutorial-card__body">
         <span class="skyborn-tutorial-card__eyebrow">Climate Constraints</span>
         <h2>Emergent Constraints Tutorial</h2>
         <p>Interactive climate-sensitivity analysis with observational constraints, uncertainty reduction, and reproducible figures.</p>
         <span class="skyborn-tutorial-card__meta">Analysis notebook</span>
       </span>
     </a>
   </div>

Tutorial Collection
-------------------

The tutorial hub is designed to keep the documentation entry points curated and
compact:

* The left sidebar keeps one ``Tutorials`` entry instead of growing one row per
  notebook.
* New notebooks can be added here as cards without making global navigation
  longer.
* Each tutorial still builds as a normal HTML documentation page and remains
  directly linkable.

Tutorial List
-------------

If you prefer a plain list of links:

* `Mann-Kendall Tutorial <mann_kendall_tutorial.html>`_
* `GridFill Tutorial <gridfill_tutorial.html>`_
* `Windspharm Tutorial <windspharm_tutorial.html>`_
* `Emergent Constraints Tutorial <ecs_emergent_constraints_analysis.html>`_

Run Tutorials Locally
---------------------

.. code-block:: bash

   pip install skyborn[dev]
   jupyter lab docs/source/notebooks/

Each notebook keeps its own data-loading notes near the top of the page.
