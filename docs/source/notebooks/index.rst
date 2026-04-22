Tutorials
=========

Skyborn keeps tutorial notebooks as full documentation pages, but this page is
the only tutorial entry shown in the global sidebar. Individual notebooks are
linked from the tutorial hub below instead of being expanded into a long
left-hand navigation tree.

.. raw:: html

   <div class="skyborn-tutorial-grid">
     <a class="skyborn-tutorial-card" href="mann_kendall_tutorial.html">
       <span class="skyborn-tutorial-card__eyebrow">Statistics</span>
       <h2>Mann-Kendall Tutorial</h2>
       <p>Real GPCP precipitation trends, four-method speed comparison, and multiple Mann-Kendall test families on global data.</p>
       <span class="skyborn-tutorial-card__meta">Notebook tutorial</span>
     </a>
     <a class="skyborn-tutorial-card" href="gridfill_tutorial.html">
       <span class="skyborn-tutorial-card__eyebrow">Interpolation</span>
       <h2>GridFill Tutorial</h2>
       <p>Gap-filling workflows for atmospheric grids, including validation, method comparison, and publication-style figures.</p>
       <span class="skyborn-tutorial-card__meta">Notebook tutorial</span>
     </a>
     <a class="skyborn-tutorial-card" href="windspharm_tutorial.html">
       <span class="skyborn-tutorial-card__eyebrow">Dynamics</span>
       <h2>Windspharm Tutorial</h2>
       <p>Spherical-harmonic wind analysis with real data, vorticity, divergence, Helmholtz decomposition, and related diagnostics.</p>
       <span class="skyborn-tutorial-card__meta">Notebook tutorial</span>
     </a>
     <a class="skyborn-tutorial-card" href="ecs_emergent_constraints_analysis.html">
       <span class="skyborn-tutorial-card__eyebrow">Climate Constraints</span>
       <h2>Emergent Constraints Tutorial</h2>
       <p>Interactive climate-sensitivity analysis with observational constraints, uncertainty reduction, and reproducible figures.</p>
       <span class="skyborn-tutorial-card__meta">Notebook tutorial</span>
     </a>
   </div>

Why This Layout
---------------

This tutorial hub keeps the global navigation compact:

* The sidebar shows a single ``Tutorials`` entry instead of one button per
  notebook.
* New notebooks can be added here as cards without making the left navigation
  tree longer.
* Notebook pages still build as normal HTML pages and remain directly linkable.

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
