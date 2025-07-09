Quick Start Guide
==================

This guide helps you get started with Skyborn quickly.

Installation
------------

Install Skyborn using pip:

.. code-block:: bash

   pip install skyborn

Basic Usage
-----------

Import and use Skyborn functions:

.. code-block:: python

   import skyborn as skb

   # Example: Gaussian PDF calculation
   pdf_values = skb.gaussian_pdf(mu=0, sigma=1, x=x_values)

   # Example: Correlation analysis
   correlation = skb.pearson_correlation(x_data, y_data)

Next Steps
----------

* Explore the :doc:`notebooks/index` for comprehensive examples
* Check the :doc:`gallery` for visualization examples
