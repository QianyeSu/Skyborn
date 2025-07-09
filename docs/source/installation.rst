Installation
============

System Requirements
-------------------

* Python 3.8 or higher
* Windows, macOS, or Linux

Installation Methods
--------------------

From PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install skyborn

Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~

For development or to get the latest features:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/your-username/skyborn.git
   cd skyborn

   # Install in development mode
   pip install -e .

From Git Repository
~~~~~~~~~~~~~~~~~~~

Install directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/your-username/skyborn.git

With Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~

Install with all optional dependencies for full functionality:

.. code-block:: bash

   pip install skyborn[all]

Or install specific dependency groups:

.. code-block:: bash

   # For visualization
   pip install skyborn[plotting]

   # For development
   pip install skyborn[dev]

   # For documentation building
   pip install skyborn[docs]

Conda Installation
~~~~~~~~~~~~~~~~~~

Skyborn can also be installed using conda:

.. code-block:: bash

   conda install -c conda-forge skyborn

Verify Installation
-------------------

Test your installation:

.. code-block:: python

   import skyborn as skb
   print(skb.__version__)

   # Test a basic function
   import numpy as np
   x = np.linspace(-3, 3, 100)
   pdf = skb.gaussian_pdf(0, 1, x)
   print("Installation successful!")

Dependencies
------------

Core dependencies:
* numpy
* pandas
* xarray
* matplotlib
* scipy

Optional dependencies:
* seaborn (for enhanced plotting)
* cartopy (for geographic plotting)
* eccodes (for GRIB file handling)
* netCDF4 (for NetCDF file operations)

Troubleshooting
---------------

If you encounter issues:

1. **Import errors**: Make sure all dependencies are installed
2. **GRIB support**: Install eccodes for GRIB file handling
3. **Plotting issues**: Install cartopy and seaborn for full plotting functionality

For support, please visit our GitHub repository.
