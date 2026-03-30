Installation
============

System Requirements
-------------------

* Python 3.9 or higher
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
   git clone https://github.com/QianyeSu/skyborn.git
   cd skyborn

   # Install in development mode
   pip install -e .

From Git Repository
~~~~~~~~~~~~~~~~~~~

Install directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/QianyeSu/skyborn.git

With Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~

Install specific dependency groups:

.. code-block:: bash

   # Optional CDO bridge on supported wheel platforms
   pip install "skyborn[cdo]"

   # Optional GRIB backends
   pip install "skyborn[grib]"

   # For development
   pip install "skyborn[dev]"

   # For documentation building
   pip install "skyborn[docs]"

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
* netCDF4
* metpy
* scikit-learn
* dask
* statsmodels
* tqdm

Optional dependencies:
* seaborn (for enhanced plotting)
* cartopy (for geographic plotting)
* cfgrib and eccodes (for GRIB file handling)
* skyborn-cdo (via ``skyborn[cdo]`` on supported platforms)

Troubleshooting
---------------

If you encounter issues:

1. **Import errors**: Make sure all dependencies are installed
2. **GRIB support**: Install ``skyborn[grib]`` or install ``cfgrib`` and ``eccodes`` manually
3. **Plotting issues**: Install cartopy and seaborn for full plotting functionality

For support, please visit our GitHub repository.
