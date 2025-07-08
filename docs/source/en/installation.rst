Installation Guide
==================

Requirements
------------

Skyborn requires Python 3.8 or later. The library depends on several scientific Python packages:

Core Dependencies
~~~~~~~~~~~~~~~~~

* **NumPy** (>=1.19.0) - Numerical computing
* **Pandas** (>=1.3.0) - Data manipulation
* **Xarray** (>=0.19.0) - N-dimensional labeled arrays
* **SciPy** (>=1.7.0) - Scientific computing
* **Matplotlib** (>=3.4.0) - Plotting

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~

* **eccodes** (>=1.4.0) - For GRIB to NetCDF conversion
* **Cartopy** (>=0.20.0) - For geospatial plotting
* **Dask** (>=2021.6.0) - For parallel computing

Installation Methods
--------------------

From PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install skyborn

This will install Skyborn and its core dependencies.

From Conda
~~~~~~~~~~

.. code-block:: bash

   conda install -c conda-forge skyborn

From Source
~~~~~~~~~~~

For the latest development version:

.. code-block:: bash

   git clone https://github.com/yourusername/skyborn.git
   cd skyborn
   pip install -e .

Installing Optional Dependencies
--------------------------------

For GRIB Conversion
~~~~~~~~~~~~~~~~~~~

To use the GRIB to NetCDF conversion functionality:

.. code-block:: bash

   # Using conda (recommended)
   conda install -c conda-forge eccodes

   # Using pip
   pip install eccodes

For Geospatial Plotting
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Using conda
   conda install -c conda-forge cartopy

   # Using pip
   pip install cartopy

Complete Installation
~~~~~~~~~~~~~~~~~~~~~

To install all optional dependencies:

.. code-block:: bash

   # Using conda
   conda install -c conda-forge skyborn eccodes cartopy dask

   # Using pip
   pip install skyborn[complete]

Verification
------------

To verify your installation:

.. code-block:: python

   import skyborn
   print(f"Skyborn version: {skyborn.__version__}")

   # Test basic functionality
   skyborn.convert_longitude_range  # Should not raise an error

   # Test GRIB conversion (if eccodes is installed)
   try:
       skyborn.grib2nc
       print("GRIB conversion available ✓")
   except AttributeError:
       print("GRIB conversion not available (eccodes not installed)")

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

1. **Import Error**: Ensure all dependencies are installed
2. **GRIB Conversion Error**: Install eccodes library
3. **Plotting Issues**: Install matplotlib and optionally cartopy

Platform-Specific Notes
~~~~~~~~~~~~~~~~~~~~~~~

**Windows**:
  - Recommend using Anaconda/Miniconda
  - Some packages may require Visual Studio Build Tools

**macOS**:
  - May need to install Xcode command line tools
  - Use conda for easier dependency management

**Linux**:
  - Usually works out of the box
  - May need development headers for some packages

Development Installation
------------------------

For developers:

.. code-block:: bash

   git clone https://github.com/yourusername/skyborn.git
   cd skyborn
   pip install -e .[dev]

   # Install pre-commit hooks
   pre-commit install

This installs additional development dependencies including testing and documentation tools.
