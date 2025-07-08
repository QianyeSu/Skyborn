Contributing
============

We welcome contributions to Skyborn!

Getting Started
---------------

1. Fork the repository
2. Create a new branch for your feature
3. Make your changes
4. Run tests
5. Submit a pull request

Development Setup
-----------------

.. code-block:: bash

   git clone https://github.com/yourusername/skyborn.git
   cd skyborn
   pip install -e .[dev]
   pre-commit install

Code Style
----------

We use:

- Black for code formatting
- Pre-commit hooks for quality checks
- Pytest for testing

Testing
-------

Run tests with:

.. code-block:: bash

   pytest tests/
