from setuptools import setup, find_packages

setup(
    name="skyborn",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        'numpy',
        'xarray',
        'matplotlib',
        'cartopy',
        'netCDF4',
        'metpy',
        'cfgrib',
        'eccodes',
        'scikit-learn'
    ],
    python_requires='>=3.10',
    author="Qianye Su",
    author_email="suqianye2000@gmail.com",
    description="Atmospheric science research utilities",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/SQYQianYe/Skyborn",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
