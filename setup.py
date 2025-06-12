import setuptools

setuptools.setup(
    name="skyborn",
    version="0.3.5",  # Keep consistent with pyproject.toml
    author="Qianye Su",
    author_email="suqianye2000@gmail.com",
    description="Atmospheric science research utilities",
    long_description="Skyborn is a tool for easy plotting ERA5 weather data.",
    license="MIT",
    license_files=("LICENSE.txt"),
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

    install_requires=[
        "numpy",
        "xarray",
        "matplotlib",
        "cartopy",
        "netCDF4",
        "metpy",
        "tqdm",
        "statsmodels",
        "scipy",
        "cfgrib",
        "eccodes",
        "scikit-learn"
    ]
)
