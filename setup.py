import setuptools

setuptools.setup(
    name="skyborn",
    version="0.1.0",  # 与pyproject.toml保持一致
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
        "numpy>=1.20.0",
        "xarray>=0.19.0",
        "matplotlib>=3.4.0",
        "cartopy>=0.20.0",
        "netCDF4>=1.5.7",
        "metpy>=1.1.0",
        "cfgrib>=0.9.9",
        "eccodes>=1.4.0",
        "scikit-learn>=1.0.0"
    ]
)
