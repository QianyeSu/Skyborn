<p align="center">
  <a href="https://github.com/QianyeSu/Skyborn" target="_blank">
    <img src="https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/_static/SkyBornLogo.svg" alt="Skyborn Logo" width="400"/>
  </a>
</p>

[![PyPI version](https://badge.fury.io/py/skyborn.svg)](https://badge.fury.io/py/skyborn)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/skyborn)](https://pypi.org/project/skyborn/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/skyborn)](https://pypi.org/project/skyborn/)
[![codecov](https://codecov.io/gh/QianyeSu/Skyborn/graph/badge.svg?token=YOUR_TOKEN_HERE)](https://codecov.io/gh/QianyeSu/Skyborn)
[![License](https://img.shields.io/github/license/QianyeSu/Skyborn)](https://github.com/QianyeSu/Skyborn/blob/main/LICENSE)
[![Tests](https://github.com/QianyeSu/Skyborn/actions/workflows/stable-ci.yml/badge.svg)](https://github.com/QianyeSu/Skyborn/actions/workflows/stable-ci.yml)
[![Platform](https://img.shields.io/badge/platform-Windows-blue)](https://github.com/QianyeSu/Skyborn)
[![Code style](https://img.shields.io/badge/code%20style-black-blue.svg)](https://github.com/psf/black)
[![Build Status](https://github.com/QianyeSu/Skyborn/actions/workflows/test-coverage.yml/badge.svg?branch=main)](https://github.com/QianyeSu/Skyborn/actions/workflows/test-coverage.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-brightgreen)](https://skyborn.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18075252.svg)](https://doi.org/10.5281/zenodo.18075252)

## System Requirements

**Operating System:** 🖥️ **Cross-Platform**

This package supports Windows, Linux, and macOS. However, it has been primarily developed and tested on Windows.

**Note:** While the package can be installed on different platforms, some Windows-specific features may not work on other operating systems.

## Installation

To install the Skyborn package, you can use pip:

```bash
pip install skyborn
```
or

```bash
pip install -U --index-url https://pypi.org/simple/ skyborn
```

## 📚 Documentation

**Full documentation is available at: [Documentation ](https://skyborn.readthedocs.io/en/latest/)**



## 🎯 Key Features & Submodules

### 📊 Spatial Trend Analysis & Climate Index Regression

Skyborn provides ultra-fast spatial trend calculation and climate index regression analysis for atmospheric data:

![Precipitation Trends Comparison](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/precipitation_trends_comparison_1979_2014.png)

**Key Capabilities:**
- **High-Speed Spatial Trends**: Calculate long-term climate trends across global grids
  - Linear trend analysis for temperature, precipitation, and other variables
  - Statistical significance testing
  - Vectorized operations for massive datasets

- **Climate Index Regression**: Rapid correlation and regression analysis with climate indices
  - NINO 3.4, PDO, NAO, AMO index integration
  - Pattern correlation analysis
  - Teleconnection mapping

**Other Applications:**
- Climate change signal detection
- Decadal variability analysis
- Teleconnection pattern identification
- Regional climate impact assessment

### 🌍 Skyborn Windspharm Submodule - Atmospheric Analysis

The Skyborn `windspharm` submodule provides powerful tools for analyzing global wind patterns through **streamfunction** and **velocity potential** calculations:

![Streamfunction and Velocity Potential](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/windspharm_sfvp_analysis.png)

**Key Capabilities:**
- **Streamfunction Analysis**: Identifies rotational (non-divergent) wind components
  - Visualizes atmospheric circulation patterns
  - Reveals jet streams and vortices
  - Essential for understanding weather systems

- **Velocity Potential Analysis**: Captures divergent wind components
  - Shows areas of convergence and divergence
  - Critical for tropical meteorology
  - Identifies monsoon circulation patterns

**Applications:**
- Climate dynamics research
- Weather pattern analysis
- Atmospheric wave propagation studies
- Tropical cyclone formation analysis

### 🔧 Skyborn Gridfill Submodule - Data Interpolation

The Skyborn `gridfill` submodule provides advanced interpolation techniques for filling missing data in atmospheric and climate datasets:

![Gridfill Missing Data Interpolation](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/gridfill_demo_result_readme.png)

**Key Features:**
- **Poisson-based Interpolation**: Physically consistent gap filling
- **Preserves Data Patterns**: Maintains spatial correlations and gradients
- **Multiple Methods Available**:
  - Basic Poisson solver
  - High-precision iterative refinement
  - Zonal initialization options
  - Relaxation parameter tuning

**Applications:**
- Satellite data gap filling
- Model output post-processing
- Climate data reanalysis
- Quality control for observational datasets

The example above demonstrates filling gaps in global precipitation data, where the algorithm successfully reconstructs missing values while preserving the underlying meteorological patterns.

### Curly Vector Plotting

Skyborn provides an NCL-like `curly_vector` renderer for atmospheric wind fields while keeping a Matplotlib and Cartopy compatible API.

![Curly Vector WRF Example](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/wrf_curly_vector_integration_direction_density08_ref10.png)

**Key Capabilities:**
- **NCL-like Curved Vectors**: Curved tails and lighter arrow heads for cleaner large-scale flow visualization
- **Matplotlib-Compatible API**: Works with standard axes and Cartopy projections without switching to a separate plotting stack
- **Flexible Inputs**: Supports regular lat-lon grids, xarray datasets, and curvilinear grids such as WRF mass-grid wind fields
- **Vertical Cross-Section Support**: Supports latitude-pressure and similar vertical profile plots through the same `curly_vector` API
- **Projection-Aware Rendering**: Supports global rendering with different `central_longitude` settings and Cartopy projections

Skyborn `curly_vector` also supports vertical cross-section rendering, for example latitude-pressure plots using meridional wind and scaled omega:

![Curly Vector Vertical Profile Example](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/hadley_msf_v_omega_ncl_curly.png)

**Minimal Example:**
```python
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from skyborn.plot import curly_vector

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))

curly_vector(
    ds,
    x="longitude",
    y="latitude",
    u="u",
    v="v",
    ax=ax,
    transform=ccrs.PlateCarree(),
    density=0.9,
    color="black",
    arrowstyle="->",
    linewidth=1.2,
    arrowsize=1.2,
    zorder=3,
    integration_direction='both',
    ref_magnitude=30.0,
    ref_length=0.1,
)
ax.coastlines()
ax.set_global()
```

## Performance Benchmarks

### 🚀 Windspharm Performance

The Skyborn `windspharm` submodule delivers **~25% performance improvement** over standard implementations through modernized Fortran code and optimized algorithms:

![Windspharm Performance Comparison](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/windspharm_performance_comparison.png)

**Key Performance Metrics:**
- **Vorticity Calculation**: ~25% faster
- **Divergence Calculation**: ~25% faster
- **Helmholtz Decomposition**: ~25% faster
- **Streamfunction/Velocity Potential**: ~25% faster

### ⚡ GPI Module Performance

The Genesis Potential Index (GPI) module achieves **dramatic speedups** through vectorized Fortran implementation and native 3D processing:

![GPI Speed Comparison](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/gpi_speed_comparison.png)

**Performance Highlights:**
- **19-25x faster** than point-by-point implementations
- Processes entire atmospheric grids in seconds
- Native multi-dimensional support (3D/4D data)

![GPI Global Distribution](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/gpi_global_distribution.png)

**Accuracy Validation:**
- Correlation coefficient > 0.99 with reference implementations
- RMSE < 1% for both VMAX and PMIN calculations

![GPI Scatter Comparison](https://raw.githubusercontent.com/QianyeSu/Skyborn/main/docs/source/images/gpi_scatter_comparison.png)

## 📖 Citation

If you use Skyborn in your research, please cite it using the following format:

```bibtex
@software{su2025skyborn,
  author = {Su, Qianye},
  title = {Skyborn: Climate Data Analysis Toolkit},
  year = {2025},
  doi = {10.5281/zenodo.18075252},
  url = {https://doi.org/10.5281/zenodo.18075252}
}
```

Or in text:
> Su, Q. (2025). Skyborn: Climate Data Analysis Toolkit. Zenodo. https://doi.org/10.5281/zenodo.18075252
