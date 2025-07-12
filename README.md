
<p align="center">
  <a href="https://github.com/QianyeSu/Skyborn" target="_blank">
    <img src="docs/source/_static/SkyBornLogo.svg" alt="Skyborn Logo" width="400"/>
  </a>
</p>

[![PyPI version](https://badge.fury.io/py/skyborn.svg)](https://badge.fury.io/py/skyborn)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/skyborn)](https://pypi.org/project/skyborn/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/skyborn)](https://pypi.org/project/skyborn/)
[![codecov](https://codecov.io/gh/QianyeSu/Skyborn/graph/badge.svg?token=YOUR_TOKEN_HERE)](https://codecov.io/gh/QianyeSu/Skyborn)
[![License](https://img.shields.io/github/license/QianyeSu/Skyborn)](https://github.com/QianyeSu/Skyborn/blob/main/LICENSE)
[![Tests](https://github.com/QianyeSu/Skyborn/actions/workflows/stable-ci.yml/badge.svg)](https://github.com/QianyeSu/Skyborn/actions/workflows/stable-ci.yml)
[![Platform](https://img.shields.io/badge/platform-Windows-blue)](https://github.com/QianyeSu/Skyborn)
[![Code style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-brightgreen)](https://skyborn.readthedocs.io/en/latest/)


## System Requirements

**Operating System:** 🖥️ **Cross-Platform**

This package supports Windows, Linux, and macOS. However, it has been primarily developed and tested on Windows.

**Note:** While the package can be installed on different platforms, some Windows-specific features may not work on other operating systems.

## 🚀 Installation

### 📦 Quick Installation (Recommended)
```bash
pip install skyborn
```

### ⚡ High-Performance Installation (Maximum Speed)

If you need **maximum performance** and have a Fortran compiler:

```bash
# Source installation with native CPU optimization
pip install --no-binary=skyborn skyborn

# Or development installation
git clone https://github.com/QianyeSu/Skyborn.git
cd Skyborn
pip install -e .
```

**Performance benefit**: Source installation can be **20-50% faster** as it compiles with `-march=native` optimization specifically for your CPU.

### 🌐 Spherical Harmonics (spharm) Module Installation

For users who need the **spherical harmonic transforms** functionality (adapted from [pyspharm](https://github.com/jswhit/pyspharm)), we recommend source installation for optimal performance:

```bash
# Recommended: Source installation for maximum performance
pip install -U --index-url https://pypi.org/simple/ --no-binary=skyborn skyborn

# Alternative: Standard installation (may have reduced performance for spharm)
pip install -U skyborn
```

The `spharm` module requires Fortran compilation and benefits significantly from CPU-specific optimizations. Source installation ensures the SPHEREPACK library is compiled with optimal flags for your system.

### 📋 Requirements for Source Installation
```bash
# Ubuntu/Debian
sudo apt-get install gfortran

# macOS
brew install gcc

# Windows (using conda)
conda install -c conda-forge fortran-compiler
```
