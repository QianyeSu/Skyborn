# Skyborn SPHARM Module

## Overview

We are pleased to announce the **Skyborn SPHARM module**, a high-performance Python interface to NCAR's SPHEREPACK Fortran library for spherical harmonic transforms. This module is specifically designed for atmospheric and oceanic modeling applications.

## Key Features

- **Modern Build System**: Updated to use Meson build system for Python 3.9+ compatibility
- **Optimized Performance**: Enhanced Fortran code with OpenMP parallelization and vectorization
- **Seamless Integration**: Fully integrated into the Skyborn ecosystem
- **Cross-Platform**: Supports Linux, macOS, and Windows

## System Requirements

### Dependencies
- **Python**: 3.9 or higher
- **NumPy**: 1.21.0 or higher
- **Fortran Compiler**: gfortran or compatible compiler supported by f2py
- **Build Tools**: Meson, Ninja (automatically installed)

### Supported Platforms
- Linux (x86_64)
- macOS (Intel and Apple Silicon)
- Windows (x86_64)

## Installation

### Quick Installation (Recommended)
```bash
pip install skyborn
```

### Development Installation
```bash
git clone https://github.com/QianyeSu/Skyborn.git
cd Skyborn
python -m pip install -e .
```

### Custom Fortran Compiler
If you need to specify a custom Fortran compiler:
```bash
# Set environment variable before installation
export FC=gfortran
python -m pip install -e .
```

## Documentation

Comprehensive documentation is available in multiple formats:

- **HTML Documentation**: Browse to `html/index.html` in your installation
- **API Reference**: Complete API documentation with examples
- **Example Programs**: Practical examples in the `examples/` directory
- **Online Documentation**: Visit [https://skyborn.readthedocs.io/](https://skyborn.readthedocs.io/)

## Quick Start Example

```python
from skyborn.spharm import Spharmt

# Create spherical harmonic transform instance
sht = Spharmt(nlon=144, nlat=73, gridtype='gaussian')

# Perform grid to spectral transform
spectral_coeffs = sht.grdtospec(grid_data)

# Perform spectral to grid transform
reconstructed_data = sht.spectogrd(spectral_coeffs)
```

## Performance Optimizations

Our SPHARM module includes several performance enhancements:

- **Optimized Fortran Code**: Modern Fortran 90+ with compiler optimizations
- **Parallel Processing**: OpenMP support for multi-core acceleration
- **Vectorization**: SIMD instructions for improved computational efficiency
- **Memory Management**: Efficient memory allocation and reuse

## License and Attribution

### Python Binding License
The Python interface and Skyborn integration are provided under the BSD-3-Clause license:

```
Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation.

THE AUTHORS DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
EVENT SHALL THE AUTHORS BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
```

### SPHEREPACK Fortran Library
The underlying SPHEREPACK Fortran library has its own licensing terms. Please refer to `LICENSE.spherepack` for complete details.

## Support and Contributing

- **Issues**: Report bugs and feature requests on [GitHub Issues](https://github.com/QianyeSu/Skyborn/issues)
- **Discussions**: Join our community discussions on [GitHub Discussions](https://github.com/QianyeSu/Skyborn/discussions)
- **Contributing**: See our contributing guidelines for development information

## Credits

**Original SPHEREPACK Development**: Jeff Whitaker <Jeffrey.S.Whitaker@noaa.gov>

**Skyborn Integration and Optimization**: Qianye Su <suqianye2000@gmail.com>

**Special Thanks**: NCAR for the original SPHEREPACK library

---

*The Skyborn SPHARM module represents our commitment to providing high-performance, accessible tools for the atmospheric and oceanic sciences community.*
