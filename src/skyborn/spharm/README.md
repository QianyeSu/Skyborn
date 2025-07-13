# Skyborn SPHARM Module

## Overview

The Skyborn SPHARM module provides high-performance spherical harmonic transforms for atmospheric and oceanic modeling applications. This module is specifically optimized for Python 3.9+ with modern Fortran backend acceleration.

## Attribution and Source

**This code is adapted from the original pyspharm project:**
- **Original Repository**: https://github.com/jswhit/pyspharm
- **Original Author**: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
- **Original License**: SPHEREPACK license (see LICENSE.spherepack)

### Modifications and Enhancements

The Skyborn version includes significant improvements over the original pyspharm:

1. **Modern Build System**: Updated to use Meson build system for better cross-platform support
2. **Performance Optimization**: Enhanced Fortran code with OpenMP parallelization and vectorization
3. **Python 3.9+ Compatibility**: Updated for modern Python environments
4. **Integration**: Seamlessly integrated into the Skyborn ecosystem
5. **Documentation**: Comprehensive documentation and examples

## Installation

The SPHARM module is automatically installed as part of the Skyborn package:

```bash
pip install skyborn
```

## Quick Example

```python
from skyborn.spharm import Spharmt

# Create spherical harmonic transform instance
sht = Spharmt(nlon=144, nlat=73, gridtype='gaussian')

# Perform grid to spectral transform
spectral_coeffs = sht.grdtospec(grid_data)

# Perform spectral to grid transform
reconstructed_data = sht.spectogrd(spectral_coeffs)
```

## Dependencies

- Python 3.9 or higher
- NumPy 1.21.0 or higher
- Fortran compiler (gfortran recommended)

## License

This module maintains the original SPHEREPACK licensing terms. See `LICENSE.spherepack` for complete details.

## Credits

- **Original SPHEREPACK Development**: NCAR (National Center for Atmospheric Research)
- **Original Python Interface**: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
- **Skyborn Integration and Optimization**: Qianye Su <suqianye2000@gmail.com>

## Support

For issues specific to the Skyborn implementation:
- GitHub Issues: https://github.com/QianyeSu/Skyborn/issues
- Documentation: https://skyborn.readthedocs.io/

For questions about the original pyspharm algorithms, please refer to:
- Original Repository: https://github.com/jswhit/pyspharm
