# Geostrophic Wind Calculation Module

High-performance Fortran implementations for computing geostrophic wind components from geopotential height fields.

## Files Overview

### z2geouv.f90 (Default - Single Precision)
**Function**: Geopotential height to geostrophic wind conversion algorithm (single precision)

**Recommended for**:
- Production use with float32 input data
- High-performance applications requiring maximum speed
- Large-scale processing where memory efficiency is critical
- Standard meteorological analysis (sufficient precision)

**Performance**: 1.4-1.8× faster than double precision version

### z2geouv_dp.f90 (High Precision - Double Precision)
**Function**: Geopotential height to geostrophic wind conversion algorithm (double precision)

**Recommended for**:
- Production applications requiring maximum precision
- Scientific research with high accuracy requirements
- Quality control and reference calculations
- Applications with double precision input data

## Detailed Description

**Core Subroutines**:
- **Z2GEOUV**: Main interface program, automatically determines latitude ordering
- **Z2GEOUV_3D**: 3D version for processing multiple pressure levels
- **Z2GEOUV_4D**: 4D version for processing pressure levels × time steps
- **ZUVNEW**: North-to-south data reordering processing program
- **Z2GUV**: Geostrophic wind core calculation program (requires south-to-north ordering)

**Input Parameters**:
- **Z**: Geopotential height field (NLAT×MLON, units: gpm)
- **NLAT,MLON**: Number of latitude and longitude grid points
- **ZMSG**: Missing value code
- **GLAT**: Latitude array (NLAT, degrees, south-to-north preferred)
- **GLON**: Longitude array (MLON, degrees)
- **IOPT**: Boundary condition option (0=non-periodic, 1=periodic)

**Output Parameters**:
- **UG,VG**: Geostrophic wind U and V components (NLAT×MLON, units: m/s)

## Mathematical Foundation

**Basic Geostrophic Wind Equations**:
- Geostrophic balance: $f \vec{V_g} = -\nabla \Phi \times \hat{k}$
- U component: $U_g = -\frac{g}{f} \frac{\partial Z}{\partial y} = -\frac{g}{f} \frac{1}{R} \frac{\partial Z}{\partial \phi}$
- V component: $V_g = \frac{g}{f} \frac{\partial Z}{\partial x} = \frac{g}{f} \frac{1}{R \cos\phi} \frac{\partial Z}{\partial \lambda}$

**Physical Constants**:
- Gravitational acceleration: $g = 9.80616 \text{ m/s}^2$
- Earth radius: $R_e = 6371220 \text{ m}$
- Earth rotation angular velocity: $\Omega = 7.292 \times 10^{-5} \text{ rad/s}$
- Coriolis parameter: $f = 2\Omega \sin\phi$

**Calculation Methods**:
- Gravity/Coriolis force ratio: $\frac{g}{f} = \frac{g}{2\Omega \sin\phi}$
- Spatial difference distances:
  - Zonal: $\Delta x = 2R \Delta\lambda \cos\phi$ (units: m)
  - Meridional: $\Delta y = R \Delta\phi$ (units: m)
- Central difference scheme: $\frac{\partial Z}{\partial x} \approx \frac{Z(i+1,j) - Z(i-1,j)}{2\Delta x}$

## Boundary Handling Strategies

**Latitude Boundaries**: Use one-sided differences with appropriate distance weighting

**Longitude Boundaries**:
- **IOPT=1**: Periodic boundary conditions $Z(0,j) = Z(\text{MLON},j)$
- **IOPT=0**: Copy adjacent point values, multiply V component by 2 to compensate for one-sided difference

## Special Region Handling

- **Near equator**: When $|\phi| < 10^{-5}$, Coriolis term set to missing value
- **Near poles**: When $|\phi| > 89.999°$, zonal distance set to 0
- **Missing data**: Automatic detection and proper handling of missing values

## Data Ordering Requirements

- **Algorithm requirement**: Latitude ordering from south to north
- **Automatic detection**: Input data ordering automatically detected
- **Flexible input**: North-to-south data automatically reordered via ZUVNEW subroutine

## Key Features

### Performance Optimizations
- **OpenMP SIMD vectorization**: Enhanced with `simdlen(8)` directives
- **Precomputed inverses**: Eliminates divisions in computational hot loops
- **Cache-friendly memory access**: Optimized loop structures
- **Multi-dimensional support**: Efficient 3D/4D processing with OpenMP parallelization

### Precision Options
- **Single precision (z2geouv.f90)**:
  - 4-byte floats, ~7 decimal digits accuracy
  - 1.4-1.8× performance improvement
  - Perfect for meteorological applications
- **Double precision (z2geouv_dp.f90)**:
  - 8-byte floats, ~15 decimal digits accuracy
  - Reference quality calculations

### Robustness
- **Automatic data ordering detection and processing**
- **Flexible boundary condition settings**
- **Complete set of geophysical parameters**
- **Accurate spherical geometry calculations**
- **Comprehensive missing value handling**

## Applications

- Numerical weather prediction systems
- Atmospheric dynamics analysis
- Wind field diagnostics and visualization
- Climate model validation
- Meteorological data processing pipelines
- Real-time weather analysis systems

## Usage Recommendation

**For Production Systems**: Use **z2geouv.f90** (single precision) as the default choice for optimal performance with meteorological data quality requirements.

**For Research/Validation**: Use **z2geouv_dp.f90** (double precision) when maximum numerical precision is required for scientific analysis or reference calculations.

## Performance Benchmark Results

Based on comprehensive testing with various data sizes:

| Data Size | Single Precision | Double Precision | Speedup |
|-----------|-----------------|------------------|---------|
| 100×160   | 0.090 ms        | 0.160 ms         | 1.78×   |
| 200×320   | 0.330 ms        | 0.470 ms         | 1.42×   |
| 300×480   | 0.860 ms        | 1.270 ms         | 1.48×   |
| 400×640   | 1.810 ms        | 2.800 ms         | 1.55×   |

Maximum precision difference: ~5.78E-03 m/s (acceptable for meteorological applications)
