# Double Precision Usage Analysis - Summary Table

| Category | File Name | Usage Count | Primary Purpose |
|----------|-----------|-------------|-----------------|
| **Heavy Usage** | `sphcom.f90` | 432 | Core spherical harmonic computations |
| | `fortran/sphcom_backup.f90` | 432 | Backup of core spherical harmonic computations |
| | `gaqd.f90` | 50+ | Gaussian quadrature points and weights |
| **Moderate Usage** | `vhags.f90` | 9 | Vector spherical harmonic analysis |
| | `vhsgs.f90` | 9 | Vector spherical harmonic synthesis |
| **Light Usage** | `shags.f90` | 12 | Spherical harmonic analysis (Gaussian grid) |
| | `shsgs.f90` | 11 | Spherical harmonic synthesis (Gaussian grid) |
| | `shagc.f90` | 9 | Spherical harmonic analysis (colatitude-longitude grid) |
| | `shsgc.f90` | 9 | Spherical harmonic synthesis (colatitude-longitude grid) |
| | `vhaes.f90` | 3 | Vector spherical harmonic analysis (equally spaced) |
| | `vhses.f90` | 3 | Vector spherical harmonic synthesis (equally spaced) |
| | `vhsgc.f90` | 3 | Vector spherical harmonic synthesis (colatitude-longitude grid) |
| | `shaec.f90` | 2 | Spherical harmonic analysis (equally spaced, colatitude-longitude grid) |
| | `shaes.f90` | 2 | Spherical harmonic analysis (equally spaced) |
| | `shsec.f90` | 2 | Spherical harmonic synthesis (equally spaced, colatitude-longitude grid) |
| | `shses.f90` | 2 | Spherical harmonic synthesis (equally spaced) |
| | `vhaec.f90` | 2 | Vector spherical harmonic analysis (equally spaced, colatitude-longitude grid) |
| | `vhagc.f90` | 2 | Vector spherical harmonic analysis (colatitude-longitude grid) |
| | `vhsec.f90` | 1 | Vector spherical harmonic synthesis (equally spaced, colatitude-longitude grid) |
| **No Usage** | `alf.f90` | 0 | Associated Legendre functions |
| | `getlegfunc.f90` | 0 | Get Legendre functions |
| | `hrfft.f90` | 0 | Real FFT routines |
| | `ihgeod.f90` | 0 | Inverse geodesic problem |
| | `invlap.f90` | 0 | Inverse Laplacian |
| | `lap.f90` | 0 | Laplacian |
| | `multsmoothfact.f90` | 0 | Multiplication smoothing factor |
| | `onedtotwod.f90` | 0 | 1D to 2D conversion |
| | `onedtotwod_vrtdiv.f90` | 0 | 1D to 2D conversion (vorticity/divergence) |
| | `specintrp.f90` | 0 | Spectral interpolation |
| | `twodtooned.f90` | 0 | 2D to 1D conversion |
| | `twodtooned_vrtdiv.f90` | 0 | 2D to 1D conversion (vorticity/divergence) |

## Double Precision Usage Patterns

### 1. **real64 kind parameter**: `integer, parameter :: real64 = kind(1.0d0)`
### 2. **real64 literals**: `1.0_real64`, `2.0_real64`
### 3. **Double precision declarations**: `double precision, intent(inout) :: dwork(*)`
### 4. **real64 type declarations**: `real(kind=real64) :: variable`
### 5. **real64 type casting**: `real(expression, kind=real64)`
### 6. **Traditional double literals**: `1.0d0`, `2.0d0`
### 7. **real(8) declarations**: `real(8), intent(out) :: array`