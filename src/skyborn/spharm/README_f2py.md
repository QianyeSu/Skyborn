# F2PY Integration for SPHEREPACK

This document explains how to generate Python interfaces for the modernized Fortran 90 SPHEREPACK files.

## Quick Start

### 1. Generate F2PY Signature File

Run the script from the spharm directory:

```bash
# From the spharm directory (containing src/ folder)
python fix_f2py_symbols.py
```

This will generate `spherepack_f90.pyf` containing interfaces for all `.f90` files.

### 2. Compile the Extension Module

```bash
# Compile all .f90 files into a Python extension
python -m numpy.f2py -c spherepack_f90.pyf src/*.f90
```

### 3. Use in Python

```python
import spherepack_f90

# Use the modernized Fortran 90 functions
# Example: spherepack_f90.shsgs(...)
```

## File Structure

```
spharm/
├── fix_f2py_symbols.py          # Script to generate .pyf files
├── spherepack_f90.pyf            # Generated F2PY signature file
├── README_f2py.md                # This documentation
└── src/
    ├── *.f90                     # Modernized Fortran 90 files
    ├── *.f                       # Original FORTRAN 77 files
    └── _spherepack.pyf           # Existing interface for .f files
```

## Comparison: .f vs .f90 Interfaces

| Feature | _spherepack.pyf (.f files) | spherepack_f90.pyf (.f90 files) |
|---------|---------------------------|----------------------------------|
| **Language** | FORTRAN 77 | Modern Fortran 2008+ |
| **Performance** | Standard | Optimized (vectorized, parallel) |
| **Maintainability** | Basic | High (explicit interfaces) |
| **Memory Safety** | Basic | Enhanced (intent declarations) |

## Troubleshooting

### Common Issues

1. **"Source directory not found"**
   - Ensure you run the script from the correct directory
   - The script looks for `src/` folder containing `.f90` files

2. **F2PY compilation errors**
   - Check that all `.f90` files have valid syntax
   - Ensure dependencies (iso_fortran_env, etc.) are available

3. **Import errors**
   - Verify the compiled extension is in your Python path
   - Check that numpy and f2py are properly installed

### Manual Path Override

If automatic path detection fails, modify the script:

```python
# In fix_f2py_symbols.py, change the main() function:
src_dir = Path("your/custom/path/to/src")
```

## Advanced Usage

### Selective Compilation

To compile only specific functions:

```bash
python -m numpy.f2py -h custom.pyf src/shsgs.f90 src/shsgc.f90 only: shsgs shsgc
python -m numpy.f2py -c custom.pyf src/shsgs.f90 src/shsgc.f90
```

### Debug Mode

Add debug flags for troubleshooting:

```bash
python -m numpy.f2py -c --debug-capi spherepack_f90.pyf src/*.f90
```

## Contributing

When adding new `.f90` files:

1. Place them in the `src/` directory
2. Re-run `python fix_f2py_symbols.py` to update the interface
3. Test the compilation and Python import

## Support

For issues related to:
- **F2PY**: Check numpy.f2py documentation
- **Fortran compilation**: Verify your Fortran compiler setup
- **SPHEREPACK algorithms**: Refer to original SPHEREPACK documentation
