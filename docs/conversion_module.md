# GRIB to NetCDF Conversion Module

The `skyborn.conversion` module provides a Python interface to the eccodes `grib_to_netcdf` tool for converting GRIB files to NetCDF format.

## Installation Requirements

Make sure you have the eccodes library installed:

```bash
# conda installation
conda install -c conda-forge eccodes

# or using pip
pip install eccodes
```

## Quick Start

### Basic Usage

```python
import skyborn

# Simple conversion
skyborn.grib2nc('input.grib', 'output.nc')

# or using full function name
skyborn.convert_grib_to_nc_simple('input.grib', 'output.nc')
```

### Advanced Usage

```python
from skyborn.conversion import convert_grib_to_nc

# High precision conversion
convert_grib_to_nc(
    grib_files='input.grib',
    output_file='output.nc',
    data_type='NC_FLOAT',           # High precision floating point
    ignore_keys=['type', 'step'],   # Ignore specific keys
    unlimited_dimension='time',     # Set time as unlimited dimension
    file_kind=4,                    # netCDF-4 format
    deflate_level=6,                # Compression level
    shuffle=True                    # Enable shuffle compression
)
```

### Batch Conversion

```python
from skyborn.conversion import batch_convert_grib_to_nc

# Convert all GRIB files in a directory
converted_files = batch_convert_grib_to_nc(
    input_directory='/path/to/grib_files',
    output_directory='/path/to/output',
    pattern='*.grib*',
    high_precision=True,
    compress=True
)

print(f"Converted {len(converted_files)} files successfully")
```

## Available Functions

### Main Functions

1. **`convert_grib_to_nc()`** - Full-featured conversion function with all options
2. **`convert_grib_to_nc_simple()`** - Simplified interface with common presets
3. **`batch_convert_grib_to_nc()`** - Batch conversion for directories

### Convenience Aliases

- **`grib2nc()`** - Alias for `convert_grib_to_nc_simple()`
- **`grib_to_netcdf()`** - Alias for `convert_grib_to_nc()`

## Function Parameters

### convert_grib_to_nc()

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `grib_files` | str/Path/List | - | Input GRIB file(s) |
| `output_file` | str/Path | - | Output NetCDF file |
| `ignore_keys` | List[str] | `['method', 'type', 'stream', 'refdate', 'hdate']` | GRIB keys to ignore |
| `split_keys` | List[str] | `['param', 'expver']` | Keys to split according to |
| `data_type` | str | `'NC_SHORT'` | NetCDF data type |
| `file_kind` | int | `2` | NetCDF file format (1-4) |
| `deflate_level` | int | `None` | Compression level (0-9) |
| `shuffle` | bool | `False` | Enable shuffle compression |
| `unlimited_dimension` | str | `None` | Set unlimited dimension |

### convert_grib_to_nc_simple()

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `grib_file` | str/Path | - | Input GRIB file |
| `output_file` | str/Path | - | Output NetCDF file |
| `high_precision` | bool | `False` | Use NC_FLOAT instead of NC_SHORT |
| `compress` | bool | `False` | Enable compression (netCDF-4) |

## Data Types

- **`NC_BYTE`** - 8-bit signed integer
- **`NC_SHORT`** - 16-bit signed integer (default)
- **`NC_INT`** - 32-bit signed integer
- **`NC_FLOAT`** - 32-bit floating point
- **`NC_DOUBLE`** - 64-bit floating point

## File Formats

- **1** - netCDF classic file format
- **2** - netCDF 64-bit classic file format (default)
- **3** - netCDF-4 file format
- **4** - netCDF-4 classic model file format

## Examples

### ERA5 Data Conversion

```python
from skyborn.conversion import convert_grib_to_nc

# Optimized settings for ERA5 data
convert_grib_to_nc(
    grib_files='era5_data.grib',
    output_file='era5_data.nc',
    ignore_keys=['method', 'type', 'stream'],  # Common ERA5 ignore keys
    split_keys=['param', 'levtype'],           # Split by parameter and level type
    data_type='NC_FLOAT',                      # Recommended precision for ERA5
    unlimited_dimension='time',                # Unlimited time dimension
    file_kind=4,                               # netCDF-4 format
    deflate_level=4                            # Moderate compression
)
```

### Multiple Files Conversion

```python
# Convert and merge multiple GRIB files
convert_grib_to_nc(
    grib_files=['file1.grib', 'file2.grib', 'file3.grib'],
    output_file='merged_output.nc',
    data_type='NC_DOUBLE',
    split_keys=['param'],
    reference_date='20240101'
)
```

### High Compression Setup

```python
# Maximum compression for storage
convert_grib_to_nc(
    grib_files='input.grib',
    output_file='compressed_output.nc',
    data_type='NC_FLOAT',
    file_kind=4,          # netCDF-4 required for compression
    deflate_level=9,      # Maximum compression
    shuffle=True          # Enable shuffle for better compression
)
```

## Error Handling

The module provides custom exception handling:

```python
from skyborn.conversion import convert_grib_to_nc, GribToNetCDFError

try:
    convert_grib_to_nc('input.grib', 'output.nc')
except GribToNetCDFError as e:
    print(f"Conversion failed: {e}")
except FileNotFoundError as e:
    print(f"Input file not found: {e}")
```

## Notes and Limitations

1. **Grid Support**: Only regular lat/lon grids and regular Gaussian grids are supported
2. **Dependencies**: Requires eccodes library with `grib_to_netcdf` command-line tool
3. **Performance**: Large files may take significant time to convert
4. **Memory**: NetCDF-4 with compression requires more memory during conversion

## Troubleshooting

### Common Issues

1. **"grib_to_netcdf command not found"**
   - Install eccodes: `conda install -c conda-forge eccodes`
   - Ensure eccodes is in your PATH

2. **"GRIB file not found"**
   - Check file path and permissions
   - Verify file exists and is readable

3. **"Conversion failed with return code 1"**
   - Check GRIB file format compatibility
   - Verify grid type is supported
   - Check available disk space

### Performance Tips

1. Use appropriate data types (`NC_SHORT` for integer data, `NC_FLOAT` for precision)
2. Enable compression for large files (`compress=True`)
3. Use batch conversion for multiple files
4. Set appropriate deflate levels (4-6 for good compression/speed balance)

## Integration with Other Skyborn Modules

The conversion module works seamlessly with other skyborn modules:

```python
import skyborn
import xarray as xr

# Convert GRIB to NetCDF
skyborn.grib2nc('input.grib', 'output.nc')

# Load with xarray and use skyborn functions
ds = xr.open_dataset('output.nc')
gradient = skyborn.calculate_gradient(ds['temperature'], 'latitude')
```
