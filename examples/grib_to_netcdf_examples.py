"""
Examples: Using skyborn.conversion module for GRIB to NetCDF conversion

This file demonstrates how to use the new conversion module in the skyborn library
to convert GRIB files to NetCDF format.
"""

import skyborn
from skyborn.conversion import (
    batch_convert_grib_to_nc,
    convert_grib_to_nc,
    convert_grib_to_nc_simple,
)


def example_basic_conversion():
    """Basic conversion example"""
    print("=== Basic GRIB to NetCDF conversion ===")

    # Basic conversion
    try:
        output_file = convert_grib_to_nc_simple(
            grib_file="input.grib", output_file="output.nc"
        )
        print(f"Conversion completed: {output_file}")
    except Exception as e:
        print(f"Conversion failed: {e}")


def example_advanced_conversion():
    """Advanced conversion example"""
    print("\n=== Advanced GRIB to NetCDF conversion ===")

    try:
        output_file = convert_grib_to_nc(
            grib_files="input.grib",
            output_file="output_advanced.nc",
            data_type="NC_FLOAT",  # High precision
            ignore_keys=["type", "step"],  # Ignore specific keys
            unlimited_dimension="time",  # Set unlimited dimension
            file_kind=4,  # netCDF-4 format
            deflate_level=6,  # Compression level
            shuffle=True,  # Enable shuffle compression
        )
        print(f"Advanced conversion completed: {output_file}")
    except Exception as e:
        print(f"Conversion failed: {e}")


def example_batch_conversion():
    """Batch conversion example"""
    print("\n=== Batch GRIB to NetCDF conversion ===")

    try:
        converted_files = batch_convert_grib_to_nc(
            input_directory="/path/to/grib_files",
            output_directory="/path/to/output",
            pattern="*.grib*",
            high_precision=True,
            compress=True,
        )
        print(f"Batch conversion completed, {len(converted_files)} files converted")
    except Exception as e:
        print(f"Batch conversion failed: {e}")


def example_multiple_files():
    """Multiple files conversion example"""
    print("\n=== Multiple files merged conversion ===")

    try:
        output_file = convert_grib_to_nc(
            grib_files=["file1.grib", "file2.grib", "file3.grib"],
            output_file="merged_output.nc",
            data_type="NC_DOUBLE",
            split_keys=["param"],  # Split by parameter
            reference_date="20240101",
        )
        print(f"Multiple files conversion completed: {output_file}")
    except Exception as e:
        print(f"Conversion failed: {e}")


def example_era5_conversion():
    """ERA5 data conversion example"""
    print("\n=== ERA5 GRIB data conversion ===")

    try:
        # Settings suitable for ERA5 data
        output_file = convert_grib_to_nc(
            grib_files="era5_data.grib",
            output_file="era5_data.nc",
            # Commonly ignored keys for ERA5
            ignore_keys=["method", "type", "stream"],
            # Split by parameter and level type
            split_keys=["param", "levtype"],
            data_type="NC_FLOAT",  # Recommended precision for ERA5
            unlimited_dimension="time",  # Unlimited time dimension
            file_kind=4,  # netCDF-4 format
            deflate_level=4,  # Moderate compression
        )
        print(f"ERA5 conversion completed: {output_file}")
    except Exception as e:
        print(f"Conversion failed: {e}")


if __name__ == "__main__":
    print("Skyborn Conversion Module Usage Examples")
    print(f"Skyborn Version: {skyborn.__version__}")

    # Run examples
    example_basic_conversion()
    example_advanced_conversion()
    example_batch_conversion()
    example_multiple_files()
    example_era5_conversion()

    print("\n=== Available Conversion Functions ===")
    print("1. convert_grib_to_nc() - Full-featured conversion function")
    print("2. convert_grib_to_nc_simple() - Simplified interface")
    print("3. batch_convert_grib_to_nc() - Batch conversion")
    print("4. grib2nc() - Alias for convert_grib_to_nc_simple")
    print("5. grib_to_netcdf() - Alias for convert_grib_to_nc")
