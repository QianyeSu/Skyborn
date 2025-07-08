"""
GRIB to NetCDF conversion utilities using eccodes grib_to_netcdf tool.

This module provides a Python interface to the grib_to_netcdf command-line tool
from the eccodes library for converting GRIB files to NetCDF format.

Author: Qianye Su
Email: suqianye2000@gmail.com

References:
    - eccodes grib_to_netcdf tool by Sandor Kertesz, modified by Shahram Najm
    - https://confluence.ecmwf.int/display/ECC/grib_to_netcdf
"""

import subprocess
import os
import sys
from pathlib import Path
from typing import Union, Optional, List, Literal
import shutil
import logging

# Set up logging
logger = logging.getLogger(__name__)

# Type definitions
PathLike = Union[str, Path]
DataType = Literal["NC_BYTE", "NC_SHORT", "NC_INT", "NC_FLOAT", "NC_DOUBLE"]
FileKind = Literal[1, 2, 3, 4]


class GribToNetCDFError(Exception):
    """Exception raised when grib_to_netcdf conversion fails."""

    pass


def _check_grib_to_netcdf_available() -> bool:
    """Check if grib_to_netcdf command is available."""
    return shutil.which("grib_to_netcdf") is not None


def _validate_grib_files(grib_files: List[PathLike]) -> List[Path]:
    """Validate that GRIB files exist and return Path objects."""
    validated_files = []
    for file in grib_files:
        file_path = Path(file)
        if not file_path.exists():
            raise FileNotFoundError(f"GRIB file not found: {file_path}")
        if not file_path.is_file():
            raise ValueError(f"Path is not a file: {file_path}")
        validated_files.append(file_path)
    return validated_files


def _build_grib_to_netcdf_command(
    output_file: PathLike,
    grib_files: List[PathLike],
    ignore_keys: Optional[List[str]] = None,
    split_keys: Optional[List[str]] = None,
    reference_date: Optional[str] = None,
    data_type: DataType = "NC_SHORT",
    no_time_validity: bool = False,
    force: bool = False,
    multi_field_off: bool = False,
    file_kind: FileKind = 2,
    deflate_level: Optional[int] = None,
    shuffle: bool = False,
    unlimited_dimension: Optional[str] = None,
) -> List[str]:
    """Build the grib_to_netcdf command with all options."""

    cmd = ["grib_to_netcdf"]

    # Ignore keys
    if ignore_keys:
        cmd.extend(["-I", ",".join(ignore_keys)])

    # Split keys
    if split_keys:
        cmd.extend(["-S", ",".join(split_keys)])

    # Reference date
    if reference_date:
        cmd.extend(["-R", reference_date])

    # Data type
    cmd.extend(["-D", data_type])

    # No time validity
    if no_time_validity:
        cmd.append("-T")

    # Force execution
    if force:
        cmd.append("-f")

    # Multi-field support off
    if multi_field_off:
        cmd.append("-M")

    # File kind
    cmd.extend(["-k", str(file_kind)])

    # Deflate compression (only for netCDF-4)
    if deflate_level is not None and file_kind in [3, 4]:
        if not (0 <= deflate_level <= 9):
            raise ValueError("Deflate level must be between 0 and 9")
        cmd.extend(["-d", str(deflate_level)])

    # Shuffle data
    if shuffle:
        cmd.append("-s")

    # Unlimited dimension
    if unlimited_dimension:
        cmd.extend(["-u", unlimited_dimension])

    # Output file
    cmd.extend(["-o", str(output_file)])

    # Input GRIB files
    cmd.extend([str(f) for f in grib_files])

    return cmd


def convert_grib_to_nc(
    grib_files: Union[PathLike, List[PathLike]],
    output_file: PathLike,
    ignore_keys: Optional[List[str]] = None,
    split_keys: Optional[List[str]] = None,
    reference_date: Optional[str] = None,
    data_type: DataType = "NC_SHORT",
    no_time_validity: bool = False,
    force: bool = False,
    multi_field_off: bool = False,
    file_kind: FileKind = 2,
    deflate_level: Optional[int] = None,
    shuffle: bool = False,
    unlimited_dimension: Optional[str] = None,
    verbose: bool = True,
) -> Path:
    """
    Convert GRIB file(s) to NetCDF format using eccodes grib_to_netcdf tool.

    This function provides a Python interface to the grib_to_netcdf command-line tool.
    Only regular lat/lon grids and regular Gaussian grids are supported.

    Parameters
    ----------
    grib_files : PathLike or List[PathLike]
        Path(s) to input GRIB file(s). Can be a single file or list of files.
    output_file : PathLike
        Path to output NetCDF file.
    ignore_keys : List[str], optional
        GRIB keys to ignore. Default: ['method', 'type', 'stream', 'refdate', 'hdate']
    split_keys : List[str], optional
        Keys to split according to. Default: ['param', 'expver']
    reference_date : str, optional
        Reference date in YYYYMMDD format. Default: '19000101'
    data_type : {'NC_BYTE', 'NC_SHORT', 'NC_INT', 'NC_FLOAT', 'NC_DOUBLE'}, default 'NC_SHORT'
        NetCDF data type for output.
    no_time_validity : bool, default False
        Don't use time of validity. Creates separate time dimensions instead.
    force : bool, default False
        Force execution not to fail on error.
    multi_field_off : bool, default False
        Turn off support for multiple fields in single GRIB message.
    file_kind : {1, 2, 3, 4}, default 2
        NetCDF file format:

        * 1 = netCDF classic file format
        * 2 = netCDF 64 bit classic file format (Default)
        * 3 = netCDF-4 file format
        * 4 = netCDF-4 classic model file format

    deflate_level : int, optional
        Compression level (0-9) for netCDF-4 output format only.
    shuffle : bool, default False
        Shuffle data before deflation compression.
    unlimited_dimension : str, optional
        Set dimension to be unlimited (e.g., 'time').
    verbose : bool, default True
        Print verbose output during conversion.

    Returns
    -------
    Path
        Path to the created NetCDF file.

    Raises
    ------
    GribToNetCDFError
        If grib_to_netcdf tool is not available or conversion fails.
    FileNotFoundError
        If input GRIB files don't exist.
    ValueError
        If invalid parameters are provided.

    Examples
    --------
    Basic usage:

    >>> convert_grib_to_nc('input.grib', 'output.nc')

    Convert multiple files with custom settings:

    >>> convert_grib_to_nc(
    ...     ['file1.grib', 'file2.grib'],
    ...     'output.nc',
    ...     data_type='NC_FLOAT',
    ...     ignore_keys=['type', 'step'],
    ...     unlimited_dimension='time'
    ... )

    High precision with compression:

    >>> convert_grib_to_nc(
    ...     'input.grib',
    ...     'output.nc',
    ...     data_type='NC_DOUBLE',
    ...     file_kind=4,  # netCDF-4 format
    ...     deflate_level=6,
    ...     shuffle=True
    ... )
    """

    # Check if grib_to_netcdf is available
    if not _check_grib_to_netcdf_available():
        raise GribToNetCDFError(
            "grib_to_netcdf command not found. "
            "Please ensure eccodes is installed and available in PATH."
        )

    # Convert single file to list
    if isinstance(grib_files, (str, Path)):
        grib_files = [grib_files]

    # Validate input files
    validated_grib_files = _validate_grib_files(grib_files)

    # Convert output path
    output_path = Path(output_file)

    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Set default values if not provided
    if ignore_keys is None:
        ignore_keys = ["method", "type", "stream", "refdate", "hdate"]
    if split_keys is None:
        split_keys = ["param", "expver"]

    # Build command
    cmd = _build_grib_to_netcdf_command(
        output_file=output_path,
        grib_files=validated_grib_files,
        ignore_keys=ignore_keys,
        split_keys=split_keys,
        reference_date=reference_date,
        data_type=data_type,
        no_time_validity=no_time_validity,
        force=force,
        multi_field_off=multi_field_off,
        file_kind=file_kind,
        deflate_level=deflate_level,
        shuffle=shuffle,
        unlimited_dimension=unlimited_dimension,
    )

    if verbose:
        print("Converting GRIB to NetCDF...")
        print(f"Input files: {[str(f) for f in validated_grib_files]}")
        print(f"Output file: {output_path}")
        print(f"Command: {' '.join(cmd)}")

    try:
        # Execute the command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if verbose and result.stdout:
            print("STDOUT:", result.stdout)

        if verbose:
            print("✅ Conversion completed successfully!")
            print(f"Output saved to: {output_path.absolute()}")

        return output_path

    except subprocess.CalledProcessError as e:
        error_msg = f"grib_to_netcdf failed with return code {e.returncode}"
        if e.stderr:
            error_msg += f"\nSTDERR: {e.stderr}"
        if e.stdout:
            error_msg += f"\nSTDOUT: {e.stdout}"
        raise GribToNetCDFError(error_msg) from e

    except Exception as e:
        raise GribToNetCDFError(f"Unexpected error during conversion: {str(e)}") from e


def convert_grib_to_nc_simple(
    grib_file: PathLike,
    output_file: PathLike,
    high_precision: bool = False,
    compress: bool = False,
) -> Path:
    """
    Simplified interface for GRIB to NetCDF conversion with common presets.

    Parameters
    ----------
    grib_file : PathLike
        Path to input GRIB file.
    output_file : PathLike
        Path to output NetCDF file.
    high_precision : bool, default False
        Use higher precision (NC_FLOAT instead of NC_SHORT).
    compress : bool, default False
        Enable compression for smaller file size (uses netCDF-4 format).

    Returns
    -------
    Path
        Path to the created NetCDF file.

    Examples
    --------
    Basic conversion:

    >>> convert_grib_to_nc_simple('data.grib', 'data.nc')

    High precision conversion:

    >>> convert_grib_to_nc_simple('data.grib', 'data.nc', high_precision=True)

    Compressed output:

    >>> convert_grib_to_nc_simple('data.grib', 'data.nc', compress=True)
    """

    # Set parameters based on presets
    data_type = "NC_FLOAT" if high_precision else "NC_SHORT"
    file_kind = 4 if compress else 2
    deflate_level = 6 if compress else None
    shuffle = compress

    return convert_grib_to_nc(
        grib_files=grib_file,
        output_file=output_file,
        data_type=data_type,
        file_kind=file_kind,
        deflate_level=deflate_level,
        shuffle=shuffle,
        unlimited_dimension="time",
    )


def batch_convert_grib_to_nc(
    input_directory: PathLike,
    output_directory: PathLike,
    pattern: str = "*.grib*",
    **kwargs,
) -> List[Path]:
    """
    Batch convert all GRIB files in a directory to NetCDF format.

    Parameters
    ----------
    input_directory : PathLike
        Directory containing GRIB files.
    output_directory : PathLike
        Directory to save converted NetCDF files.
    pattern : str, default "*.grib*"
        File pattern to match GRIB files.
    **kwargs
        Additional arguments passed to convert_grib_to_nc.

    Returns
    -------
    List[Path]
        List of paths to created NetCDF files.

    Examples
    --------
    Convert all GRIB files in a directory:

    >>> batch_convert_grib_to_nc('/path/to/grib_files', '/path/to/output')

    Convert with custom pattern and settings:

    >>> batch_convert_grib_to_nc(
    ...     '/path/to/grib_files',
    ...     '/path/to/output',
    ...     pattern="*.grb2",
    ...     high_precision=True
    ... )
    """
    input_dir = Path(input_directory)
    output_dir = Path(output_directory)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all GRIB files
    grib_files = list(input_dir.glob(pattern))

    if not grib_files:
        raise FileNotFoundError(
            f"No GRIB files found in {input_dir} with pattern {pattern}"
        )

    converted_files = []

    for grib_file in grib_files:
        # Create output filename
        output_file = output_dir / f"{grib_file.stem}.nc"

        try:
            result_path = convert_grib_to_nc(grib_file, output_file, **kwargs)
            converted_files.append(result_path)
        except Exception as e:
            logger.error(f"Failed to convert {grib_file}: {e}")
            continue

    print(f"✅ Batch conversion completed! Converted {len(converted_files)} files.")
    return converted_files


# Convenience function aliases
grib2nc = convert_grib_to_nc_simple  # Short alias
grib_to_netcdf = convert_grib_to_nc  # Alternative name
