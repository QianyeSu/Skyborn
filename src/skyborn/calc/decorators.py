"""
Decorators for the Skyborn calc module.

This module provides decorators for data validation, caching, and file handling
in statistical calculations.
"""

import functools
import numpy as np
import os
import inspect
import warnings
from typing import Any, Dict, Callable, Optional, Union

__all__ = [
    "validate_dimensions",
    "cache_result",
    "handle_missing_data",
    "dims_test",
    "check_if_file_exists",
]

try:
    import xarray as xr

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False
    xr = None


def validate_dimensions(
    required_dims: Optional[list] = None, max_dims: int = 3
) -> Callable:
    """
    Decorator to validate xarray DataArray dimensions.

    Parameters
    ----------
    required_dims : list, optional
        List of required dimension names
    max_dims : int, default 3
        Maximum number of dimensions allowed

    Returns
    -------
    decorator : callable
        Decorated function
    """

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            if not hasattr(self, "DataArray"):
                return func(self, *args, **kwargs)

            data = self.DataArray
            dims = list(data.dims)

            # Check maximum dimensions
            if len(dims) > max_dims:
                raise ValueError(
                    f"Function {func.__name__} supports maximum {max_dims} dimensions, "
                    f"got {len(dims)}: {dims}"
                )

            # Check minimum dimensions
            if len(dims) < 2:
                raise ValueError(
                    f"Function {func.__name__} requires at least 2 dimensions, "
                    f"got {len(dims)}: {dims}"
                )

            # Check required dimensions
            if required_dims:
                # Handle coordinate name mapping if provided
                if hasattr(self, "coords_name") and self.coords_name:
                    coords_name = self.coords_name
                    rename_dict = {}

                    # Create rename dictionary
                    for key, value in coords_name.items():
                        if key in required_dims and value in dims:
                            rename_dict[value] = key
                        elif value in required_dims and key in dims:
                            rename_dict[key] = value

                    if rename_dict:
                        self.DataArray = data.rename(rename_dict)
                        dims = list(self.DataArray.dims)

                # Check if required dimensions are present
                missing_dims = set(required_dims) - set(dims)
                if missing_dims:
                    raise ValueError(
                        f"Missing required dimensions: {missing_dims}. "
                        f"Available dimensions: {dims}. "
                        f"Use 'coords_name' parameter to map coordinate names."
                    )

            # Store ordered dimensions for later use
            if hasattr(self, "DataArray"):
                self.ordered_dims = list(self.DataArray.dims)

            return func(self, *args, **kwargs)

        return wrapper

    return decorator


def cache_result(
    cache_dir: str = ".cache", file_prefix: str = "result", recompute: bool = False
) -> Callable:
    """
    Decorator to cache computation results to disk.

    Parameters
    ----------
    cache_dir : str, default '.cache'
        Directory to store cached files
    file_prefix : str, default 'result'
        Prefix for cached filenames
    recompute : bool, default False
        Force recomputation even if cache exists

    Returns
    -------
    decorator : callable
        Decorated function
    """

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Extract path from kwargs if provided
            save_path = kwargs.get("path", None)
            force_recompute = kwargs.get("rewrite_trends", recompute)

            # If path is provided and file exists and not forcing recompute
            if save_path and os.path.exists(save_path) and not force_recompute:
                warnings.warn(
                    f"Loading existing results from {save_path}. "
                    f"Set rewrite_trends=False to force recomputation."
                )
                if HAS_XARRAY:
                    return xr.open_dataset(save_path)
                else:
                    raise ImportError("xarray required for loading cached results")

            # Compute result
            result = func(*args, **kwargs)

            # Save result if path is provided
            if save_path and hasattr(result, "to_netcdf"):
                os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
                result.to_netcdf(save_path)

            return result

        return wrapper

    return decorator


def handle_missing_data(fill_method: str = "drop", min_valid: int = 3) -> Callable:
    """
    Decorator to handle missing data in time series analysis.

    Parameters
    ----------
    fill_method : str, default 'drop'
        Method to handle missing data ('drop', 'interpolate', 'mean')
    min_valid : int, default 3
        Minimum number of valid observations required

    Returns
    -------
    decorator : callable
        Decorated function
    """

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(data, *args, **kwargs):
            if hasattr(data, "values"):
                # xarray DataArray
                data_values = data.values
                if hasattr(data, "dims"):
                    time_dim = kwargs.get("dim", "time")
                    if time_dim in data.dims:
                        # Handle missing data along time dimension
                        valid_count = np.sum(
                            np.isfinite(data_values), axis=data.get_axis_num(time_dim)
                        )
                        if np.any(valid_count < min_valid):
                            warnings.warn(
                                f"Some locations have fewer than {min_valid} "
                                f"valid observations"
                            )
            else:
                # numpy array
                data_values = np.asarray(data)
                valid_count = np.sum(np.isfinite(data_values), axis=0)
                if np.any(valid_count < min_valid):
                    warnings.warn(
                        f"Some locations have fewer than {min_valid} "
                        f"valid observations"
                    )

            return func(data, *args, **kwargs)

        return wrapper

    return decorator


# Legacy compatibility decorators
def dims_test(func: Callable) -> Callable:
    """
    Legacy compatibility decorator for dimension testing.

    This is a simplified version of the original dims_test decorator
    that provides basic dimension validation for backward compatibility.
    """

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if hasattr(self, "DataArray"):
            dims = self.DataArray.dims

            if len(dims) == 1:
                raise ValueError("requires at least a 2D dataarray (x,t)")
            elif len(dims) > 3:
                raise ValueError(
                    "Currently only supports 2D (x,t) and 3D dataarray (x,y,t)"
                )

            # Store ordered dims for compatibility
            self.ordered_dims = np.flipud(sorted(dims))

        return func(self, *args, **kwargs)

    return wrapper


def check_if_file_exists(func: Callable) -> Callable:
    """
    Legacy compatibility decorator for file existence checking.

    This provides backward compatibility with the original decorator.
    """

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        # Get function defaults
        signature = inspect.signature(func)
        defaults = {
            k: v.default
            for k, v in signature.parameters.items()
            if v.default is not inspect.Parameter.empty
        }

        # Update with provided kwargs
        params = defaults.copy()
        params.update(kwargs)

        # Check if should reload existing file
        if (
            params.get("rewrite_trends", True)
            and params.get("path") is not None
            and os.path.exists(params["path"])
        ):

            warnings.warn(
                "Loading file. To rewrite, set rewrite_trends=False or change filename."
            )
            if HAS_XARRAY:
                return xr.open_dataset(params["path"])
            else:
                raise ImportError("xarray required for loading files")

        return func(self, *args, **kwargs)

    return wrapper


# Modern decorator aliases for new code
validate_xarray_dims = validate_dimensions
cache_computation = cache_result
validate_time_series = handle_missing_data
