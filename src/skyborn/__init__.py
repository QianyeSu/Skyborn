# skyborn/__init__.py

# Version info
__version__ = "0.3.11"

# Define what should be available at top level for lazy import
_CALC_FUNCTIONS = {
    "linear_regression",
    "convert_longitude_range",
    "pearson_correlation",
    "spearman_correlation",
    "kendall_correlation",
    "calculate_potential_temperature",
    "mann_kendall_test",
    "mann_kendall_multidim",
    "mann_kendall_xarray",
    "trend_analysis",
    "mk_test",
    "mk_multidim",
    "gaussian_pdf",
    "emergent_constraint_posterior",
    "emergent_constraint_prior",
    "calc_GAUSSIAN_PDF",
    "calc_PDF_EC",
    "find_std_from_PDF",
    "calc_PDF_EC_PRIOR",
}

_GRADIENT_FUNCTIONS = {
    "calculate_gradient",
    "calculate_meridional_gradient",
    "calculate_zonal_gradient",
    "calculate_vertical_gradient",
}

_CAUSALITY_FUNCTIONS = {"liang_causality", "granger_causality"}

_CONVERSION_FUNCTIONS = {
    "convert_grib_to_nc",
    "convert_grib_to_nc_simple",
    "batch_convert_grib_to_nc",
    "grib2nc",
    "grib_to_netcdf",
    "GribToNetCDFError",
}

_SUBMODULES = {
    "plot",
    "interp",
    "ROF",
    "conversion",
    "calc",
    "spharm",
    "windspharm",
    "gridfill",
}

# Cache for already imported modules/functions
_imported_cache = {}


def __getattr__(name):
    """
    Lazy import mechanism to avoid loading heavy dependencies at startup.

    This significantly improves import time by only loading modules when they're actually used.
    """
    # Check cache first
    if name in _imported_cache:
        return _imported_cache[name]

    # Handle calculation functions
    if name in _CALC_FUNCTIONS:
        from . import calc

        func = getattr(calc, name)
        _imported_cache[name] = func
        return func

    # Handle gradient functions
    elif name in _GRADIENT_FUNCTIONS:
        from . import gradients

        func = getattr(gradients, name)
        _imported_cache[name] = func
        return func

    # Handle causality functions
    elif name in _CAUSALITY_FUNCTIONS:
        from . import causality

        func = getattr(causality, name)
        _imported_cache[name] = func
        return func

    # Handle conversion functions
    elif name in _CONVERSION_FUNCTIONS:
        from . import conversion

        func = getattr(conversion, name)
        _imported_cache[name] = func
        return func

    # Handle submodules
    elif name in _SUBMODULES:
        import importlib

        module = importlib.import_module(f".{name}", package="skyborn")
        _imported_cache[name] = module
        return module

    # Handle special cases
    elif name == "fill":
        from .gridfill import fill

        _imported_cache[name] = fill
        return fill
    elif name == "fill_cube":
        from .gridfill import fill_cube

        _imported_cache[name] = fill_cube
        return fill_cube
    elif name == "gridfill_fill":
        from .gridfill import fill

        _imported_cache[name] = fill
        return fill
    elif name == "gridfill_fill_cube":
        from .gridfill import fill_cube

        _imported_cache[name] = fill_cube
        return fill_cube

    # Not found
    raise AttributeError(f"module 'skyborn' has no attribute '{name}'")
