"""
Shared test fixtures and configuration for skyborn test suite.

This file provides common fixtures and pytest configuration for testing
the skyborn library across all modules.
"""

import pytest
import numpy as np
import xarray as xr
import pandas as pd
import tempfile
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for testing
import matplotlib.pyplot as plt
from pathlib import Path

# Seed random number generator for reproducible tests
np.random.seed(42)


@pytest.fixture
def sample_climate_data():
    """
    Create standard climate test data as xarray Dataset.

    Returns
    -------
    xr.Dataset
        Dataset with temperature and precipitation data
    """
    time = pd.date_range("2020-01-01", periods=12, freq="M")
    lat = np.linspace(-90, 90, 73)
    lon = np.linspace(0, 357.5, 144)

    # Create realistic temperature data (seasonal cycle + noise)
    temp_base = 273.15 + 15 * np.cos(2 * np.pi * np.arange(12) / 12)  # Seasonal cycle
    temp_spatial = 20 * np.cos(np.deg2rad(lat))[:, np.newaxis]  # Latitude gradient
    temp_data = (
        temp_base[:, np.newaxis, np.newaxis]
        + temp_spatial[np.newaxis, :, :]
        + np.random.randn(12, 73, 144) * 2
    )

    # Create precipitation data (positive values only)
    precip_data = np.random.exponential(2, (12, 73, 144))

    return xr.Dataset(
        {
            "temperature": (
                ["time", "lat", "lon"],
                temp_data,
                {"units": "K", "long_name": "Surface Air Temperature"},
            ),
            "precipitation": (
                ["time", "lat", "lon"],
                precip_data,
                {"units": "mm/day", "long_name": "Daily Precipitation Rate"},
            ),
        },
        coords={"time": time, "lat": lat, "lon": lon},
    )


@pytest.fixture
def sample_2d_field():
    """
    Create a simple 2D field for testing plotting functions.

    Returns
    -------
    xr.DataArray
        2D temperature field
    """
    lat = np.linspace(-90, 90, 37)
    lon = np.linspace(0, 357.5, 72)

    # Create a realistic temperature pattern
    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing="ij")
    temp = (
        273.15
        + 25 * np.cos(np.deg2rad(lat_grid))
        + 5 * np.sin(2 * np.deg2rad(lon_grid))
        + np.random.randn(37, 72) * 2
    )

    return xr.DataArray(
        temp,
        coords={"lat": lat, "lon": lon},
        dims=["lat", "lon"],
        attrs={"units": "K", "long_name": "Surface Air Temperature"},
    )


@pytest.fixture
def sample_wind_data():
    """
    Create sample wind field data for vector analysis.

    Returns
    -------
    tuple
        Tuple of (u, v) wind components as xarray DataArrays
    """
    lat = np.linspace(-90, 90, 37)
    lon = np.linspace(0, 357.5, 72)

    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing="ij")

    # Create realistic wind patterns
    u = 20 * np.cos(np.deg2rad(lat_grid)) + np.random.randn(37, 72) * 3
    v = 5 * np.sin(2 * np.deg2rad(lon_grid)) + np.random.randn(37, 72) * 2

    u_da = xr.DataArray(
        u,
        coords={"lat": lat, "lon": lon},
        dims=["lat", "lon"],
        attrs={"units": "m/s", "long_name": "Zonal Wind Component"},
    )

    v_da = xr.DataArray(
        v,
        coords={"lat": lat, "lon": lon},
        dims=["lat", "lon"],
        attrs={"units": "m/s", "long_name": "Meridional Wind Component"},
    )

    return u_da, v_da


@pytest.fixture
def sample_grid():
    """
    Create standard latitude-longitude grid coordinates.

    Returns
    -------
    dict
        Dictionary with 'lat' and 'lon' coordinate arrays
    """
    return {"lat": np.linspace(-90, 90, 73), "lon": np.linspace(0, 357.5, 144)}


@pytest.fixture
def temp_dir():
    """
    Create a temporary directory for file I/O tests.

    Yields
    ------
    pathlib.Path
        Path to temporary directory
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_1d_timeseries():
    """
    Create a 1D time series for line plotting tests.

    Returns
    -------
    xr.DataArray
        1D time series data
    """
    time = pd.date_range("2020-01-01", periods=100, freq="D")
    # Create time series with trend and noise
    trend = 0.1 * np.arange(100)
    seasonal = 2 * np.sin(2 * np.pi * np.arange(100) / 30)
    noise = np.random.randn(100) * 0.5
    data = trend + seasonal + noise

    return xr.DataArray(
        data,
        coords={"time": time},
        dims=["time"],
        attrs={"units": "units", "long_name": "Sample Time Series"},
    )


@pytest.fixture
def sample_regression_data():
    """
    Create sample data for regression testing.

    Returns
    -------
    tuple
        Tuple of (data, predictor) arrays
    """
    # Create 3D data array (time, lat, lon)
    n_time, n_lat, n_lon = 50, 20, 30
    predictor = np.random.randn(n_time)

    # Create data with known relationship to predictor
    slope_field = np.random.randn(n_lat, n_lon)
    intercept_field = np.random.randn(n_lat, n_lon)

    data = np.zeros((n_time, n_lat, n_lon))
    for i in range(n_time):
        data[i] = (
            slope_field * predictor[i]
            + intercept_field
            + np.random.randn(n_lat, n_lon) * 0.1
        )

    return data, predictor


# Configure matplotlib for testing
@pytest.fixture(autouse=True)
def setup_matplotlib():
    """Automatically configure matplotlib for all tests."""
    plt.ioff()  # Turn off interactive mode
    yield
    plt.close("all")  # Close all figures after each test


# pytest markers
def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "mpl_image_compare: marks tests as matplotlib image comparison tests"
    )
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line(
        "markers", "requires_data: marks tests that require external data files"
    )


# Skip tests based on available dependencies
def pytest_collection_modifyitems(config, items):
    """Modify test collection to add skip markers based on dependencies."""
    try:
        import iris
    except ImportError:
        skip_iris = pytest.mark.skip(reason="iris not available")
        for item in items:
            if "iris" in str(item.fspath):
                item.add_marker(skip_iris)

    try:
        import eccodes
    except ImportError:
        skip_eccodes = pytest.mark.skip(reason="eccodes not available")
        for item in items:
            if "grib" in str(item.fspath).lower():
                item.add_marker(skip_eccodes)
