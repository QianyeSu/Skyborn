"""Focused tests for Chemke-style growth-rate diagnostics."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

import skyborn.calc.growth_rate.core as growth_rate_core
from skyborn.calc.growth_rate import baroc_growth_rate, barot_growth_rate

OMEGA = 7.292e-5
RADIUS = 6_371_000.0
GAS_CONSTANT_DRY = 287.04
HEAT_CAPACITY = 1004.7
KAPPA = GAS_CONSTANT_DRY / HEAT_CAPACITY


def _theta_from_temperature(
    temperature: np.ndarray, pressure_pa: np.ndarray
) -> np.ndarray:
    """Return potential temperature from temperature and pressure."""

    return (
        np.asarray(temperature, dtype=np.float64)
        * (100000.0 / np.asarray(pressure_pa, dtype=np.float64)) ** KAPPA
    )


def _barot_growth_reference(lat: np.ndarray, u: np.ndarray) -> float:
    """Pure NumPy port of the Chemke barotropic reference algorithm."""

    lat = np.asarray(lat, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    y = np.deg2rad(lat) * RADIUS
    dy1 = (y[1:] + y[:-1]) / 2.0
    dy2 = dy1[1:] - dy1[:-1]
    dy = (y[1:] - y[:-1]) * np.cos(np.deg2rad((lat[1:] + lat[:-1]) / 2.0))
    y1 = (y[1] - y[0]) * np.cos(np.deg2rad((lat[1] + lat[0]) / 2.0))
    y2 = (y[-1] - y[-2]) * np.cos(np.deg2rad((lat[-1] + lat[-2]) / 2.0))
    beta = 2.0 * OMEGA * np.cos(np.deg2rad(lat)) / RADIUS
    waven = np.arange(0.0, 2.0e-5, 1.0e-7, dtype=np.float64)
    growth = np.full(waven.shape, np.nan, dtype=np.float64)

    for k, kval in enumerate(waven):
        a = np.zeros((lat.size, lat.size), dtype=np.float64)
        b = np.zeros((lat.size, lat.size), dtype=np.float64)
        for n in range(lat.size):
            if n == lat.size - 1:
                a[n, n] = (
                    kval * beta[n]
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) / y2 * (-1.0 / dy[n - 1] - 1.0 / y2)
                    - kval / y2 * ((u[n - 1] - u[n]) / dy[n - 1] - u[n] / y2)
                )
                a[n, n - 1] = (kval * u[n]) / (y2 * dy[n - 1])
                b[n, n] = -(kval**2) + (1.0 / y2) * (-1.0 / dy[n - 1] - 1.0 / y2)
                b[n, n - 1] = 1.0 / (y2 * dy[n - 1])
            elif n == 0:
                a[n, n] = (
                    kval * beta[n]
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) / y1 * (-1.0 / y1 - 1.0 / dy[n])
                    - kval / y1 * ((-u[n]) / y1 - (u[n] - u[n + 1]) / dy[n])
                )
                a[n, n + 1] = (kval * u[n]) / (y1 * dy[n])
                b[n, n] = -(kval**2) + (1.0 / y1) * (-1.0 / y1 - 1.0 / dy[n])
                b[n, n + 1] = 1.0 / (y1 * dy[n])
            else:
                a[n, n - 1] = (kval * u[n]) / (dy2[n - 1] * dy[n - 1])
                a[n, n + 1] = (kval * u[n]) / (dy2[n - 1] * dy[n])
                a[n, n] = (
                    kval * beta[n]
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) / dy2[n - 1] * (-1.0 / dy[n - 1] - 1.0 / dy[n])
                    - kval
                    / dy2[n - 1]
                    * ((u[n - 1] - u[n]) / dy[n - 1] - (u[n] - u[n + 1]) / dy[n])
                )
                b[n, n - 1] = 1.0 / (dy2[n - 1] * dy[n - 1])
                b[n, n + 1] = 1.0 / (dy2[n - 1] * dy[n])
                b[n, n] = -(kval**2) + (1.0 / dy2[n - 1]) * (
                    -1.0 / dy[n - 1] - 1.0 / dy[n]
                )

        eigvals = np.linalg.eigvals(np.linalg.solve(b, a))
        growth[k] = np.max(np.imag(eigvals))

    return float(np.nanmax(growth))


def _baroc_growth_reference(
    u: np.ndarray,
    temperature: np.ndarray,
    pressure_pa: np.ndarray,
    lat: float,
) -> float:
    """Pure NumPy port of the Chemke baroclinic reference algorithm."""

    u = np.asarray(u, dtype=np.float64)
    temperature = np.asarray(temperature, dtype=np.float64)
    pressure_pa = np.asarray(pressure_pa, dtype=np.float64)
    theta = _theta_from_temperature(temperature, pressure_pa)

    dz1 = (pressure_pa[1:] + pressure_pa[:-1]) / 2.0
    dz2 = np.empty_like(pressure_pa)
    dz2[0] = pressure_pa[1] - pressure_pa[0]
    dz2[1:-1] = dz1[1:] - dz1[:-1]
    dz2[-1] = pressure_pa[-1] - pressure_pa[-2]

    f = 2.0 * OMEGA * np.sin(np.deg2rad(lat))
    beta = 2.0 * OMEGA * np.cos(np.deg2rad(lat)) / RADIUS
    rho = pressure_pa / (GAS_CONSTANT_DRY * temperature)
    nz = -((theta[1:] - theta[:-1]) * (1.0 / (rho[:-1] * theta[:-1])))
    waven = np.arange(0.0, 2.0e-5, 1.0e-7, dtype=np.float64)
    growth = np.full(waven.shape, np.nan, dtype=np.float64)

    for k, kval in enumerate(waven):
        a = np.zeros((pressure_pa.size, pressure_pa.size), dtype=np.float64)
        b = np.zeros((pressure_pa.size, pressure_pa.size), dtype=np.float64)
        for n in range(pressure_pa.size):
            if n == pressure_pa.size - 1:
                a[n, n] = (
                    kval * beta
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) * (f**2) / dz2[n] * (-1.0 / nz[n - 1])
                    - kval * (f**2) / dz2[n] * ((u[n - 1] - u[n]) / nz[n - 1])
                )
                a[n, n - 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n - 1])
                b[n, n] = -(kval**2) + (f**2) / dz2[n] * (-1.0 / nz[n - 1])
                b[n, n - 1] = (f**2) / (dz2[n] * nz[n - 1])
            elif n == 0:
                a[n, n] = (
                    kval * beta
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) * (f**2) / dz2[n] * (-1.0 / nz[n])
                    - kval * (f**2) / dz2[n] * ((u[n + 1] - u[n]) / nz[n])
                )
                a[n, n + 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n])
                b[n, n] = -(kval**2) + (f**2) / dz2[n] * (-1.0 / nz[n])
                b[n, n + 1] = (f**2) / (dz2[n] * nz[n])
            else:
                a[n, n - 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n - 1])
                a[n, n + 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n])
                a[n, n] = (
                    kval * beta
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) * (f**2) / dz2[n] * (-1.0 / nz[n - 1] - 1.0 / nz[n])
                    - kval
                    * (f**2)
                    / dz2[n]
                    * ((u[n - 1] - u[n]) / nz[n - 1] - (u[n] - u[n + 1]) / nz[n])
                )
                b[n, n - 1] = (f**2) / (dz2[n] * nz[n - 1])
                b[n, n + 1] = (f**2) / (dz2[n] * nz[n])
                b[n, n] = -(kval**2) + (f**2) / dz2[n] * (
                    -1.0 / nz[n - 1] - 1.0 / nz[n]
                )

        eigvals = np.linalg.eigvals(np.linalg.solve(b, a))
        growth[k] = np.max(np.imag(eigvals))

    turning_points = np.where(growth[:-1] - growth[1:] > 0.0)[0]
    if turning_points.size == 0:
        return float("nan")

    return float(np.max(growth[: turning_points[-1] + 2]))


class TestGrowthRate:
    """Regression checks for the first-stage growth-rate API."""

    def test_barot_growth_rate_matches_reference(self):
        """Compiled barotropic growth should match the Chemke-style NumPy port."""

        lat = np.linspace(30.0, 60.0, 7)
        u = np.array([5.0, 12.0, 28.0, 42.0, 26.0, 11.0, 4.0], dtype=np.float64)

        expected = _barot_growth_reference(lat, u)
        result = barot_growth_rate(u, lat, output_units="s^-1")

        assert_allclose(result, expected, rtol=1e-8, atol=1e-11)

    def test_barot_growth_rate_constant_flow_is_nearly_neutral(self):
        """A constant barotropic flow should not produce meaningful growth."""

        lat = np.linspace(25.0, 55.0, 7)
        u = np.full(lat.shape, 20.0)

        result = barot_growth_rate(u, lat, output_units="s^-1")

        assert np.isfinite(result)
        assert abs(result) < 1e-10

    def test_baroc_growth_rate_matches_reference_on_explicit_solver_grid(self):
        """Compiled baroclinic growth should match the Chemke-style NumPy port."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0])
        u = np.array([8.0, 14.0, 24.0, 32.0, 38.0])
        lat = 45.0

        expected = _baroc_growth_reference(u, temperature, pressure, lat)
        result = baroc_growth_rate(
            u,
            temperature,
            pressure,
            lat=lat,
            tropopause_pressure=pressure[0],
            target_pressure=pressure,
            output_units="s^-1",
        )

        assert_allclose(result, expected, rtol=1e-8, atol=1e-11)

    def test_baroc_growth_rate_requires_latitude(self):
        """Latitude is required for the Coriolis-dependent baroclinic solver."""

        with pytest.raises(ValueError, match="`lat` is required"):
            baroc_growth_rate(
                np.array([10.0, 20.0, 30.0]),
                np.array([240.0, 260.0, 280.0]),
                np.array([30000.0, 60000.0, 100000.0]),
            )

    def test_output_units_conversion_to_day_inverse(self):
        """The public API should convert from per-second to per-day units explicitly."""

        lat = np.linspace(30.0, 60.0, 7)
        u = np.array([5.0, 12.0, 28.0, 42.0, 26.0, 11.0, 4.0], dtype=np.float64)

        growth_s = barot_growth_rate(u, lat, output_units="s^-1")
        growth_day = barot_growth_rate(u, lat, output_units="day^-1")

        assert_allclose(growth_day, growth_s * 86400.0, rtol=1e-12, atol=1e-12)

    def test_barot_growth_rate_defaults_to_per_second_units(self):
        """The default public output unit should now be per-second growth."""

        lat = np.linspace(30.0, 60.0, 7)
        u = np.array([5.0, 12.0, 28.0, 42.0, 26.0, 11.0, 4.0], dtype=np.float64)

        default_result = barot_growth_rate(u, lat)
        explicit_result = barot_growth_rate(u, lat, output_units="s^-1")

        assert_allclose(default_result, explicit_result, rtol=1e-12, atol=1e-12)

    def test_baroc_growth_rate_infers_tropopause_pressure_when_omitted(
        self, monkeypatch
    ):
        """When no explicit tropopause pressure is given, WMO tropopause should be used."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0])
        u = np.array([8.0, 14.0, 24.0, 32.0, 38.0])
        captured = {}

        monkeypatch.setattr(
            growth_rate_core,
            "_trop_wmo_profile",
            lambda temperature_profile, pressure_profile, pressure_unit="Pa": {
                "pressure": 230.0,
                "height": 11000.0,
                "level_index": 1,
                "lapse_rate": 1.8,
                "success": True,
            },
        )

        def fake_baroc_backend(
            u_solver, theta_solver, pressure_solver, temperature_solver, lat
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["temperature_solver"] = np.array(
                temperature_solver, dtype=np.float64, copy=True
            )
            captured["lat"] = float(lat)
            return 1.25e-6, 0

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            fake_baroc_backend,
        )

        result = baroc_growth_rate(u, temperature, pressure, lat=45.0)

        assert_allclose(result, 1.25e-6, rtol=0.0, atol=0.0)
        assert captured["lat"] == 45.0
        assert_allclose(captured["pressure_solver"][0], 23000.0, rtol=0.0, atol=1e-12)
        assert_allclose(captured["pressure_solver"][-1], 100000.0, rtol=0.0, atol=1e-12)
        assert captured["temperature_solver"].shape == captured["pressure_solver"].shape

    def test_baroc_growth_rate_raises_when_auto_tropopause_detection_fails(
        self, monkeypatch
    ):
        """A failed WMO tropopause diagnosis should raise a clear error."""

        monkeypatch.setattr(
            growth_rate_core,
            "_trop_wmo_profile",
            lambda temperature_profile, pressure_profile, pressure_unit="Pa": {
                "pressure": -999.0,
                "height": -999.0,
                "level_index": -1,
                "lapse_rate": -999.0,
                "success": False,
            },
        )

        with pytest.raises(RuntimeError, match="Automatic WMO tropopause detection"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0, 32.0, 38.0]),
                np.array([220.0, 235.0, 255.0, 272.0, 288.0]),
                np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0]),
                lat=45.0,
            )
