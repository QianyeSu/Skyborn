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
    smooth_window: int = 1,
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

    if smooth_window > 1:
        half_window = smooth_window // 2
        growth_smoothed = np.full_like(growth, np.nan)
        for i in range(growth.size):
            start = max(0, i - half_window)
            stop = min(growth.size, i + half_window + 1)
            window_values = growth[start:stop]
            finite_mask = np.isfinite(window_values)
            if np.any(finite_mask):
                growth_smoothed[i] = np.mean(window_values[finite_mask])
        growth = growth_smoothed

    turning_points = np.where(growth[:-1] - growth[1:] > 0.0)[0]
    if turning_points.size == 0:
        return float("nan")

    return float(np.max(growth[: turning_points[-1] + 2]))


class TestGrowthRateHelpers:
    """Focused coverage for the growth-rate Python wrapper helpers."""

    def test_helper_missing_value_and_unit_normalization_branches(self):
        """Helper functions should normalize missing values and pressure units."""

        assert growth_rate_core._is_nan_missing_value(None) is True
        assert growth_rate_core._is_nan_missing_value(object()) is False
        assert_allclose(growth_rate_core._return_missing_value(-999.0), -999.0)
        assert (
            growth_rate_core._contains_missing(np.array([1.0, np.nan, 3.0]), np.nan)
            is True
        )
        assert (
            growth_rate_core._contains_missing(np.array([1.0, -999.0, 3.0]), -999.0)
            is True
        )
        assert_allclose(
            growth_rate_core._pressure_to_pa(np.array([1000.0, 850.0])),
            np.array([100000.0, 85000.0]),
        )
        assert_allclose(
            growth_rate_core._target_pressure_to_pa(np.array([300.0, 700.0])),
            np.array([30000.0, 70000.0]),
        )
        assert growth_rate_core._normalize_vertical_interp("pressure") == "linear"
        assert growth_rate_core._normalize_vertical_interp("logp") == "log"

    def test_helper_validation_errors(self):
        """Helper validation should reject malformed shapes, coordinates, and units."""

        with pytest.raises(ValueError, match="must be one-dimensional"):
            growth_rate_core._as_1d_float64(np.ones((2, 2)), "values")

        with pytest.raises(ValueError, match="dimension mismatch"):
            growth_rate_core._same_shape_or_raise(
                ("a", np.array([1.0, 2.0])),
                ("b", np.array([1.0])),
            )

        with pytest.raises(ValueError, match="strictly monotonic"):
            growth_rate_core._require_monotonic(
                np.array([1000.0, 900.0, 950.0]), "pressure"
            )

        with pytest.raises(ValueError, match="strictly positive"):
            growth_rate_core._pressure_to_pa(np.array([1000.0, 0.0]))

        with pytest.raises(ValueError, match="strictly positive"):
            growth_rate_core._target_pressure_to_pa(np.array([300.0, -1.0]))

        with pytest.raises(ValueError, match="output_units"):
            growth_rate_core._normalize_output_units("month^-1")

        with pytest.raises(ValueError, match="vertical_interp"):
            growth_rate_core._normalize_vertical_interp("cubic")

        with pytest.raises(TypeError, match="solver_levels"):
            growth_rate_core._normalize_solver_levels(3.5)

        with pytest.raises(ValueError, match="at least 2"):
            growth_rate_core._normalize_solver_levels(1)

        with pytest.raises(TypeError, match="smooth_window"):
            growth_rate_core._normalize_smooth_window(2.5)

        with pytest.raises(ValueError, match="at least 1"):
            growth_rate_core._normalize_smooth_window(0)

        with pytest.raises(ValueError, match="odd integer"):
            growth_rate_core._normalize_smooth_window(4)

        assert growth_rate_core._normalize_smooth_window(1) == 1
        assert growth_rate_core._normalize_smooth_window(np.int64(5)) == 5

        with pytest.raises(TypeError, match="smooth_window"):
            growth_rate_core._normalize_smooth_window(True)

    def test_infer_tropopause_reorders_descending_pressure_for_wmo_solver(
        self, monkeypatch
    ):
        """The WMO tropopause helper should reorder descending pressure inputs."""

        captured = {}

        def fake_tropopause(temperature_profile, pressure_profile, pressure_unit="Pa"):
            captured["temperature"] = np.array(temperature_profile, copy=True)
            captured["pressure"] = np.array(pressure_profile, copy=True)
            captured["pressure_unit"] = pressure_unit
            return {
                "pressure": 250.0,
                "height": 11000.0,
                "level_index": 2,
                "lapse_rate": 1.8,
                "success": True,
            }

        monkeypatch.setattr(growth_rate_core, "_trop_wmo_profile", fake_tropopause)

        result = growth_rate_core._infer_tropopause_pressure_pa(
            np.array([285.0, 260.0, 230.0]),
            np.array([100000.0, 70000.0, 30000.0]),
        )

        assert_allclose(result, 25000.0)
        assert captured["pressure_unit"] == "Pa"
        assert_allclose(captured["pressure"], np.array([30000.0, 70000.0, 100000.0]))
        assert_allclose(captured["temperature"], np.array([230.0, 260.0, 285.0]))

    def test_build_solver_pressure_grid_rejects_non_tropospheric_bounds(
        self, monkeypatch
    ):
        """The solver grid should reject diagnosed tropopauses at or below the lower bound."""

        monkeypatch.setattr(
            growth_rate_core,
            "_infer_tropopause_pressure_pa",
            lambda temperature, pressure_pa: 100000.0,
        )

        with pytest.raises(ValueError, match="diagnosed tropopause pressure"):
            growth_rate_core._build_solver_pressure_grid_pa(
                np.array([30000.0, 70000.0, 100000.0]),
                np.array([230.0, 260.0, 285.0]),
                target_pressure=None,
                solver_levels=45,
            )


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

    def test_barot_growth_rate_returns_missing_value_when_input_is_missing(self):
        """Missing barotropic inputs should short-circuit to the public missing marker."""

        result = barot_growth_rate(
            np.array([10.0, -999.0, 20.0]),
            np.array([20.0, 30.0, 40.0]),
            missing_value=-999.0,
        )

        assert result == -999.0

    def test_barot_growth_rate_reorders_descending_latitudes_for_backend(
        self, monkeypatch
    ):
        """Descending latitude input should be reordered before the backend call."""

        captured = {}

        def fake_barot_backend(lat_values, u_values):
            captured["lat"] = np.array(lat_values, dtype=np.float64, copy=True)
            captured["u"] = np.array(u_values, dtype=np.float64, copy=True)
            return 2.5e-6, 0

        monkeypatch.setattr(
            growth_rate_core, "_dbarot_growth_rate_1d", fake_barot_backend
        )

        result = barot_growth_rate(
            np.array([30.0, 20.0, 10.0]),
            np.array([50.0, 40.0, 30.0]),
            output_units="s^-1",
        )

        assert_allclose(result, 2.5e-6)
        assert_allclose(captured["lat"], np.array([30.0, 40.0, 50.0]))
        assert_allclose(captured["u"], np.array([10.0, 20.0, 30.0]))

    def test_barot_growth_rate_raises_when_backend_reports_error(self, monkeypatch):
        """Barotropic backend errors should raise a clear runtime error."""

        monkeypatch.setattr(
            growth_rate_core,
            "_dbarot_growth_rate_1d",
            lambda lat_values, u_values: (0.0, 17),
        )

        with pytest.raises(RuntimeError, match="ier=17"):
            barot_growth_rate(
                np.array([5.0, 12.0, 28.0]),
                np.array([30.0, 45.0, 60.0]),
            )

    def test_baroc_growth_rate_uses_wmo_tropopause_when_target_grid_is_omitted(
        self, monkeypatch
    ):
        """When no explicit solver grid is given, WMO tropopause should be used."""

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
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            lat,
            smooth_window,
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["temperature_solver"] = np.array(
                temperature_solver, dtype=np.float64, copy=True
            )
            captured["lat"] = float(lat)
            captured["smooth_window"] = int(smooth_window)
            return 1.25e-6, 0

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            fake_baroc_backend,
        )

        result = baroc_growth_rate(u, temperature, pressure, lat=45.0)

        assert_allclose(result, 1.25e-6, rtol=0.0, atol=0.0)
        assert captured["lat"] == 45.0
        assert captured["smooth_window"] == 1
        assert_allclose(captured["pressure_solver"][0], 23000.0, rtol=0.0, atol=1e-12)
        assert_allclose(captured["pressure_solver"][-1], 100000.0, rtol=0.0, atol=1e-12)
        assert captured["temperature_solver"].shape == captured["pressure_solver"].shape

    def test_baroc_growth_rate_exposes_solver_levels_for_auto_grid(self, monkeypatch):
        """The public API should let callers control the automatic solver-grid size."""

        captured = {}

        monkeypatch.setattr(
            growth_rate_core,
            "_infer_tropopause_pressure_pa",
            lambda temperature, pressure_pa: 25000.0,
        )

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            lat,
            smooth_window,
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["smooth_window"] = int(smooth_window)
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0, 32.0, 38.0]),
            np.array([220.0, 235.0, 255.0, 272.0, 288.0]),
            np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0]),
            lat=45.0,
            solver_levels=11,
        )

        assert_allclose(result, 1.0e-6, rtol=0.0, atol=0.0)
        assert captured["pressure_solver"].shape == (11,)
        assert captured["smooth_window"] == 1
        assert_allclose(captured["pressure_solver"][0], 25000.0, atol=1e-12)
        assert_allclose(captured["pressure_solver"][-1], 100000.0, atol=1e-12)

    def test_baroc_growth_rate_ignores_solver_levels_when_target_grid_is_explicit(
        self, monkeypatch
    ):
        """Explicit target_pressure should override the automatic solver-level setting."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            lat,
            smooth_window,
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["smooth_window"] = int(smooth_window)
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([30000.0, 55000.0, 80000.0, 100000.0]),
            solver_levels=99,
        )

        assert_allclose(result, 1.0e-6, rtol=0.0, atol=0.0)
        assert captured["smooth_window"] == 1
        assert_allclose(
            captured["pressure_solver"],
            np.array([30000.0, 55000.0, 80000.0, 100000.0]),
        )

    def test_baroc_growth_rate_passes_smooth_window_to_fortran_backend(
        self, monkeypatch
    ):
        """The public API should forward smooth_window to the compiled backend."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            lat,
            smooth_window,
        ):
            captured["smooth_window"] = int(smooth_window)
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([30000.0, 60000.0, 100000.0]),
            smooth_window=5,
        )

        assert_allclose(result, 1.0e-6, rtol=0.0, atol=0.0)
        assert captured["smooth_window"] == 5

    def test_baroc_growth_rate_matches_reference_with_spectrum_smoothing(self):
        """Compiled spectrum smoothing should match the NumPy reference path."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0])
        u = np.array([8.0, 14.0, 24.0, 32.0, 38.0])
        lat = 45.0

        expected = _baroc_growth_reference(
            u,
            temperature,
            pressure,
            lat,
            smooth_window=5,
        )
        result = baroc_growth_rate(
            u,
            temperature,
            pressure,
            lat=lat,
            target_pressure=pressure,
            smooth_window=5,
            output_units="s^-1",
        )

        assert_allclose(result, expected, rtol=1e-8, atol=1e-11)

    def test_baroc_growth_rate_defaults_to_log_pressure_interpolation(
        self, monkeypatch
    ):
        """The default baroclinic interpolation mode should be log-pressure."""

        methods = []

        monkeypatch.setattr(
            growth_rate_core,
            "_infer_tropopause_pressure_pa",
            lambda temperature, pressure_pa: 30000.0,
        )

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            methods.append(method)
            return np.asarray(x, dtype=np.float64)

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, lat, smooth_window: (
                1.0e-6,
                0,
            ),
        )

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([230.0, 250.0, 280.0]),
            np.array([30000.0, 65000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([30000.0, 65000.0, 100000.0]),
        )

        assert_allclose(result, 1.0e-6)
        assert methods == ["log", "log"]

    def test_baroc_growth_rate_accepts_logp_alias(self, monkeypatch):
        """The legacy ``logp`` alias should still map to log-pressure interpolation."""

        methods = []

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            methods.append(method)
            return np.asarray(x, dtype=np.float64)

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, lat, smooth_window: (
                1.0e-6,
                0,
            ),
        )

        baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([230.0, 250.0, 280.0]),
            np.array([30000.0, 65000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([30000.0, 65000.0, 100000.0]),
            vertical_interp="logp",
        )

        assert methods == ["log", "log"]

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

    def test_baroc_growth_rate_returns_missing_value_when_input_is_missing(self):
        """Missing baroclinic inputs should short-circuit to the public missing marker."""

        result = baroc_growth_rate(
            np.array([8.0, -999.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([30000.0, 60000.0, 100000.0]),
            missing_value=-999.0,
        )

        assert result == -999.0

    def test_baroc_growth_rate_rejects_nonpositive_temperature(self):
        """Baroclinic temperature profiles must remain strictly positive."""

        with pytest.raises(ValueError, match="strictly positive"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 0.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                target_pressure=np.array([30000.0, 60000.0, 100000.0]),
            )

    def test_baroc_growth_rate_returns_missing_when_interpolation_has_nan(
        self, monkeypatch
    ):
        """Interpolation failures should return the public missing marker."""

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            if float(np.asarray(x)[0]) == 8.0:
                return np.array([np.nan, np.nan, np.nan])
            return np.array([220.0, 235.0, 255.0])

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([30000.0, 60000.0, 100000.0]),
            missing_value=-999.0,
        )

        assert result == -999.0

    def test_baroc_growth_rate_reorders_descending_target_pressure_for_backend(
        self, monkeypatch
    ):
        """Descending explicit solver grids should be reordered before the backend call."""

        captured = {}

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            first_value = float(np.asarray(x)[0])
            if first_value == 8.0:
                return np.array([101.0, 102.0, 103.0])
            return np.array([201.0, 202.0, 203.0])

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            lat,
            smooth_window,
        ):
            captured["u_solver"] = np.array(u_solver, dtype=np.float64, copy=True)
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["temperature_solver"] = np.array(
                temperature_solver, dtype=np.float64, copy=True
            )
            captured["smooth_window"] = int(smooth_window)
            return 2.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)
        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            target_pressure=np.array([100000.0, 60000.0, 30000.0]),
        )

        assert_allclose(result, 2.0e-6)
        assert_allclose(
            captured["pressure_solver"], np.array([30000.0, 60000.0, 100000.0])
        )
        assert captured["smooth_window"] == 1
        assert_allclose(captured["u_solver"], np.array([103.0, 102.0, 101.0]))
        assert_allclose(captured["temperature_solver"], np.array([203.0, 202.0, 201.0]))

    def test_baroc_growth_rate_raises_when_backend_reports_error(self, monkeypatch):
        """Baroclinic backend errors should raise a clear runtime error."""

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, lat, smooth_window: (
                0.0,
                19,
            ),
        )

        with pytest.raises(RuntimeError, match="ier=19"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                target_pressure=np.array([30000.0, 60000.0, 100000.0]),
            )

    def test_baroc_growth_rate_rejects_even_smooth_window(self):
        """Centered spectrum smoothing should require an odd window length."""

        with pytest.raises(ValueError, match="odd integer"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                target_pressure=np.array([30000.0, 60000.0, 100000.0]),
                smooth_window=4,
            )
