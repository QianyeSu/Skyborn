"""Regression tests for GPI ``error_flag`` initialization and aggregation."""

import numpy as np
import pytest

from skyborn.calc.GPI.core import potential_intensity

PRESSURE_LEVELS = np.array([1000.0, 850.0, 700.0, 500.0, 300.0, 200.0, 100.0])
TEMPERATURE_PROFILE = np.array([300.0, 290.0, 280.0, 265.0, 240.0, 220.0, 200.0])
MIXING_RATIO_PROFILE = np.array([0.018, 0.012, 0.007, 0.002, 2e-4, 2e-5, 1e-6])


def _inputs(sst_value, shape=(4, 5)):
    nlev = PRESSURE_LEVELS.size
    sst = np.full(shape, sst_value)
    psl = np.full(shape, 101000.0)
    temp = np.broadcast_to(
        TEMPERATURE_PROFILE.reshape((nlev,) + (1,) * len(shape)), (nlev,) + shape
    ).copy()
    mixr = np.broadcast_to(
        MIXING_RATIO_PROFILE.reshape((nlev,) + (1,) * len(shape)), (nlev,) + shape
    ).copy()
    return sst, psl, temp, mixr


def test_all_masked_cold_grid_reports_success_not_garbage_flag():
    """A grid with no computable point must not leak an unset status flag."""
    sst, psl, temp, mixr = _inputs(270.0)  # below the 5 degC TC threshold

    min_pressure, max_wind, error_flag = potential_intensity(
        sst, psl, PRESSURE_LEVELS, temp, mixr
    )

    assert error_flag == 1
    assert np.all(np.isnan(min_pressure))
    assert np.all(np.isnan(max_wind))


def test_all_masked_cold_4d_grid_reports_success():
    shape = (2, 3, 4)  # (time, y, x)
    nlev = PRESSURE_LEVELS.size
    sst = np.full(shape, 270.0)
    psl = np.full(shape, 101000.0)
    temp = np.broadcast_to(
        TEMPERATURE_PROFILE.reshape(1, -1, 1, 1), (shape[0], nlev) + shape[1:]
    ).copy()
    mixr = np.broadcast_to(
        MIXING_RATIO_PROFILE.reshape(1, -1, 1, 1), (shape[0], nlev) + shape[1:]
    ).copy()

    _, _, error_flag = potential_intensity(sst, psl, PRESSURE_LEVELS, temp, mixr)

    assert error_flag == 1


@pytest.mark.parametrize("repeat", range(5))
def test_warm_grid_error_flag_is_deterministic(repeat):
    sst, psl, temp, mixr = _inputs(302.0, shape=(16, 24))

    min_pressure, max_wind, error_flag = potential_intensity(
        sst, psl, PRESSURE_LEVELS, temp, mixr
    )

    assert error_flag == 1
    assert np.all(np.isfinite(min_pressure))
    assert np.all(max_wind > 0.0)
