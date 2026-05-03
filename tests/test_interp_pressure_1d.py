"""Focused tests for one-dimensional pressure-coordinate interpolation."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose, assert_array_equal

import skyborn.interp.interpolation as interpolation_mod
from skyborn.interp.interpolation import interp_pressure_1d


def _log_interp(x0: float, x1: float, p0: float, p1: float, p: float) -> float:
    """Return the expected legacy log-pressure interpolation result."""

    slope = (x0 - x1) / (np.log(p0) - np.log(p1))
    return x1 + slope * (np.log(p) - np.log(p1))


def _compat_low_end_log_extrap(
    x_low: float, x_high: float, p_low: float, p_high: float, p: float
) -> float:
    """Return the archived low-end log extrapolation result from ``int2p_dp.f``."""

    slope = (x_low - x_high) / (np.log(p_low) - np.log(p_high))
    return x_low + slope * (np.log(p) - np.log(p_high))


class TestInterpPressure1D:
    """Behavior checks for `interp_pressure_1d`."""

    def test_descriptive_keyword_names_are_supported(self):
        """The public API should accept descriptive pressure-interpolation keywords."""

        result = interp_pressure_1d(
            values=np.array([10.0, 20.0, 40.0]),
            source_pressure=np.array([1000.0, 850.0, 700.0]),
            target_pressure=np.array([925.0, 775.0]),
            method="linear",
        )

        assert_allclose(result, np.array([15.0, 30.0]), rtol=0.0, atol=1e-12)

    def test_legacy_keyword_aliases_remain_supported(self):
        """The legacy keyword aliases should remain valid for compatibility."""

        result = interp_pressure_1d(
            x=np.array([10.0, 20.0, 40.0]),
            p_in=np.array([1000.0, 850.0, 700.0]),
            p_out=np.array([925.0, 775.0]),
            method="linear",
        )

        assert_allclose(result, np.array([15.0, 30.0]), rtol=0.0, atol=1e-12)

    def test_linear_interpolation_matches_hand_calculation(self):
        """Linear-in-pressure mode should match a small hand-computed example."""

        x = np.array([10.0, 20.0, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([925.0, 775.0])

        result = interp_pressure_1d(x, p_in, p_out, method="linear")

        assert_allclose(result, np.array([15.0, 30.0]), rtol=0.0, atol=1e-12)

    def test_log_interpolation_matches_hand_calculation(self):
        """Log-pressure mode should match the legacy log interpolation formula."""

        x = np.array([10.0, 20.0, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([925.0, 775.0])

        result = interp_pressure_1d(x, p_in, p_out, method="log")
        expected = np.array(
            [
                _log_interp(10.0, 20.0, 1000.0, 850.0, 925.0),
                _log_interp(20.0, 40.0, 850.0, 700.0, 775.0),
            ]
        )

        assert_allclose(result, expected, rtol=0.0, atol=1e-12)

    def test_input_and_output_pressure_order_do_not_change_results(self):
        """Increasing and decreasing pressure axes should behave equivalently."""

        x = np.array([40.0, 20.0, 10.0])
        p_in = np.array([700.0, 850.0, 1000.0])
        p_out = np.array([775.0, 925.0])

        result = interp_pressure_1d(x, p_in, p_out, method="linear")

        assert_allclose(result, np.array([30.0, 15.0]), rtol=0.0, atol=1e-12)

    def test_out_of_range_targets_stay_missing_without_extrapolation(self):
        """Targets outside the source range should remain missing by default."""

        x = np.array([10.0, 20.0, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([1050.0, 650.0])

        result = interp_pressure_1d(x, p_in, p_out, method="linear")

        assert np.isnan(result[0])
        assert np.isnan(result[1])

    def test_extrapolation_uses_the_end_slopes(self):
        """Extrapolation should use the first and last valid bracketing slopes."""

        x = np.array([10.0, 20.0, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([1050.0, 650.0])

        result = interp_pressure_1d(
            x,
            p_in,
            p_out,
            method="linear",
            extrapolate=True,
        )

        expected = np.array([10.0 - 50.0 / 15.0, 40.0 + 100.0 / 3.0 / 5.0])
        assert_allclose(result, expected, rtol=0.0, atol=1e-12)

    def test_low_end_log_extrapolation_matches_archived_kernel(self):
        """Low-end log extrapolation should preserve the archived F77 behavior."""

        x = np.array([35.0, -9999.0, 20.0, 14.0, 10.0])
        p_in = np.array([700.0, 775.0, 850.0, 925.0, 1000.0])
        p_out = np.array([650.0, 700.0, 900.0, 1025.0])

        result = interp_pressure_1d(
            x,
            p_in,
            p_out,
            method="log",
            extrapolate=True,
            missing_value=-9999.0,
        )

        expected = np.array(
            [
                _compat_low_end_log_extrap(35.0, 20.0, 700.0, 850.0, 650.0),
                35.0,
                _log_interp(14.0, 20.0, 925.0, 850.0, 900.0),
                _log_interp(10.0, 14.0, 1000.0, 925.0, 1025.0),
            ]
        )

        assert_allclose(result, expected, rtol=0.0, atol=1e-12)

    def test_missing_input_levels_are_ignored_before_interpolation(self):
        """NaN or sentinel-marked source levels should be dropped cleanly."""

        x = np.array([10.0, np.nan, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([850.0])

        result = interp_pressure_1d(x, p_in, p_out, method="linear")

        assert_allclose(result, np.array([25.0]), rtol=0.0, atol=1e-12)

    def test_exact_pressure_matches_round_trip_exactly(self):
        """Exact source-pressure matches should copy the source value exactly."""

        x = np.array([10.0, 20.0, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([700.0, 850.0, 1000.0])

        result = interp_pressure_1d(x, p_in, p_out, method="log")

        assert_array_equal(result, np.array([40.0, 20.0, 10.0]))

    def test_custom_missing_value_is_preserved_in_output(self):
        """Finite custom missing values should round-trip through the wrapper."""

        x = np.array([10.0, -9999.0, 40.0])
        p_in = np.array([1000.0, 850.0, 700.0])
        p_out = np.array([850.0, 1100.0])

        result = interp_pressure_1d(
            x,
            p_in,
            p_out,
            method="linear",
            extrapolate=False,
            missing_value=-9999.0,
        )

        assert_allclose(result[:1], np.array([25.0]), rtol=0.0, atol=1e-12)
        assert result[1] == -9999.0

    def test_xarray_input_preserves_name_attrs_and_target_coordinate(self):
        """A 1D DataArray input should return a 1D DataArray with pressure coords."""

        x = xr.DataArray(
            [10.0, 20.0, 40.0],
            dims=["plev"],
            coords={"plev": [1000.0, 850.0, 700.0]},
            name="u",
            attrs={"units": "m s-1"},
        )
        p_out = xr.DataArray(
            [925.0, 775.0],
            dims=["target_plev"],
            coords={"target_plev": [925.0, 775.0]},
            attrs={"units": "hPa"},
        )

        result = interp_pressure_1d(x, x["plev"], p_out, method="linear")

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("target_plev",)
        assert_array_equal(result["target_plev"].values, np.array([925.0, 775.0]))
        assert result.name == "u"
        assert result.attrs["units"] == "m s-1"
        assert_allclose(result.values, np.array([15.0, 30.0]), rtol=0.0, atol=1e-12)

    def test_wrapper_fails_clearly_without_compiled_backend(self, monkeypatch):
        """The public wrapper should explain when the compiled backend is unavailable."""

        monkeypatch.setattr(interpolation_mod, "_dinterp_pressure_1d", None)

        with pytest.raises(RuntimeError, match="compiled Fortran backend"):
            interp_pressure_1d(
                np.array([10.0, 20.0]),
                np.array([1000.0, 900.0]),
                np.array([950.0]),
            )
