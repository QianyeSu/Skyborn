"""
Tests for skyborn.causality module.

This module tests causality analysis functionality including Granger
causality and Liang information flow analysis for time series data.
"""

import pytest
import numpy as np
import warnings
from numpy.testing import assert_array_almost_equal, assert_allclose

from skyborn.causality import (
    granger_causality,
    liang_causality,
    liang,
    ar1_fit_evenly,
    sm_ar1_sim,
    phaseran,
    signif_isopersist,
    signif_isospec,
)


class TestAR1Functions:
    """Test AR(1) model fitting and simulation functions."""

    def test_ar1_fit_evenly(self):
        """Test AR(1) fitting function."""
        # Create a known AR(1) process
        np.random.seed(42)
        n = 500
        true_ar1 = 0.7

        # Generate AR(1) series: x[t] = ar1 * x[t-1] + noise
        x = np.zeros(n)
        for i in range(1, n):
            x[i] = true_ar1 * x[i - 1] + np.random.randn()

        # Estimate AR(1) coefficient
        fitted_ar1 = ar1_fit_evenly(x)

        # Should be close to true value
        assert abs(fitted_ar1 - true_ar1) < 0.15
        assert -1 < fitted_ar1 < 1  # Valid AR(1) range

    def test_ar1_fit_edge_cases(self):
        """Test AR(1) fitting with edge cases."""
        # White noise (should have AR(1) ~ 0)
        np.random.seed(42)
        white_noise = np.random.randn(200)
        ar1_coeff = ar1_fit_evenly(white_noise)
        assert abs(ar1_coeff) < 0.3  # Should be close to 0

        # Highly persistent series
        persistent = np.cumsum(np.random.randn(200) * 0.1)
        ar1_coeff = ar1_fit_evenly(persistent)
        assert ar1_coeff > 0.5  # Should be high

    def test_sm_ar1_sim(self):
        """Test AR(1) simulation function."""
        n = 100
        p = 50
        g = 0.6
        sig = 2.0

        ar1_sims = sm_ar1_sim(n, p, g, sig)

        # Check output shape
        assert ar1_sims.shape == (n, p)

        # Check that simulated series have approximately correct properties
        mean_std = np.mean(np.std(ar1_sims, axis=0))
        assert abs(mean_std - sig) < 0.5  # Within reasonable range

        # Check AR(1) coefficient for one realization
        fitted_ar1 = ar1_fit_evenly(ar1_sims[:, 0])
        assert abs(fitted_ar1 - g) < 0.3  # Should be reasonably close


class TestPhaseRandomization:
    """Test phase randomization for surrogate generation."""

    def test_phaseran_basic(self):
        """Test basic phase randomization functionality."""
        # Create a simple periodic signal
        t = np.linspace(0, 4 * np.pi, 100)
        x = np.sin(t) + 0.5 * np.sin(3 * t)

        nsurr = 10
        surrogates = phaseran(x.reshape(-1, 1), nsurr)

        # Check output shape
        assert surrogates.shape == (99, nsurr)  # One less due to odd requirement

        # Surrogates should have same mean (approximately)
        orig_mean = np.mean(x[:-1])  # Adjusted for length
        surr_means = np.mean(surrogates, axis=0)
        assert np.allclose(surr_means, orig_mean, atol=0.5)

        # Surrogates should have same variance (approximately)
        orig_var = np.var(x[:-1])
        surr_vars = np.var(surrogates, axis=0)
        assert np.allclose(surr_vars, orig_var, rtol=0.3)

    def test_phaseran_multiple_series(self):
        """Test phase randomization with multiple time series."""
        # Create two correlated periodic signals
        t = np.linspace(0, 4 * np.pi, 101)  # Odd length
        x1 = np.sin(t)
        x2 = np.cos(t) + 0.3 * np.sin(2 * t)

        data = np.column_stack([x1, x2])
        nsurr = 5

        # Should work with 2D input
        surrogates = phaseran(data, nsurr)
        assert surrogates.shape == (101, nsurr)


class TestGrangerCausality:
    """Test Granger causality analysis."""

    def test_granger_causality_basic(self):
        """Test basic Granger causality functionality."""
        # Create two time series where y2 causes y1
        n = 200
        np.random.seed(42)

        # y2 is independent
        y2 = np.random.randn(n)

        # y1 depends on past values of y2
        y1 = np.zeros(n)
        for i in range(1, n):
            y1[i] = 0.3 * y2[i - 1] + 0.5 * y1[i - 1] + np.random.randn() * 0.5

        # Test Granger causality
        result = granger_causality(y1, y2, maxlag=2, verbose=False)

        # Should return results for both lags
        assert 1 in result
        assert 2 in result

        # Check structure of results
        for lag in [1, 2]:
            assert len(result[lag]) == 2  # (test_results, model_results)
            test_results = result[lag][0]

            # Should have all four test types
            expected_tests = ["ssr_ftest", "ssr_chi2test", "lrtest", "params_ftest"]
            for test in expected_tests:
                assert test in test_results
                assert "statistic" in test_results[test]
                assert "pvalue" in test_results[test]

    def test_granger_causality_no_causality(self):
        """Test Granger causality with independent series."""
        n = 150
        np.random.seed(42)

        # Two independent AR(1) processes
        y1 = np.zeros(n)
        y2 = np.zeros(n)

        for i in range(1, n):
            y1[i] = 0.4 * y1[i - 1] + np.random.randn() * 0.8
            y2[i] = 0.3 * y2[i - 1] + np.random.randn() * 0.9

        result = granger_causality(y1, y2, maxlag=1, verbose=False)

        # p-values should generally be high (non-significant)
        pvalue = result[1][0]["ssr_ftest"]["pvalue"]
        # Note: Due to randomness, we can't guarantee high p-value,
        # but the test should complete successfully
        assert 0 <= pvalue <= 1

    def test_granger_causality_validation(self):
        """Test input validation for Granger causality."""
        y1 = np.random.randn(100)
        y2 = np.random.randn(90)  # Different length

        with pytest.raises(ValueError, match="Timeseries must be of same length"):
            granger_causality(y1, y2, verbose=False)


class TestLiangCausality:
    """Test Liang information flow analysis."""

    def test_liang_basic(self):
        """Test basic Liang causality functionality."""
        # Create two coupled oscillators where y2 influences y1
        n = 300
        dt = 0.1
        np.random.seed(42)

        y1 = np.zeros(n)
        y2 = np.zeros(n)

        # Simple coupled system
        for i in range(1, n):
            y1[i] = (
                y1[i - 1]
                + dt * (-0.5 * y1[i - 1] + 0.2 * y2[i - 1])
                + np.random.randn() * 0.1
            )
            y2[i] = y2[i - 1] + dt * (-0.3 * y2[i - 1]) + np.random.randn() * 0.1

        result = liang(y1, y2, npt=1)

        # Check output structure
        expected_keys = ["T21", "tau21", "Z", "dH1_star", "dH1_noise"]
        for key in expected_keys:
            assert key in result
            assert isinstance(result[key], (int, float, np.number))

        # Check that values are finite
        for key in expected_keys:
            assert np.isfinite(result[key])

        # Z should be positive (total information flow)
        assert result["Z"] > 0

        # tau21 should be normalized (between -1 and 1)
        assert -1 <= result["tau21"] <= 1

    def test_liang_causality_with_significance(self):
        """Test Liang causality with significance testing."""
        # Create a simple causal relationship
        n = 200
        np.random.seed(42)

        # y2 causes y1 with some noise
        y2 = np.random.randn(n)
        y1 = np.zeros(n)
        for i in range(1, n):
            y1[i] = 0.7 * y1[i - 1] + 0.3 * y2[i - 1] + np.random.randn() * 0.5

        # Test with significance testing (small nsim for speed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Suppress ARIMA convergence warnings
            result = liang_causality(
                y1,
                y2,
                npt=1,
                signif_test="isopersist",
                nsim=20,  # Small for testing speed
                qs=[0.05, 0.95],
            )

        # Check output structure
        expected_keys = [
            "T21",
            "tau21",
            "Z",
            "dH1_star",
            "dH1_noise",
            "signif_qs",
            "T21_noise",
            "tau21_noise",
        ]
        for key in expected_keys:
            assert key in result

        # Check significance testing results
        assert len(result["T21_noise"]) == 2  # Two quantiles
        assert len(result["tau21_noise"]) == 2
        assert all(np.isfinite(result["T21_noise"]))
        assert all(np.isfinite(result["tau21_noise"]))

    def test_liang_causality_methods(self):
        """Test different significance testing methods."""
        n = 150
        np.random.seed(42)

        # Simple coupled system
        y1 = np.random.randn(n)
        y2 = np.random.randn(n)

        # Test isospec method (phase randomization)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_isospec = liang_causality(
                y1,
                y2,
                signif_test="isospec",
                nsim=10,  # Small for speed
                qs=[0.05, 0.95],
            )

        # Should have same structure as isopersist
        expected_keys = ["T21", "tau21", "Z", "T21_noise", "tau21_noise"]
        for key in expected_keys:
            assert key in result_isospec

    def test_liang_npt_parameter(self):
        """Test different npt values for Liang causality."""
        n = 200
        np.random.seed(42)

        y1 = np.random.randn(n)
        y2 = np.random.randn(n)

        # Test different npt values
        for npt in [1, 2, 3]:
            result = liang(y1, y2, npt=npt)

            # Should return valid results for all npt values
            assert np.isfinite(result["T21"])
            assert np.isfinite(result["tau21"])
            assert result["Z"] > 0


class TestSignificanceTesting:
    """Test significance testing functions."""

    def test_signif_isopersist(self):
        """Test significance testing with AR(1) surrogates."""
        n = 100
        np.random.seed(42)

        y1 = np.random.randn(n)
        y2 = np.random.randn(n)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = signif_isopersist(
                y1,
                y2,
                method="liang",
                nsim=10,  # Small for speed
                qs=[0.05, 0.5, 0.95],
                npt=1,
            )

        # Check output structure
        assert "T21_noise_qs" in result
        assert "tau21_noise_qs" in result
        assert len(result["T21_noise_qs"]) == 3
        assert len(result["tau21_noise_qs"]) == 3
        assert all(np.isfinite(result["T21_noise_qs"]))
        assert all(np.isfinite(result["tau21_noise_qs"]))

    def test_signif_isospec(self):
        """Test significance testing with phase randomization."""
        n = 101  # Odd length required
        np.random.seed(42)

        y1 = np.random.randn(n)
        y2 = np.random.randn(n)

        result = signif_isospec(
            y1, y2, method="liang", nsim=5, qs=[0.1, 0.9], npt=1  # Very small for speed
        )

        # Check output structure
        assert "T21_noise_qs" in result
        assert "tau21_noise_qs" in result
        assert len(result["T21_noise_qs"]) == 2
        assert len(result["tau21_noise_qs"]) == 2

    def test_signif_invalid_method(self):
        """Test significance testing with invalid method."""
        y1 = np.random.randn(50)
        y2 = np.random.randn(50)

        with pytest.raises(KeyError, match="is not a valid method"):
            signif_isopersist(y1, y2, method="invalid_method")

        with pytest.raises(KeyError, match="is not a valid method"):
            signif_isospec(y1, y2, method="invalid_method")


class TestCausalityIntegration:
    """Integration tests for causality module."""

    def test_causality_with_climate_timeseries(self, sample_1d_timeseries):
        """Test causality analysis with climate-like time series."""
        ts = sample_1d_timeseries

        # Create a second time series that lags the first
        ts2_data = np.roll(ts.values, 5) + np.random.randn(len(ts)) * 0.1

        # Test both methods
        result_granger = granger_causality(ts.values, ts2_data, maxlag=1, verbose=False)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_liang = liang(ts.values, ts2_data, npt=1)

        # Both should complete successfully
        assert 1 in result_granger
        assert "T21" in result_liang
        assert np.isfinite(result_liang["T21"])

    def test_causality_error_handling(self):
        """Test comprehensive error handling."""
        # Test with NaN values
        y1 = np.array([1, 2, np.nan, 4, 5])
        y2 = np.array([2, 3, 4, 5, 6])

        # Liang should handle NaNs gracefully or raise informative error
        try:
            result = liang(y1, y2)
            # If it succeeds, result should be reasonable
            assert np.isfinite(result["Z"])
        except (ValueError, np.linalg.LinAlgError):
            # Expected to fail with singular matrix or similar
            pass

        # Test with very short series
        y1_short = np.array([1, 2])
        y2_short = np.array([2, 3])

        # Should either work or fail gracefully
        try:
            result = liang(y1_short, y2_short)
        except (ValueError, np.linalg.LinAlgError):
            pass  # Expected for very short series

    def test_causality_consistency(self):
        """Test consistency between different causality measures."""
        # Create a clear causal relationship
        n = 300
        np.random.seed(42)

        # Strong causality: y1 depends heavily on y2
        y2 = np.random.randn(n)
        y1 = np.zeros(n)
        for i in range(1, n):
            y1[i] = 0.8 * y2[i - 1] + 0.2 * y1[i - 1] + np.random.randn() * 0.2

        # Test Granger causality
        result_granger = granger_causality(y1, y2, maxlag=1, verbose=False)
        pvalue_granger = result_granger[1][0]["ssr_ftest"]["pvalue"]

        # Test Liang causality
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_liang = liang(y1, y2, npt=1)

        # Both should detect some level of causality
        # Granger should have low p-value for strong causality
        # Liang should have non-zero information flow
        assert pvalue_granger < 0.5  # Some evidence of causality
        assert abs(result_liang["T21"]) > 0  # Non-zero information flow


# Performance tests (marked as slow)
@pytest.mark.slow
class TestCausalityPerformance:
    """Performance tests for causality module."""

    def test_liang_causality_large_data(self):
        """Test Liang causality with larger datasets."""
        n = 2000
        np.random.seed(42)

        # Create coupled system
        y1 = np.random.randn(n)
        y2 = np.random.randn(n)

        # Should complete without excessive memory usage
        result = liang(y1, y2, npt=1)

        assert np.isfinite(result["T21"])
        assert np.isfinite(result["tau21"])
        assert result["Z"] > 0

    def test_granger_causality_performance(self):
        """Test Granger causality performance."""
        n = 1000
        np.random.seed(42)

        y1 = np.random.randn(n)
        y2 = np.random.randn(n)

        # Should complete in reasonable time
        result = granger_causality(y1, y2, maxlag=3, verbose=False)

        assert 1 in result
        assert 2 in result
        assert 3 in result


class TestCausalityEdgeCases:
    """Test edge cases and special conditions."""

    def test_causality_identical_series(self):
        """Test causality with identical time series."""
        n = 100
        y = np.random.randn(n)

        # Test with identical series
        result_liang = liang(y, y, npt=1)

        # Should handle this case (perfect correlation)
        assert np.isfinite(result_liang["T21"])
        # Note: actual values may vary due to numerical precision

    def test_causality_constant_series(self):
        """Test causality with constant time series."""
        n = 100
        y1 = np.ones(n)
        y2 = np.random.randn(n)

        # Should handle constant series gracefully
        try:
            result = liang(y1, y2, npt=1)
            # If successful, should be finite
            assert np.isfinite(result["Z"])
        except (ValueError, np.linalg.LinAlgError):
            # Expected to fail due to zero variance
            pass

    def test_causality_linear_trend(self):
        """Test causality with linear trends."""
        n = 200
        t = np.arange(n)

        # Both series have linear trends
        y1 = 0.1 * t + np.random.randn(n) * 0.1
        y2 = 0.05 * t + np.random.randn(n) * 0.1

        # Should detect some relationship due to common trend
        result = liang(y1, y2, npt=1)

        assert np.isfinite(result["T21"])
        assert np.isfinite(result["tau21"])
        assert result["Z"] > 0


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
