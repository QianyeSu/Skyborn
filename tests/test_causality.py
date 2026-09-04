"""
Tests for skyborn.calc.causality module.

This module tests causality analysis functionality including Granger
causality and Liang information flow analysis for time series data.
"""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_almost_equal
from statsmodels.tsa.arima_process import arma_generate_sample

import skyborn.calc.causality.core as causality_core
from skyborn.calc.causality import (
    Surrogate,
    ar1_fit_evenly,
    granger_causality,
    liang,
    liang_causality,
    phaseran,
    signif_isopersist,
    signif_isospec,
    sm_ar1_sim,
)
from skyborn.calc.causality.core import (
    _ar1_native_available,
    _liang_backend_available,
    _liang_batch,
    _liang_batch_backend,
    _liang_python,
    _load_liang_backend,
    _require_liang_backend,
)


class TestBackendLoading:
    """Test compiled backend discovery and failure handling."""

    @staticmethod
    def _fake_path(exists):
        class FakePath:
            def resolve(self):
                return self

            @property
            def parent(self):
                return self

            def __truediv__(self, _other):
                return self

            def exists(self):
                return exists

        return FakePath()

    def test_load_backend_uses_import_fallback(self, monkeypatch):
        sentinel = object()
        monkeypatch.setattr(
            causality_core,
            "Path",
            lambda _path: self._fake_path(exists=False),
        )
        monkeypatch.setattr(causality_core, "import_module", lambda _name: sentinel)

        assert _load_liang_backend() is sentinel

    def test_load_backend_returns_none_when_import_fails(self, monkeypatch):
        monkeypatch.setattr(
            causality_core,
            "Path",
            lambda _path: self._fake_path(exists=False),
        )

        def fail_import(_name):
            raise ImportError("backend unavailable")

        monkeypatch.setattr(causality_core, "import_module", fail_import)
        assert _load_liang_backend() is None

    def test_load_backend_skips_invalid_spec(self, monkeypatch):
        monkeypatch.setattr(
            causality_core,
            "Path",
            lambda _path: self._fake_path(exists=True),
        )
        monkeypatch.setattr(
            causality_core.importlib_util,
            "spec_from_file_location",
            lambda *_args: None,
        )
        monkeypatch.setattr(causality_core, "import_module", lambda _name: None)

        assert _load_liang_backend() is None

    def test_load_backend_skips_loader_exception(self, monkeypatch):
        class FailingLoader:
            def exec_module(self, _module):
                raise RuntimeError("load failed")

        class FakeSpec:
            loader = FailingLoader()

        monkeypatch.setattr(
            causality_core,
            "Path",
            lambda _path: self._fake_path(exists=True),
        )
        monkeypatch.setattr(
            causality_core.importlib_util,
            "spec_from_file_location",
            lambda *_args: FakeSpec(),
        )

        def fail_import(_name):
            raise ImportError("backend unavailable")

        monkeypatch.setattr(causality_core, "import_module", fail_import)
        assert _load_liang_backend() is None

    def test_require_backend_raises_when_missing(self, monkeypatch):
        monkeypatch.setattr(causality_core, "_LIANG_BACKEND", None)
        with pytest.raises(ImportError, match="backend is unavailable"):
            _require_liang_backend()


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

    def test_ar1_fit_caps_coefficient_above_one(self, monkeypatch, capsys):
        """AR(1) fitting should cap an out-of-range fitted coefficient."""

        class FakeFit:
            params = np.array([0.0, 0.0, 1.5])

        class FakeARIMA:
            def __init__(self, *args, **kwargs):
                pass

            def fit(self):
                return FakeFit()

        monkeypatch.setattr(causality_core, "ARIMA", FakeARIMA)
        result = ar1_fit_evenly(np.arange(10.0))

        assert result == 1.0 - np.spacing(1.0) ** (1 / 4)
        assert "greater than 1" in capsys.readouterr().out

    def test_sm_ar1_sim(self):
        """Test AR(1) simulation function."""
        n = 100
        p = 50
        g = 0.6
        sig = 2.0

        ar1_sims = sm_ar1_sim(n, p, g, sig)

        # Check output shape
        assert ar1_sims.shape == (n, p)
        assert ar1_sims.flags.f_contiguous

        # Check that simulated series have approximately correct properties
        mean_std = np.mean(np.std(ar1_sims, axis=0))
        assert abs(mean_std - sig) < 0.5  # Within reasonable range

        # Check AR(1) coefficient for one realization
        fitted_ar1 = ar1_fit_evenly(ar1_sims[:, 0])
        assert abs(fitted_ar1 - g) < 0.3  # Should be reasonably close

    def test_sm_ar1_sim_matches_legacy_random_stream(self):
        """Batch AR(1) generation should preserve the historical output."""
        n = 64
        p = 5
        g = 0.6
        sig = 2.0
        ar = np.r_[1, -g]
        ma = np.r_[1, 0.0]
        sig_n = sig * np.sqrt(1 - g**2)

        np.random.seed(123)
        actual = sm_ar1_sim(n, p, g, sig)

        np.random.seed(123)
        expected = np.empty((n, p))
        for index in range(p):
            expected[:, index] = arma_generate_sample(
                ar=ar,
                ma=ma,
                nsample=n,
                burnin=50,
                scale=sig_n,
            )

        assert_allclose(actual, expected, rtol=0.0, atol=0.0)

    def test_sm_ar1_sim_zero_surrogates(self):
        """Zero requested surrogates should return an empty matrix."""
        result = sm_ar1_sim(n=32, p=0, g=0.6, sig=2.0)
        assert result.shape == (32, 0)
        assert result.flags.f_contiguous

    def test_sm_ar1_sim_rejects_invalid_backend(self):
        with pytest.raises(ValueError, match="backend must be one of"):
            sm_ar1_sim(32, 4, 0.6, 2.0, backend="invalid")

    def test_sm_ar1_sim_rejects_invalid_dimensions(self):
        with pytest.raises(ValueError, match="n must be positive"):
            sm_ar1_sim(0, 4, 0.6, 2.0)
        with pytest.raises(ValueError, match="n must be positive"):
            sm_ar1_sim(32, -1, 0.6, 2.0)

    @pytest.mark.skipif(
        not _ar1_native_available(),
        reason="compiled AR(1) backend is not available",
    )
    def test_sm_ar1_sim_fortran_matches_legacy_random_stream(self):
        n = 64
        p = 48
        g = 0.6
        sig = 2.0
        ar = np.r_[1, -g]
        ma = np.r_[1, 0.0]
        sig_n = sig * np.sqrt(1 - g**2)

        np.random.seed(123)
        actual = sm_ar1_sim(n, p, g, sig, backend="fortran")

        np.random.seed(123)
        expected = np.empty((n, p))
        for index in range(p):
            expected[:, index] = arma_generate_sample(
                ar=ar,
                ma=ma,
                nsample=n,
                burnin=50,
                scale=sig_n,
            )

        assert_allclose(actual, expected, rtol=0.0, atol=0.0)

    def test_sm_ar1_sim_fortran_requires_compiled_backend(self, monkeypatch):
        monkeypatch.setattr(causality_core, "_LIANG_BACKEND", None)
        with pytest.raises(ImportError, match="compiled AR"):
            sm_ar1_sim(32, 4, 0.6, 2.0, backend="fortran")

    def test_sm_ar1_sim_auto_falls_back_without_backend(self, monkeypatch):
        monkeypatch.setattr(causality_core, "_LIANG_BACKEND", None)
        np.random.seed(123)
        result = sm_ar1_sim(32, 32, 0.6, 2.0, backend="auto")
        assert result.shape == (32, 32)

    def test_ar1_filter_native_requires_entry_point(self, monkeypatch):
        monkeypatch.setattr(causality_core, "_LIANG_BACKEND", object())
        with pytest.raises(ImportError, match="AR\\(1\\) filtering backend"):
            causality_core._ar1_filter_native(
                np.ones((51, 1), dtype=np.float64),
                0.2,
                1,
            )


class TestPhaseRandomization:
    """Test phase randomization for surrogate generation."""

    @staticmethod
    def _legacy_phaseran(recblk, nsurr):
        recblk = np.asarray(recblk, dtype=np.float64)
        was_1d = recblk.ndim == 1
        if was_1d:
            recblk = recblk[:, np.newaxis]

        nfrms = recblk.shape[0]
        if nfrms % 2 == 0:
            nfrms -= 1
            recblk = recblk[:nfrms]

        len_ser = int((nfrms - 1) / 2)
        fft_recblk = np.fft.fft(recblk, axis=0)
        surrblk = np.zeros((nfrms, recblk.shape[1], nsurr))
        for k in range(nsurr):
            ph_rnd = np.random.rand(len_ser, 1)
            ph_interv1 = np.exp(2 * np.pi * 1j * ph_rnd)
            ph_interv2 = np.conj(np.flipud(ph_interv1))
            fft_recblk_surr = np.copy(fft_recblk)
            fft_recblk_surr[1 : len_ser + 1] = fft_recblk[1 : len_ser + 1] * ph_interv1
            fft_recblk_surr[len_ser + 1 :] = fft_recblk[len_ser + 1 :] * ph_interv2
            surrblk[:, :, k] = np.real(np.fft.ifft(fft_recblk_surr, axis=0))

        if was_1d:
            return surrblk[:, 0, :]
        return surrblk

    def test_phaseran_matches_legacy_random_stream(self):
        """Batch FFT generation should preserve the historical phase stream."""
        rng = np.random.default_rng(7)
        data = rng.normal(size=(101, 2))

        np.random.seed(123)
        expected = self._legacy_phaseran(data, 17)
        np.random.seed(123)
        actual = phaseran(data, 17)

        assert_allclose(actual, expected, rtol=0.0, atol=1e-14)

    def test_phaseran_multiple_batches(self):
        """Surrogate counts above the block size should use all blocks."""
        rng = np.random.default_rng(8)
        data = rng.normal(size=(33, 2))
        nsurr = 257

        np.random.seed(321)
        expected = self._legacy_phaseran(data, nsurr)
        np.random.seed(321)
        actual = phaseran(data, nsurr)

        assert actual.shape == (33, 2, nsurr)
        assert_allclose(actual, expected, rtol=0.0, atol=1e-14)

    def test_phaseran_basic(self):
        """Test basic phase randomization functionality."""
        # Create a simple periodic signal
        t = np.linspace(0, 4 * np.pi, 100)
        x = np.sin(t) + 0.5 * np.sin(3 * t)

        nsurr = 10
        surrogates = phaseran(x.reshape(-1, 1), nsurr)

        # Check output shape
        assert surrogates.shape == (99, 1, nsurr)  # One less due to odd requirement

        # Surrogates should have same mean (approximately)
        orig_mean = np.mean(x[:-1])  # Adjusted for length
        surr_means = np.mean(surrogates[:, 0, :], axis=0)
        assert np.allclose(surr_means, orig_mean, atol=0.5)

        # Surrogates should have same variance (approximately)
        orig_var = np.var(x[:-1])
        surr_vars = np.var(surrogates[:, 0, :], axis=0)
        assert np.allclose(surr_vars, orig_var, rtol=0.3)

    def test_phaseran_1d_output_is_fortran_contiguous(self):
        """The 1D batch layout should be ready for the native Liang kernel."""
        surrogates = phaseran(np.arange(33.0), nsurr=4)

        assert surrogates.shape == (33, 4)
        assert surrogates.flags.f_contiguous

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
        assert surrogates.shape == (101, 2, nsurr)

    def test_phaseran_input_validation(self):
        """Phase randomization should reject unsupported inputs."""
        with pytest.raises(ValueError, match="1D or 2D"):
            phaseran(np.zeros((2, 2, 2)), 3)

        with pytest.raises(ValueError, match="positive integer"):
            phaseran(np.zeros(9), 0)


class TestSurrogate:
    """Test the unified surrogate-generation facade."""

    def test_phase_randomized_alias_preserves_output(self):
        data = np.random.default_rng(41).normal(size=33)
        np.random.seed(123)
        expected = phaseran(data, 4)
        np.random.seed(123)
        actual = Surrogate(data).generate("phase_randomized", 4)
        assert_allclose(actual, expected, rtol=0.0, atol=1e-14)

    def test_ar1_generation_and_alias(self):
        data = np.random.default_rng(42).normal(size=64)
        surrogate = Surrogate(data)
        result = surrogate.generate(
            "ar1",
            4,
            backend="python",
            g=0.4,
            sig=1.5,
        )
        assert result.shape == (64, 4)
        assert result.flags.f_contiguous

    def test_ar1_generation_fits_defaults(self):
        data = np.random.default_rng(43).normal(size=64)
        result = Surrogate(data).ar1(2, backend="python")
        assert result.shape == (64, 2)

    def test_surrogate_validation(self):
        with pytest.raises(ValueError, match="1D or 2D"):
            Surrogate(np.zeros((2, 2, 2)))
        with pytest.raises(ValueError, match="at least one"):
            Surrogate(np.empty((0,)))
        with pytest.raises(ValueError, match="finite"):
            Surrogate(np.array([0.0, np.nan]))
        with pytest.raises(ValueError, match="one-dimensional"):
            Surrogate(np.zeros((8, 2))).ar1(2, g=0.2, sig=1.0)

    def test_surrogate_method_and_keyword_validation(self):
        surrogate = Surrogate(np.arange(9.0))
        with pytest.raises(KeyError, match="valid surrogate"):
            surrogate.generate("unknown", 2)
        with pytest.raises(TypeError, match="Unexpected keyword"):
            surrogate.generate("isospec", 2, g=0.2)


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
                assert len(test_results[test]) >= 2
                assert np.isfinite(test_results[test][0])
                assert 0 <= test_results[test][1] <= 1

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
        pvalue = result[1][0]["ssr_ftest"][1]
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

    @pytest.mark.skipif(
        not _liang_backend_available(),
        reason="compiled Liang backend is not available",
    )
    def test_compiled_liang_matches_python_reference(self):
        """Compiled single and batch kernels should match the reference path."""
        rng = np.random.default_rng(1234)
        y1 = rng.normal(size=256)
        y2 = 0.4 * y1 + rng.normal(size=256)

        compiled = liang(y1, y2, npt=2)
        reference = _liang_python(y1, y2, npt=2)
        for key in ("T21", "tau21", "Z", "dH1_star", "dH1_noise"):
            assert_allclose(compiled[key], reference[key], rtol=1e-10, atol=1e-12)

        batch_y1 = np.column_stack([y1, rng.normal(size=256), rng.normal(size=256)])
        batch_y2 = np.column_stack([y2, rng.normal(size=256), rng.normal(size=256)])
        t21, tau21 = _liang_batch_backend(batch_y1, batch_y2, npt=2)
        expected_t21 = np.array(
            [
                _liang_python(batch_y1[:, i], batch_y2[:, i], npt=2)["T21"]
                for i in range(3)
            ]
        )
        expected_tau21 = np.array(
            [
                _liang_python(batch_y1[:, i], batch_y2[:, i], npt=2)["tau21"]
                for i in range(3)
            ]
        )
        assert_allclose(t21, expected_t21, rtol=1e-10, atol=1e-12)
        assert_allclose(tau21, expected_tau21, rtol=1e-10, atol=1e-12)

    def test_python_liang_input_validation(self):
        """The Python Liang reference should validate all input constraints."""
        with pytest.raises(ValueError, match="1D arrays"):
            _liang_python(np.zeros((4, 1)), np.zeros(4))
        with pytest.raises(ValueError, match="share a length"):
            _liang_python(np.zeros(4), np.zeros(5))
        with pytest.raises(ValueError, match="at least two"):
            _liang_python(np.zeros(2), np.zeros(2), npt=1)
        with pytest.raises(ValueError, match="finite"):
            _liang_python(np.array([1.0, np.nan, 2.0]), np.arange(3.0))

    def test_python_liang_singular_covariance(self):
        """The Python Liang reference should reject singular covariance."""
        y = np.arange(10.0)
        with pytest.raises(ValueError, match="singular covariance"):
            _liang_python(y, y)

    def test_liang_input_validation(self):
        """The public Liang function should validate all input constraints."""
        with pytest.raises(ValueError, match="1D arrays"):
            liang(np.zeros((4, 1)), np.zeros(4))
        with pytest.raises(ValueError, match="share a length"):
            liang(np.zeros(4), np.zeros(5))
        with pytest.raises(ValueError, match="at least two"):
            liang(np.zeros(2), np.zeros(2), npt=1)
        with pytest.raises(ValueError, match="finite"):
            liang(np.array([1.0, np.nan, 2.0]), np.arange(3.0))

    def test_liang_python_fallback_and_batch_fallback(self, monkeypatch):
        """The uncompiled Liang paths should remain functional and validated."""
        rng = np.random.default_rng(11)
        y1 = rng.normal(size=32)
        y2 = rng.normal(size=32)
        batch_y1 = np.column_stack([y1, rng.normal(size=32)])
        batch_y2 = np.column_stack([y2, rng.normal(size=32)])

        monkeypatch.setattr(causality_core, "_LIANG_BACKEND", None)
        result = liang(y1, y2)
        reference = _liang_python(y1, y2)
        assert_allclose(result["T21"], reference["T21"])

        t21, tau21 = _liang_batch(batch_y1, batch_y2, npt=1)
        assert t21.shape == (2,)
        assert tau21.shape == (2,)

        with pytest.raises(ValueError, match="2D arrays"):
            _liang_batch(batch_y1[:, :1], batch_y2, npt=1)


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

    def test_liang_causality_invalid_significance_method(self):
        """Liang causality should reject unknown significance methods."""
        y1 = np.random.default_rng(12).normal(size=32)
        y2 = np.random.default_rng(13).normal(size=32)
        with pytest.raises(KeyError, match="not a valid significance"):
            liang_causality(y1, y2, signif_test="invalid")


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
        pvalue_granger = result_granger[1][0]["ssr_ftest"][1]

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

        # Perfectly correlated inputs make the Liang covariance inversion
        # undefined and should fail with a clear validation error.
        with pytest.raises(ValueError, match="singular covariance"):
            liang(y, y, npt=1)

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
