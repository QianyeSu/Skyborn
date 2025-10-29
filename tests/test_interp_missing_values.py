"""
Comprehensive tests for skyborn.interp.missing_values module.

Tests py2fort_msg and fort2py_msg functions with all supported data types.
Target: 95%+ coverage for missing_values.py
"""

import numpy as np
import pytest

from skyborn.interp.missing_values import (
    complex_dtypes,
    float_dtypes,
    fort2py_msg,
    int_dtypes,
    msg_dtype,
    py2fort_msg,
    string_dtypes,
    supported_dtypes,
    uint_dtypes,
)


class TestMissingValueConstants:
    """Test module-level constants and dictionaries."""

    def test_msg_dtype_keys(self):
        """Test msg_dtype contains expected types."""
        expected_types = [
            np.complex64,
            np.complex128,
            np.float16,
            np.float32,
            np.float64,
            np.int8,
            np.int16,
            np.int32,
            np.int64,
            np.uint8,
            np.uint16,
            np.uint32,
            np.uint64,
            str,
        ]
        for dtype in expected_types:
            assert dtype in msg_dtype

    def test_msg_dtype_values_are_correct_type(self):
        """Test msg_dtype values match their key types."""
        for dtype, value in msg_dtype.items():
            if dtype == str:
                assert isinstance(value, str)
            else:
                assert isinstance(value, dtype)

    def test_complex_dtypes_list(self):
        """Test complex_dtypes list."""
        assert np.complex64 in complex_dtypes
        assert np.complex128 in complex_dtypes
        assert len(complex_dtypes) == 2

    def test_float_dtypes_list(self):
        """Test float_dtypes list."""
        assert np.float16 in float_dtypes
        assert np.float32 in float_dtypes
        assert np.float64 in float_dtypes
        assert len(float_dtypes) == 3

    def test_int_dtypes_list(self):
        """Test int_dtypes list."""
        assert np.int8 in int_dtypes
        assert np.int16 in int_dtypes
        assert np.int32 in int_dtypes
        assert np.int64 in int_dtypes
        assert len(int_dtypes) == 4

    def test_uint_dtypes_list(self):
        """Test uint_dtypes list."""
        assert np.uint8 in uint_dtypes
        assert np.uint16 in uint_dtypes
        assert np.uint32 in uint_dtypes
        assert np.uint64 in uint_dtypes
        assert len(uint_dtypes) == 4

    def test_string_dtypes_list(self):
        """Test string_dtypes list."""
        assert str in string_dtypes
        assert len(string_dtypes) == 1

    def test_supported_dtypes(self):
        """Test supported_dtypes matches msg_dtype keys."""
        assert set(supported_dtypes) == set(msg_dtype.keys())


class TestPy2FortMsg:
    """Test py2fort_msg function."""

    def test_float32_default_missing_value(self):
        """Test float32 array with default NaN handling."""
        data = np.array([1.0, 2.0, np.nan, 4.0], dtype=np.float32)
        result, msg_py, msg_fort = py2fort_msg(data)

        assert np.isnan(msg_py)
        assert msg_fort == msg_dtype[np.float32]
        assert result[2] == msg_fort

    def test_float64_explicit_missing_value(self):
        """Test float64 with explicit missing value."""
        data = np.array([1.0, -999.0, 3.0], dtype=np.float64)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=-999.0)

        assert msg_py == -999.0
        assert msg_fort == msg_dtype[np.float64]
        assert result[1] == msg_fort

    def test_int32_default_missing_value(self):
        """Test int32 with default missing value."""
        data = np.array([1, 2, 127, 4], dtype=np.int32)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=127)

        assert msg_py == 127
        assert msg_fort == msg_dtype[np.int32]
        assert result[2] == msg_fort

    def test_int8_with_explicit_fortran_value(self):
        """Test int8 with explicit Fortran missing value."""
        data = np.array([1, -128, 3], dtype=np.int8)
        custom_fort = np.int8(99)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=-128, msg_fort=custom_fort)

        assert msg_py == -128
        assert msg_fort == custom_fort
        assert result[1] == 99

    def test_uint16_array(self):
        """Test uint16 array."""
        data = np.array([100, 200, 65535, 400], dtype=np.uint16)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=65535)

        assert msg_py == 65535
        assert result[2] == msg_dtype[np.uint16]

    def test_complex64_with_nan(self):
        """Test complex64 with NaN."""
        data = np.array([1 + 2j, np.nan + np.nan * 1j, 3 + 4j], dtype=np.complex64)
        result, msg_py, msg_fort = py2fort_msg(data)

        assert np.isnan(msg_py)
        assert msg_fort == msg_dtype[np.complex64]

    def test_complex128_custom_missing(self):
        """Test complex128 with custom missing value."""
        data = np.array([1 + 1j, -999 - 999j, 2 + 2j], dtype=np.complex128)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=-999 - 999j)

        assert msg_py == -999 - 999j
        assert result[1] == msg_dtype[np.complex128]

    def test_no_missing_values(self):
        """Test array with no missing values."""
        data = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=-999.0)

        # Should not modify data if no missing values found
        assert not np.any(result == msg_fort)

    def test_all_missing_values(self):
        """Test array with all missing values."""
        data = np.array([np.nan, np.nan, np.nan], dtype=np.float64)
        result, msg_py, msg_fort = py2fort_msg(data)

        assert np.all(result == msg_fort)

    def test_unsupported_dtype_raises_error(self):
        """Test unsupported dtype raises exception."""
        # Use a structured dtype which is not supported
        data = np.array([(1, 2.0), (3, 4.0)], dtype=[("a", "i4"), ("b", "f4")])

        with pytest.raises(Exception, match="not a supported type"):
            py2fort_msg(data)

    def test_2d_array(self):
        """Test 2D array handling."""
        data = np.array([[1.0, np.nan], [3.0, 4.0]], dtype=np.float32)
        result, msg_py, msg_fort = py2fort_msg(data)

        assert result.shape == (2, 2)
        assert result[0, 1] == msg_fort

    def test_preserves_non_missing_values(self):
        """Test non-missing values are preserved."""
        data = np.array([1.5, -999.0, 3.7], dtype=np.float64)
        result, msg_py, msg_fort = py2fort_msg(data, msg_py=-999.0)

        assert result[0] == 1.5
        assert result[2] == 3.7


class TestFort2PyMsg:
    """Test fort2py_msg function."""

    def test_float32_default_conversion(self):
        """Test float32 Fortran to Python missing value."""
        fort_msg = msg_dtype[np.float32]
        data = np.array([1.0, 2.0, fort_msg, 4.0], dtype=np.float32)
        result, msg_fort, msg_py = fort2py_msg(data)

        assert np.isnan(msg_py)
        assert msg_fort == fort_msg
        assert np.isnan(result[2])

    def test_float64_custom_python_value(self):
        """Test float64 with custom Python missing value."""
        fort_msg = msg_dtype[np.float64]
        data = np.array([1.0, fort_msg, 3.0], dtype=np.float64)
        result, msg_fort, msg_py = fort2py_msg(data, msg_py=-999.0)

        assert msg_py == -999.0
        assert result[1] == -999.0

    def test_int32_conversion(self):
        """Test int32 Fortran to Python conversion."""
        fort_msg = msg_dtype[np.int32]
        data = np.array([1, fort_msg, 3], dtype=np.int32)
        result, msg_fort, msg_py = fort2py_msg(data)

        assert msg_py == fort_msg
        assert result[1] == msg_py

    def test_int16_explicit_fortran_value(self):
        """Test int16 with explicit Fortran value."""
        custom_fort = np.int16(9999)
        data = np.array([1, custom_fort, 3], dtype=np.int16)
        result, msg_fort, msg_py = fort2py_msg(data, msg_fort=custom_fort)

        assert msg_fort == custom_fort
        assert result[1] == msg_dtype[np.int16]

    def test_uint8_array(self):
        """Test uint8 array."""
        fort_msg = msg_dtype[np.uint8]
        data = np.array([10, fort_msg, 30], dtype=np.uint8)
        result, msg_fort, msg_py = fort2py_msg(data)

        assert result[1] == msg_py

    def test_complex64_conversion(self):
        """Test complex64 conversion."""
        fort_msg = msg_dtype[np.complex64]
        data = np.array([1 + 2j, fort_msg, 3 + 4j], dtype=np.complex64)
        result, msg_fort, msg_py = fort2py_msg(data)

        assert np.isnan(msg_py)
        assert np.isnan(result[1])

    def test_complex128_custom_python(self):
        """Test complex128 with custom Python value."""
        fort_msg = msg_dtype[np.complex128]
        data = np.array([1 + 1j, fort_msg, 2 + 2j], dtype=np.complex128)
        result, msg_fort, msg_py = fort2py_msg(data, msg_py=-999 - 999j)

        assert msg_py == -999 - 999j
        assert result[1] == -999 - 999j

    def test_no_fortran_missing_values(self):
        """Test array with no Fortran missing values."""
        data = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        result, msg_fort, msg_py = fort2py_msg(data)

        # Should not modify data
        assert np.allclose(result, data, equal_nan=True)

    def test_all_fortran_missing_values(self):
        """Test array with all Fortran missing values."""
        fort_msg = msg_dtype[np.float64]
        data = np.array([fort_msg, fort_msg, fort_msg], dtype=np.float64)
        result, msg_fort, msg_py = fort2py_msg(data)

        assert np.all(np.isnan(result))

    def test_unsupported_dtype_raises_error(self):
        """Test unsupported dtype raises exception."""
        data = np.array([(1, 2.0), (3, 4.0)], dtype=[("a", "i4"), ("b", "f4")])

        with pytest.raises(Exception, match="not a supported type"):
            fort2py_msg(data)

    def test_2d_array_conversion(self):
        """Test 2D array conversion."""
        fort_msg = msg_dtype[np.float32]
        data = np.array([[1.0, fort_msg], [3.0, 4.0]], dtype=np.float32)
        result, msg_fort, msg_py = fort2py_msg(data)

        assert result.shape == (2, 2)
        assert np.isnan(result[0, 1])

    def test_preserves_non_missing_values(self):
        """Test non-missing values are preserved."""
        fort_msg = msg_dtype[np.float64]
        data = np.array([1.5, fort_msg, 3.7], dtype=np.float64)
        result, msg_fort, msg_py = fort2py_msg(data, msg_py=-999.0)

        assert result[0] == 1.5
        assert result[2] == 3.7


class TestRoundTripConversion:
    """Test round-trip Python -> Fortran -> Python conversion."""

    def test_float32_roundtrip(self):
        """Test float32 round-trip preserves data."""
        original = np.array([1.0, np.nan, 3.0, np.nan, 5.0], dtype=np.float32)

        # Python to Fortran
        fort_data, msg_py1, msg_fort1 = py2fort_msg(original.copy())

        # Fortran back to Python
        py_data, msg_fort2, msg_py2 = fort2py_msg(fort_data)

        # Check non-NaN values are preserved
        assert np.allclose(py_data[~np.isnan(py_data)], [1.0, 3.0, 5.0])
        # Check NaN positions are preserved
        assert np.isnan(py_data[1])
        assert np.isnan(py_data[3])

    def test_int64_roundtrip(self):
        """Test int64 round-trip."""
        original = np.array([10, -999, 30], dtype=np.int64)

        fort_data, _, _ = py2fort_msg(original.copy(), msg_py=-999)
        py_data, _, _ = fort2py_msg(fort_data)

        assert py_data[0] == 10
        assert py_data[2] == 30

    def test_complex128_roundtrip(self):
        """Test complex128 round-trip."""
        original = np.array([1 + 2j, np.nan + np.nan * 1j, 3 + 4j], dtype=np.complex128)

        fort_data, _, _ = py2fort_msg(original.copy())
        py_data, _, _ = fort2py_msg(fort_data)

        assert py_data[0] == 1 + 2j
        assert np.isnan(py_data[1])
        assert py_data[2] == 3 + 4j


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_array(self):
        """Test empty array handling."""
        data = np.array([], dtype=np.float32)
        result, _, _ = py2fort_msg(data)
        assert len(result) == 0

    def test_single_element_array(self):
        """Test single element array."""
        data = np.array([np.nan], dtype=np.float64)
        result, _, _ = py2fort_msg(data)
        assert result[0] == msg_dtype[np.float64]

    def test_large_array_performance(self):
        """Test large array doesn't crash."""
        data = np.random.randn(10000).astype(np.float32)
        data[::100] = np.nan
        result, _, _ = py2fort_msg(data)
        assert len(result) == 10000

    def test_float16_handling(self):
        """Test float16 (half precision) handling."""
        data = np.array([1.0, np.nan, 3.0], dtype=np.float16)
        result, msg_py, msg_fort = py2fort_msg(data)

        assert np.isnan(msg_py)
        assert result[1] == msg_fort


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
