"""
Comprehensive tests for skyborn.interp.errors module.

Tests all custom exception classes to ensure proper error handling.
Target: 100% coverage for errors.py
"""

import pytest

from skyborn.interp.errors import ChunkError, CoordinateError, DimensionError


class TestChunkError:
    """Test ChunkError exception class."""

    def test_chunk_error_instantiation(self):
        """Test ChunkError can be instantiated."""
        error = ChunkError("Test chunk error message")
        assert isinstance(error, Exception)
        assert str(error) == "Test chunk error message"

    def test_chunk_error_raise(self):
        """Test ChunkError can be raised and caught."""
        with pytest.raises(ChunkError) as exc_info:
            raise ChunkError("Chunking issue detected")
        assert "Chunking issue detected" in str(exc_info.value)

    def test_chunk_error_inheritance(self):
        """Test ChunkError inherits from Exception."""
        error = ChunkError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, ChunkError)

    def test_chunk_error_empty_message(self):
        """Test ChunkError with empty message."""
        error = ChunkError("")
        assert str(error) == ""

    def test_chunk_error_multiline_message(self):
        """Test ChunkError with multiline message."""
        msg = "Line 1\nLine 2\nLine 3"
        error = ChunkError(msg)
        assert str(error) == msg


class TestCoordinateError:
    """Test CoordinateError exception class."""

    def test_coordinate_error_instantiation(self):
        """Test CoordinateError can be instantiated."""
        error = CoordinateError("Test coordinate error message")
        assert isinstance(error, Exception)
        assert str(error) == "Test coordinate error message"

    def test_coordinate_error_raise(self):
        """Test CoordinateError can be raised and caught."""
        with pytest.raises(CoordinateError) as exc_info:
            raise CoordinateError("Coordinate mismatch detected")
        assert "Coordinate mismatch detected" in str(exc_info.value)

    def test_coordinate_error_inheritance(self):
        """Test CoordinateError inherits from Exception."""
        error = CoordinateError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, CoordinateError)

    def test_coordinate_error_with_details(self):
        """Test CoordinateError with detailed message."""
        msg = "Expected longitude in range [-180, 180], got 270"
        error = CoordinateError(msg)
        assert "270" in str(error)

    def test_coordinate_error_unicode(self):
        """Test CoordinateError with unicode characters."""
        error = CoordinateError("坐标错误: 经度范围不正确")
        assert "坐标错误" in str(error)


class TestDimensionError:
    """Test DimensionError exception class."""

    def test_dimension_error_instantiation(self):
        """Test DimensionError can be instantiated."""
        error = DimensionError("Test dimension error message")
        assert isinstance(error, Exception)
        assert str(error) == "Test dimension error message"

    def test_dimension_error_raise(self):
        """Test DimensionError can be raised and caught."""
        with pytest.raises(DimensionError) as exc_info:
            raise DimensionError("Dimension mismatch: expected 2D, got 3D")
        assert "Dimension mismatch" in str(exc_info.value)

    def test_dimension_error_inheritance(self):
        """Test DimensionError inherits from Exception."""
        error = DimensionError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, DimensionError)

    def test_dimension_error_with_shape_info(self):
        """Test DimensionError with shape information."""
        msg = "Shape mismatch: expected (100, 200), got (100, 150)"
        error = DimensionError(msg)
        assert "(100, 200)" in str(error)
        assert "(100, 150)" in str(error)

    def test_dimension_error_formatting(self):
        """Test DimensionError with formatted string."""
        expected_shape = (10, 20)
        actual_shape = (10, 30)
        error = DimensionError(
            f"Shape mismatch: expected {expected_shape}, got {actual_shape}"
        )
        assert "10" in str(error)
        assert "20" in str(error)
        assert "30" in str(error)


class TestErrorsInteraction:
    """Test interactions between different error types."""

    def test_catch_specific_errors(self):
        """Test catching specific error types."""
        with pytest.raises(ChunkError):
            raise ChunkError("Chunk error")

        with pytest.raises(CoordinateError):
            raise CoordinateError("Coordinate error")

        with pytest.raises(DimensionError):
            raise DimensionError("Dimension error")

    def test_catch_as_exception(self):
        """Test all errors can be caught as Exception."""
        for error_class, message in [
            (ChunkError, "Chunk"),
            (CoordinateError, "Coordinate"),
            (DimensionError, "Dimension"),
        ]:
            with pytest.raises(Exception) as exc_info:
                raise error_class(message)
            assert message in str(exc_info.value)

    def test_error_type_differentiation(self):
        """Test errors can be differentiated by type."""
        chunk_err = ChunkError("chunk")
        coord_err = CoordinateError("coord")
        dim_err = DimensionError("dim")

        assert type(chunk_err) != type(coord_err)
        assert type(coord_err) != type(dim_err)
        assert type(chunk_err) != type(dim_err)

    def test_error_comparison(self):
        """Test error equality and comparison."""
        err1 = ChunkError("message")
        err2 = ChunkError("message")
        err3 = ChunkError("different")

        # Errors are distinct objects even with same message
        assert err1 is not err2
        assert str(err1) == str(err2)
        assert str(err1) != str(err3)


class TestErrorsEdgeCases:
    """Test edge cases for error handling."""

    def test_error_with_none_message(self):
        """Test errors with None as message."""
        # Python's Exception handles None gracefully
        error = ChunkError(None)
        assert error is not None

    def test_error_repr(self):
        """Test error representation."""
        error = CoordinateError("test message")
        repr_str = repr(error)
        assert "CoordinateError" in repr_str or "test message" in repr_str

    def test_multiple_raises(self):
        """Test raising same error multiple times."""
        error = DimensionError("persistent error")

        with pytest.raises(DimensionError):
            raise error

        # Can raise the same error instance again
        with pytest.raises(DimensionError):
            raise error

    def test_error_in_function(self):
        """Test errors raised from functions."""

        def raise_chunk_error():
            raise ChunkError("Function chunk error")

        def raise_coord_error():
            raise CoordinateError("Function coord error")

        def raise_dim_error():
            raise DimensionError("Function dim error")

        with pytest.raises(ChunkError):
            raise_chunk_error()

        with pytest.raises(CoordinateError):
            raise_coord_error()

        with pytest.raises(DimensionError):
            raise_dim_error()

    def test_error_context_preservation(self):
        """Test error context is preserved in exception chain."""
        try:
            try:
                raise ValueError("Original error")
            except ValueError as e:
                raise CoordinateError("Wrapped error") from e
        except CoordinateError as ce:
            assert ce.__cause__ is not None
            assert isinstance(ce.__cause__, ValueError)
            assert "Original error" in str(ce.__cause__)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
