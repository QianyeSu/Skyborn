"""
Tests for skyborn.conversion.grib_to_netcdf module.

This module contains comprehensive tests for the GRIB to NetCDF conversion
functionality, including both unit tests and integration tests.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np

from skyborn.conversion.grib_to_netcdf import (
    convert_grib_to_nc,
    convert_grib_to_nc_simple,
    batch_convert_grib_to_nc,
    _check_grib_to_netcdf_available,
    _validate_grib_files,
    _build_grib_to_netcdf_command,
    GribToNetCDFError,
)


class TestGribToNetCDFValidation:
    """Test validation functions for GRIB to NetCDF conversion."""

    def test_check_grib_to_netcdf_availability(self):
        """Test checking if grib_to_netcdf command is available."""
        # This test checks if the function returns a boolean
        result = _check_grib_to_netcdf_available()
        assert isinstance(result, bool)

    def test_validate_grib_files_success(self, temp_dir):
        """Test validating existing GRIB files."""
        # Create temporary GRIB files
        grib_file1 = temp_dir / "test1.grib"
        grib_file2 = temp_dir / "test2.grib"

        # Create dummy files
        grib_file1.write_text("dummy grib content")
        grib_file2.write_text("dummy grib content")

        # Test validation
        result = _validate_grib_files([grib_file1, grib_file2])

        assert len(result) == 2
        assert all(isinstance(path, Path) for path in result)
        assert grib_file1 in result
        assert grib_file2 in result

    def test_validate_grib_files_not_found(self):
        """Test validating non-existent GRIB files."""
        with pytest.raises(FileNotFoundError, match="GRIB file.*not found"):
            _validate_grib_files(["nonexistent.grib"])

    def test_validate_grib_files_mixed(self, temp_dir):
        """Test validating mix of existing and non-existing files."""
        # Create one existing file
        existing_file = temp_dir / "existing.grib"
        existing_file.write_text("dummy content")

        # Test with mix of existing and non-existing
        with pytest.raises(FileNotFoundError):
            _validate_grib_files([existing_file, "nonexistent.grib"])


class TestGribToNetCDFCommandBuilder:
    """Test command building for grib_to_netcdf."""

    def test_build_basic_command(self):
        """Test building basic grib_to_netcdf command."""
        cmd = _build_grib_to_netcdf_command(
            output_file="output.nc", grib_files=["input.grib"]
        )

        expected_parts = [
            "grib_to_netcdf",
            "-D",
            "NC_SHORT",
            "-k",
            "2",
            "-o",
            "output.nc",
            "input.grib",
        ]

        for part in expected_parts:
            assert part in cmd, f"Expected '{part}' in command: {' '.join(cmd)}"

    def test_build_advanced_command(self):
        """Test building advanced grib_to_netcdf command with all options."""
        cmd = _build_grib_to_netcdf_command(
            output_file="output.nc",
            grib_files=["input1.grib", "input2.grib"],
            ignore_keys=["type", "step"],
            split_keys=["param"],
            data_type="NC_FLOAT",
            file_kind=4,
            deflate_level=6,
            shuffle=True,
            unlimited_dimension="time",
            force=True,
        )

        expected_parts = [
            "grib_to_netcdf",
            "-I",
            "type,step",
            "-S",
            "param",
            "-D",
            "NC_FLOAT",
            "-f",
            "-k",
            "4",
            "-d",
            "6",
            "-s",
            "-u",
            "time",
            "-o",
            "output.nc",
            "input1.grib",
            "input2.grib",
        ]

        for part in expected_parts:
            assert part in cmd, f"Expected '{part}' in command: {' '.join(cmd)}"

    def test_build_command_invalid_deflate_level(self):
        """Test command building with invalid deflate level."""
        with pytest.raises(ValueError, match="Deflate level must be between 0 and 9"):
            _build_grib_to_netcdf_command(
                output_file="output.nc", grib_files=["input.grib"], deflate_level=10
            )

        with pytest.raises(ValueError, match="Deflate level must be between 0 and 9"):
            _build_grib_to_netcdf_command(
                output_file="output.nc", grib_files=["input.grib"], deflate_level=-1
            )

    def test_build_command_data_types(self):
        """Test command building with different data types."""
        valid_types = ["NC_BYTE", "NC_SHORT", "NC_INT", "NC_FLOAT", "NC_DOUBLE"]

        for data_type in valid_types:
            cmd = _build_grib_to_netcdf_command(
                output_file="output.nc", grib_files=["input.grib"], data_type=data_type
            )
            assert data_type in cmd

    def test_build_command_file_kinds(self):
        """Test command building with different file kinds."""
        valid_kinds = [1, 2, 3, 4]

        for kind in valid_kinds:
            cmd = _build_grib_to_netcdf_command(
                output_file="output.nc", grib_files=["input.grib"], file_kind=kind
            )
            assert str(kind) in cmd


class TestGribToNetCDFConversion:
    """Test the main conversion functions."""

    @patch("skyborn.conversion.grib_to_netcdf._check_grib_to_netcdf_available")
    def test_convert_tool_not_available(self, mock_check):
        """Test conversion when grib_to_netcdf tool is not available."""
        mock_check.return_value = False

        with pytest.raises(GribToNetCDFError, match="grib_to_netcdf command not found"):
            convert_grib_to_nc("input.grib", "output.nc")

    @patch("skyborn.conversion.grib_to_netcdf.subprocess.run")
    @patch("skyborn.conversion.grib_to_netcdf._check_grib_to_netcdf_available")
    @patch("skyborn.conversion.grib_to_netcdf._validate_grib_files")
    def test_convert_success(self, mock_validate, mock_check, mock_run, temp_dir):
        """Test successful GRIB to NetCDF conversion."""
        # Setup mocks
        mock_check.return_value = True
        mock_validate.return_value = [Path("input.grib")]
        mock_run.return_value = MagicMock(
            returncode=0, stdout="Conversion successful", stderr=""
        )

        output_file = temp_dir / "output.nc"
        result = convert_grib_to_nc("input.grib", str(output_file), verbose=False)

        assert isinstance(result, Path)
        assert result.name == "output.nc"
        mock_run.assert_called_once()

    @patch("skyborn.conversion.grib_to_netcdf.subprocess.run")
    @patch("skyborn.conversion.grib_to_netcdf._check_grib_to_netcdf_available")
    @patch("skyborn.conversion.grib_to_netcdf._validate_grib_files")
    def test_convert_failure(self, mock_validate, mock_check, mock_run):
        """Test conversion failure handling."""
        # Setup mocks
        mock_check.return_value = True
        mock_validate.return_value = [Path("input.grib")]

        # Mock subprocess failure
        from subprocess import CalledProcessError

        mock_run.side_effect = CalledProcessError(
            returncode=1, cmd=["grib_to_netcdf"], stderr="Error: Invalid GRIB file"
        )

        with pytest.raises(GribToNetCDFError, match="grib_to_netcdf failed"):
            convert_grib_to_nc("input.grib", "output.nc", verbose=False)

    @patch("skyborn.conversion.grib_to_netcdf.convert_grib_to_nc")
    def test_convert_simple_interface(self, mock_convert):
        """Test the simplified conversion interface."""
        mock_convert.return_value = Path("output.nc")

        result = convert_grib_to_nc_simple(
            "input.grib", "output.nc", high_precision=True, compress=True
        )

        # Verify the underlying function was called with correct parameters
        mock_convert.assert_called_once()
        call_args = mock_convert.call_args

        # Check that high_precision=True resulted in NC_FLOAT
        assert call_args[1]["data_type"] == "NC_FLOAT"

        # Check that compress=True resulted in appropriate settings
        assert call_args[1]["file_kind"] == 4
        assert call_args[1]["deflate_level"] == 6
        assert call_args[1]["shuffle"] is True

    @patch("skyborn.conversion.grib_to_netcdf.convert_grib_to_nc")
    def test_convert_simple_default_options(self, mock_convert):
        """Test simplified interface with default options."""
        mock_convert.return_value = Path("output.nc")

        result = convert_grib_to_nc_simple("input.grib", "output.nc")

        mock_convert.assert_called_once()
        call_args = mock_convert.call_args

        # Check defaults
        assert call_args[1]["data_type"] == "NC_SHORT"  # high_precision=False
        assert call_args[1]["file_kind"] == 2  # compress=False


class TestBatchConversion:
    """Test batch conversion functionality."""

    def test_batch_convert_no_files(self, temp_dir):
        """Test batch conversion when no GRIB files are found."""
        input_dir = temp_dir / "input"
        output_dir = temp_dir / "output"
        input_dir.mkdir()

        with pytest.raises(FileNotFoundError, match="No GRIB files found"):
            batch_convert_grib_to_nc(input_dir, output_dir)

    @patch("skyborn.conversion.grib_to_netcdf.convert_grib_to_nc")
    def test_batch_convert_success(self, mock_convert, temp_dir):
        """Test successful batch conversion."""
        # Setup directories
        input_dir = temp_dir / "input"
        output_dir = temp_dir / "output"
        input_dir.mkdir()
        output_dir.mkdir()

        # Create dummy GRIB files
        grib_files = ["file1.grib", "file2.grib", "file3.grib"]
        for filename in grib_files:
            (input_dir / filename).write_text("dummy grib content")

        # Mock successful conversions
        def mock_convert_side_effect(input_file, output_file, **kwargs):
            return Path(output_file)

        mock_convert.side_effect = mock_convert_side_effect

        # Run batch conversion
        results = batch_convert_grib_to_nc(input_dir, output_dir, pattern="*.grib")

        # Check results
        assert len(results) == 3
        assert mock_convert.call_count == 3

        # Verify all files were processed
        processed_files = [call[0][0] for call in mock_convert.call_args_list]
        expected_files = [str(input_dir / f) for f in grib_files]

        for expected in expected_files:
            assert any(expected in str(processed) for processed in processed_files)

    def test_batch_convert_custom_pattern(self, temp_dir):
        """Test batch conversion with custom file pattern."""
        input_dir = temp_dir / "input"
        output_dir = temp_dir / "output"
        input_dir.mkdir()

        # Create files with different extensions
        (input_dir / "file1.grib").write_text("content")
        (input_dir / "file2.grb").write_text("content")
        (input_dir / "file3.txt").write_text("content")

        # Should find no .grb2 files
        with pytest.raises(FileNotFoundError):
            batch_convert_grib_to_nc(input_dir, output_dir, pattern="*.grb2")


class TestConversionIntegration:
    """Integration tests for the conversion module."""

    def test_conversion_error_types(self):
        """Test that proper error types are raised."""
        # Test GribToNetCDFError is properly defined
        assert issubclass(GribToNetCDFError, Exception)

        # Test creating error instances
        error = GribToNetCDFError("Test error message")
        assert str(error) == "Test error message"

    def test_conversion_with_pathlib(self, temp_dir):
        """Test conversion functions work with pathlib.Path objects."""
        grib_file = temp_dir / "input.grib"
        netcdf_file = temp_dir / "output.nc"

        # Create dummy file
        grib_file.write_text("dummy content")

        # Test that Path objects are handled correctly in validation
        result = _validate_grib_files([grib_file])
        assert len(result) == 1
        assert isinstance(result[0], Path)

    @pytest.mark.slow
    def test_conversion_performance(self):
        """Test conversion functions handle reasonable file counts."""
        # Test command building with many files
        many_files = [f"file{i}.grib" for i in range(100)]

        cmd = _build_grib_to_netcdf_command(
            output_file="output.nc", grib_files=many_files
        )

        # Should complete without issues
        assert "grib_to_netcdf" in cmd
        assert "output.nc" in cmd
        assert len([f for f in cmd if f.endswith(".grib")]) == 100


class TestConversionErrorHandling:
    """Test comprehensive error handling."""

    def test_invalid_input_types(self):
        """Test handling of invalid input types."""
        with pytest.raises(TypeError):
            _validate_grib_files("not_a_list")

        with pytest.raises(TypeError):
            _validate_grib_files([123, 456])  # Numbers instead of paths

    def test_empty_file_lists(self):
        """Test handling of empty file lists."""
        result = _validate_grib_files([])
        assert result == []

        # Empty grib_files should raise error in command building
        with pytest.raises((ValueError, IndexError)):
            _build_grib_to_netcdf_command(output_file="output.nc", grib_files=[])


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
