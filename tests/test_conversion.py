"""
Test skyborn.conversion module GRIB to NetCDF conversion functionality

This file contains basic tests for the conversion module.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock

from skyborn.conversion import (
    convert_grib_to_nc,
    convert_grib_to_nc_simple,
    batch_convert_grib_to_nc,
    _check_grib_to_netcdf_available,
    _validate_grib_files,
    _build_grib_to_netcdf_command,
    GribToNetCDFError,
)


# class TestGribToNetCDF:
#     """Test GRIB to NetCDF conversion functionality"""

#     def test_check_grib_to_netcdf_available(self):
#         """Test checking if grib_to_netcdf is available"""
#         # This test depends on whether eccodes is installed on the system
#         result = _check_grib_to_netcdf_available()
#         assert isinstance(result, bool)

#     def test_validate_grib_files_success(self):
#         """Test validating existing GRIB files"""
#         with tempfile.NamedTemporaryFile(suffix=".grib", delete=False) as tmp:
#             tmp_path = Path(tmp.name)

#         try:
#             result = _validate_grib_files([tmp_path])
#             assert len(result) == 1
#             assert isinstance(result[0], Path)
#         finally:
#             tmp_path.unlink()

#     def test_validate_grib_files_not_found(self):
#         """Test validating non-existent GRIB files"""
#         with pytest.raises(FileNotFoundError):
#             _validate_grib_files(["nonexistent.grib"])

#     def test_build_grib_to_netcdf_command_basic(self):
#         """Test building basic command"""
#         cmd = _build_grib_to_netcdf_command(
#             output_file="output.nc", grib_files=["input.grib"]
#         )

#         expected_parts = [
#             "grib_to_netcdf",
#             "-D",
#             "NC_SHORT",
#             "-k",
#             "2",
#             "-o",
#             "output.nc",
#             "input.grib",
#         ]

#         for part in expected_parts:
#             assert part in cmd

#     def test_build_grib_to_netcdf_command_advanced(self):
#         """Test building advanced command"""
#         cmd = _build_grib_to_netcdf_command(
#             output_file="output.nc",
#             grib_files=["input1.grib", "input2.grib"],
#             ignore_keys=["type", "step"],
#             split_keys=["param"],
#             data_type="NC_FLOAT",
#             file_kind=4,
#             deflate_level=6,
#             shuffle=True,
#             unlimited_dimension="time",
#             force=True,
#         )

#         expected_parts = [
#             "grib_to_netcdf",
#             "-I",
#             "type,step",
#             "-S",
#             "param",
#             "-D",
#             "NC_FLOAT",
#             "-f",
#             "-k",
#             "4",
#             "-d",
#             "6",
#             "-s",
#             "-u",
#             "time",
#             "-o",
#             "output.nc",
#             "input1.grib",
#             "input2.grib",
#         ]

#         for part in expected_parts:
#             assert part in cmd

#     def test_build_grib_to_netcdf_command_invalid_deflate(self):
#         """Test invalid compression level"""
#         with pytest.raises(ValueError, match="Deflate level must be between 0 and 9"):
#             _build_grib_to_netcdf_command(
#                 output_file="output.nc", grib_files=["input.grib"], deflate_level=10
#             )

#     @patch("skyborn.conversion._check_grib_to_netcdf_available")
#     def test_convert_grib_to_nc_tool_not_available(self, mock_check):
#         """Test when grib_to_netcdf tool is not available"""
#         mock_check.return_value = False

#         with pytest.raises(GribToNetCDFError, match="grib_to_netcdf command not found"):
#             convert_grib_to_nc("input.grib", "output.nc")

#     @patch("skyborn.conversion.subprocess.run")
#     @patch("skyborn.conversion._check_grib_to_netcdf_available")
#     @patch("skyborn.conversion._validate_grib_files")
#     def test_convert_grib_to_nc_success(self, mock_validate, mock_check, mock_run):
#         """Test successful conversion"""
#         # Setup mocks
#         mock_check.return_value = True
#         mock_validate.return_value = [Path("input.grib")]
#         mock_run.return_value = MagicMock(
#             returncode=0, stdout="Conversion successful", stderr=""
#         )

#         result = convert_grib_to_nc("input.grib", "output.nc", verbose=False)

#         assert isinstance(result, Path)
#         assert result.name == "output.nc"
#         mock_run.assert_called_once()

#     @patch("skyborn.conversion.subprocess.run")
#     @patch("skyborn.conversion._check_grib_to_netcdf_available")
#     @patch("skyborn.conversion._validate_grib_files")
#     def test_convert_grib_to_nc_failure(self, mock_validate, mock_check, mock_run):
#         """Test conversion failure"""
#         # Setup mocks
#         mock_check.return_value = True
#         mock_validate.return_value = [Path("input.grib")]

#         # Mock command failure
#         from subprocess import CalledProcessError

#         mock_run.side_effect = CalledProcessError(
#             returncode=1, cmd=["grib_to_netcdf"], stderr="Error message"
#         )

#         with pytest.raises(GribToNetCDFError, match="grib_to_netcdf failed"):
#             convert_grib_to_nc("input.grib", "output.nc", verbose=False)

#     @patch("skyborn.conversion.convert_grib_to_nc")
#     def test_convert_grib_to_nc_simple(self, mock_convert):
#         """Test simplified interface"""
#         mock_convert.return_value = Path("output.nc")

#         result = convert_grib_to_nc_simple(
#             "input.grib", "output.nc", high_precision=True, compress=True
#         )

#         # Verify call arguments
#         mock_convert.assert_called_once()
#         call_args = mock_convert.call_args

#         assert call_args[1]["data_type"] == "NC_FLOAT"  # high_precision=True
#         assert call_args[1]["file_kind"] == 4  # compress=True
#         assert call_args[1]["deflate_level"] == 6  # compress=True
#         assert call_args[1]["shuffle"] is True  # compress=True

#     def test_batch_convert_no_files(self):
#         """Test batch conversion when no files are found"""
#         with tempfile.TemporaryDirectory() as tmpdir:
#             input_dir = Path(tmpdir) / "input"
#             output_dir = Path(tmpdir) / "output"
#             input_dir.mkdir()

#             with pytest.raises(FileNotFoundError, match="No GRIB files found"):
#                 batch_convert_grib_to_nc(input_dir, output_dir)


def test_module_imports():
    """Test module imports"""
    import skyborn

    # Test direct import from skyborn
    assert hasattr(skyborn, "convert_grib_to_nc")
    assert hasattr(skyborn, "convert_grib_to_nc_simple")
    assert hasattr(skyborn, "batch_convert_grib_to_nc")
    assert hasattr(skyborn, "grib2nc")
    assert hasattr(skyborn, "grib_to_netcdf")

    # Test import from skyborn.conversion
    from skyborn.conversion import convert_grib_to_nc

    assert callable(convert_grib_to_nc)


if __name__ == "__main__":
    # Run simple tests
    print("Running skyborn.conversion module tests...")

    # Test module imports
    test_module_imports()
    print("✅ Module import tests passed")

    # Test tool checking
    # available = _check_grib_to_netcdf_available()
    # print(f"grib_to_netcdf tool availability: {available}")

    # if available:
    #     print("✅ grib_to_netcdf tool is available, actual conversion is possible")
    # else:
    #     print("⚠️  grib_to_netcdf tool is not available, please install eccodes")

    print("Tests completed!")
