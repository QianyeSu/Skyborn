"""Compatibility namespace for the optional ``skyborn-cdo`` package."""

try:
    from skyborn_cdo import Cdo, CdoError, get_cdo_path, get_cdo_version
except ImportError as exc:
    raise ImportError(
        "skyborn.cdo requires the optional package 'skyborn-cdo'. "
        "Install it with `pip install skyborn-cdo` or "
        '`pip install "skyborn[cdo]"` on supported platforms '
        "(Windows x86_64, Linux x86_64, macOS arm64)."
    ) from exc


__all__ = ["Cdo", "CdoError", "get_cdo_path", "get_cdo_version"]
