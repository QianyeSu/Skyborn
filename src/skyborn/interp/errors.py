class Error(Exception):
    """Base class for exceptions in this module."""

    pass


class AttributeError(Error):
    """Exception raised when the arguments of Skyborn interpolation functions
    have a mismatch of attributes with other arguments."""

    pass


class ChunkError(Error):
    """Exception raised when a Dask array is chunked in a way that is
    incompatible with a Skyborn interpolation function."""

    pass


class CoordinateError(Error):
    """Exception raised when a Skyborn interpolation function is passed a NumPy
    array as an argument without a required coordinate array being passed
    separately."""

    pass


class DimensionError(Error):
    """Exception raised when the arguments of Skyborn interpolation functions
    have a mismatch of the necessary dimensionality."""

    pass


class MetaError(Error):
    """Exception raised when the support for the retention of metadata is not
    supported."""

    pass
