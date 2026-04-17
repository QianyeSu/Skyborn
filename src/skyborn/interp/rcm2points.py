import typing

import numpy as np
import xarray as xr

from .errors import ChunkError, CoordinateError, DimensionError
from .fortran.rcm2points import drcm2points
from .missing_values import (
    complex_dtypes,
    float_dtypes,
    fort2py_msg,
    msg_dtype,
    py2fort_msg,
)

supported_types = typing.Union[xr.DataArray, np.ndarray]


def _default_missing_value_for_dtype(ndtype):
    if ndtype in float_dtypes:
        return np.nan
    if ndtype in complex_dtypes:
        return np.nan + np.nan * 1j
    return msg_dtype[ndtype]


# Fortran Wrapper _<funcname>()
# This wrapper is executed within dask processes (if any), and could/should
# do anything that can benefit from parallel execution.


def _rcm2points(lat2d, lon2d, fi, lat1d, lon1d, msg_py, opt):
    fi = np.transpose(fi, axes=(2, 1, 0))
    lat2d = np.transpose(lat2d, axes=(1, 0))
    lon2d = np.transpose(lon2d, axes=(1, 0))

    ndtype = fi.dtype.type
    msg_fort = msg_dtype[ndtype]
    if msg_py is None:
        msg_py = _default_missing_value_for_dtype(ndtype)

    needs_restore = False
    if ndtype in float_dtypes or ndtype in complex_dtypes:
        msg_is_nan = bool(np.isnan(msg_py))
    else:
        msg_is_nan = False

    # Skip the generic missing-value conversion when the transposed source
    # field already matches the Fortran-side sentinel contract. This avoids an
    # extra full scan plus the restore pass on the input cube.
    if msg_is_nan:
        if np.isnan(fi).any():
            fi, msg_py, msg_fort = py2fort_msg(fi, msg_py=msg_py, msg_fort=msg_fort)
            needs_restore = True
    elif msg_py != msg_fort:
        # For explicit custom sentinels, fall back to the historical single-pass
        # conversion path instead of pre-scanning for `.any()`, which would add
        # a second full read of the input cube.
        fi, msg_py, msg_fort = py2fort_msg(fi, msg_py=msg_py, msg_fort=msg_fort)
        needs_restore = True

    fo = drcm2points(lat2d, lon2d, fi, lat1d, lon1d, xmsg=msg_fort, opt=opt)
    fo = np.asarray(fo).T

    if needs_restore:
        fort2py_msg(fi, msg_fort=msg_fort, msg_py=msg_py)
    fort2py_msg(fo, msg_fort=msg_fort, msg_py=msg_py)

    return fo


# Helper function that checks the validity of the input


def _validity_check(lat2d, lon2d, fi, lat1d, lon1d):
    if (lon2d is None) | (lat2d is None):
        raise CoordinateError("rcm2points: lon2d and lat2d should always be provided")

    if lat2d.shape[0] != lon2d.shape[0] or lat2d.shape[1] != lon2d.shape[1]:
        raise DimensionError(
            "rcm2points: The input lat/lon grids must be the same size !"
        )

    if lat1d.shape[0] != lon1d.shape[0]:
        raise DimensionError("rcm2points: The output lat/lon grids must be same size !")

    if (
        lat2d.shape[0] < 2
        or lon2d.shape[0] < 2
        or lat2d.shape[1] < 2
        or lon2d.shape[1] < 2
    ):
        raise DimensionError(
            "rcm2points: The input/output lat/lon grids must have at least 2 elements !"
        )

    if fi.ndim < 2:
        raise DimensionError("rcm2points: fi must be at least two dimensions !\n")

    if (
        fi.shape[fi.ndim - 2] != lat2d.shape[0]
        or fi.shape[fi.ndim - 1] != lon2d.shape[1]
    ):
        raise DimensionError(
            "rcm2points: The rightmost dimensions of fi must be (nlat2d x nlon2d),"
            "where nlat2d and nlon2d are the size of the lat2d/lon2d arrays !"
        )


# Outer Wrapper <funcname>()
# This wrapper is excecuted in the __main__ python process, and should be
# used for any tasks which would not benefit from parallel execution.

# TODO: Even though this function is advertised to work on multi-dimensional arrays,
#  it is currently only applicable to 3D-arrays due to the implementation of `_rcm2points`


# TODO: This function requires the input to have the coordinates in the rightmost two dimensions,
#  but xarray.DataArrray inputs with coordinates anywhere could/should actually be fine
def rcm2points(
    lat2d: supported_types,
    lon2d: supported_types,
    fi: supported_types,
    lat1d: supported_types,
    lon1d: supported_types,
    opt: np.number = 0,
    msg: np.number = None,
    meta: bool = False,
) -> supported_types:
    """Interpolates data on a curvilinear grid (i.e. RCM, WRF, NARR) to an
    unstructured grid.

    Interpolates data on a curvilinear grid, such as those used by the RCM
    (Regional Climate Model), WRF (Weather Research and Forecasting) and NARR
    (North American Regional Reanalysis) models/datasets to an unstructured grid.
    All of these have latitudes that are oriented south-to-north. An inverse
    distance squared algorithm is used to perform the interpolation. Missing
    values are allowed and no extrapolation is performed.

    Parameters
    -----------

    lat2d : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        A two-dimensional array that specifies the latitudes locations of ``fi``. The latitude
        order must be south-to-north. Because this array is two-dimensional it is not an
        associated coordinate variable of ``fi``, so it should always be explicitly provided.

    lon2d : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        A two-dimensional array that specifies the longitude locations of ``fi``. The longitude
        order must be west-to-east. Because this array is two-dimensional it is not an
        associated coordinate variable of ``fi``, so it should always be explicitly provided.

    fi : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        A multi-dimensional array to be interpolated. The rightmost two dimensions (latitude,
        longitude) are the dimensions to be interpolated.

    lat1d : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        A one-dimensional array that specifies the latitude coordinates of the output locations.

    lon1d : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        A one-dimensional array that specifies the longitude coordinates of the output locations.

    opt : :obj:`numpy.number`
        ``opt=0`` or ``1`` means use an inverse distance weight interpolation.
        ``opt=2`` means use a bilinear interpolation.

    msg : :obj:`numpy.number`
        A numpy scalar value that represent a missing value in ``fi``. This argument allows a user to
        use a missing value scheme other than NaN or masked arrays, similar to what NCL allows.

    meta : :obj:`bool`
        If set to True and the input array is an Xarray, the metadata from the input array will be
        copied to the output array; default is False.
        Warning: This option is not currently supported.

    Returns
    -------

        fo : :class:`xarray.DataArray`, :class:`numpy.ndarray`
            The interpolated grid. A multi-dimensional array of the same size as ``fi`` except that the
            rightmost dimension sizes have been replaced by the number of coordinate pairs (``lat1d``,
            ``lon1d``).
    """

    # Basic validity checks
    _validity_check(lat2d, lon2d, fi, lat1d, lon1d)

    is_input_xr = isinstance(fi, xr.DataArray)
    fi_data = fi.data if is_input_xr else np.asarray(fi)
    lat2d_data = lat2d.data if isinstance(lat2d, xr.DataArray) else np.asarray(lat2d)
    lon2d_data = lon2d.data if isinstance(lon2d, xr.DataArray) else np.asarray(lon2d)
    lat1d_data = lat1d.data if isinstance(lat1d, xr.DataArray) else np.asarray(lat1d)
    lon1d_data = lon1d.data if isinstance(lon1d, xr.DataArray) else np.asarray(lon1d)

    # If input data is already chunked
    if is_input_xr and fi.chunks is not None:
        # Ensure last two dimensions of `fi` are not chunked
        if list(fi.chunks)[-2:] != [(lat2d.shape[0],), (lat2d.shape[1],)]:
            # [(lon2d.shape[0]), (lon2d.shape[1])] could also be used
            raise ChunkError(
                "rcm2points: fi must be unchunked along the rightmost two dimensions"
            )

    fo = _rcm2points(lat2d_data, lon2d_data, fi_data, lat1d_data, lon1d_data, msg, opt)

    # If input was xarray.DataArray, convert output to xarray.DataArray as well
    if is_input_xr:
        fo = xr.DataArray(fo, attrs=fi.attrs)

    return fo
