import numpy as np

# all missing values are represented by the maximum value in their dtype,
# int8 would be bin_01111111 or dec_127
# unint8 would be bin_11111111 or dec_255
# floats and complex types use max per IEEE-754

msg_dtype = {
    np.complex64: np.complex64(np.finfo(np.complex64).max),
    np.complex128: np.complex128(np.finfo(np.complex128).max),
    np.float16: np.float16(np.finfo(np.float16).max),
    np.float32: np.float32(np.finfo(np.float32).max),
    np.float64: np.float64(np.finfo(np.float64).max),
    # np.float128 removed - not available in Windows numpy
    np.int8: np.int8(np.iinfo(np.int8).max),
    np.int16: np.int16(np.iinfo(np.int16).max),
    np.int32: np.int32(np.iinfo(np.int32).max),
    np.int64: np.int64(np.iinfo(np.int64).max),
    np.uint8: np.uint8(np.iinfo(np.uint8).max),
    np.uint16: np.uint16(np.iinfo(np.uint16).max),
    np.uint32: np.uint32(np.iinfo(np.uint32).max),
    np.uint64: np.uint64(np.iinfo(np.uint64).max),
    str: str(""),
}

# lists of classes of dtypes
complex_dtypes = [np.complex64, np.complex128]
float_dtypes = [np.float16, np.float32, np.float64]  # float128 removed
int_dtypes = [np.int8, np.int16, np.int32, np.int64]
uint_dtypes = [np.uint8, np.uint16, np.uint32, np.uint64]
string_dtypes = [str]
supported_dtypes = msg_dtype.keys()


# python to fortran
def py2fort_msg(ndarray, msg_py=None, msg_fort=None):
    msg_indices = None
    ndtype = ndarray.dtype.type

    if ndtype not in supported_dtypes:
        raise Exception(
            "The ndarray.dtype.type of "
            + np.dtype(ndtype).name
            + " is not a supported type"
        )

    if msg_py is None:
        if ndtype in float_dtypes:
            msg_py = np.nan
        elif ndtype in complex_dtypes:
            msg_py = np.nan + np.nan * 1j
        else:
            msg_py = msg_dtype[ndtype]

    if msg_fort is None:
        msg_fort = msg_dtype[ndtype]

    if np.isnan(msg_py):
        msg_indices = np.isnan(ndarray)
    else:
        msg_indices = ndarray == msg_py

    if msg_indices.any():
        ndarray[msg_indices] = msg_fort

    return ndarray, msg_py, msg_fort


# todo: Should we force output missing value to (1) always be np.nan or (2) whatever it was given in input
#      Skyborn Fortran-backed functions currently implement (2), while some legacy tools implement (1)
def fort2py_msg(ndarray, msg_fort=None, msg_py=None):
    msg_indices = None
    ndtype = ndarray.dtype.type

    if ndtype not in supported_dtypes:
        raise Exception(
            "The ndarray.dtype.type of "
            + np.dtype(ndtype).name
            + " is not a supported type"
        )

    if msg_fort is None:
        msg_fort = msg_dtype[ndtype]

    if msg_py is None:
        if ndtype in float_dtypes:
            msg_py = np.nan
        elif ndtype in complex_dtypes:
            msg_py = np.nan + np.nan * 1j
        else:
            msg_py = msg_dtype[ndtype]

    msg_indices = ndarray == msg_fort

    if msg_indices.any():
        ndarray[msg_indices] = msg_py

    return ndarray, msg_fort, msg_py
