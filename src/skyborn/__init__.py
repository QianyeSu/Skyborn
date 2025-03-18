"""
Skyborn 大气科学研究工具库

提供大气科学常用数据处理、分析和可视化功能模块
"""

from .io import read_netcdf, read_grib
from .calculations import (
    convert_longitude_range,
    linearRegression
)
from .plotting import add_equal_axes, createFigure
from .gradients import (
    calculate_gradient,
    calculate_meridional_gradient,
    calculate_zonal_gradient,
    calculate_vertical_gradient
)

from .causality import (
    liang_causality,
    granger_causality
)

__version__ = "0.2.0"
__all__ = ['io', 'calculations', 'gradients', 'plotting', 'causality']
