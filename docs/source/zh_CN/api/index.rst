API 参考
========

本节包含所有 Skyborn 模块的完整 API 文档。

.. toctree::
   :maxdepth: 2

   conversion
   calculations
   gradients
   causality
   interpolation
   plotting

核心函数
--------

.. currentmodule:: skyborn

主要函数
~~~~~~~~

.. autosummary::
   :toctree: generated/

   convert_grib_to_nc
   convert_grib_to_nc_simple
   batch_convert_grib_to_nc
   grib2nc
   grib_to_netcdf

数据处理
~~~~~~~~

.. autosummary::
   :toctree: generated/

   convert_longitude_range
   linear_regression

梯度计算
~~~~~~~~

.. autosummary::
   :toctree: generated/

   calculate_gradient
   calculate_meridional_gradient
   calculate_zonal_gradient
   calculate_vertical_gradient

因果关系分析
~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   liang_causality
   granger_causality
