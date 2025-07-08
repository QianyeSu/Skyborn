.. Skyborn 文档主文件

欢迎使用 Skyborn 文档！
=======================

Skyborn 是一个专为大气数据处理、分析和可视化设计的综合性 Python 库。
它提供了处理气象和气候数据集的工具，包括 GRIB 到 NetCDF 转换、
数据插值、梯度计算和专业绘图功能。

.. toctree::
   :maxdepth: 2
   :caption: 目录:
   :numbered:

   installation
   quickstart
   api/index
   modules/index
   examples/index
   changelog
   contributing

.. toctree::
   :maxdepth: 1
   :caption: 语言选择:

   中文 <self>
   English <../en/index>

主要功能
========

* **数据转换**: 使用 eccodes 进行 GRIB 到 NetCDF 转换
* **插值处理**: 高级插值和重网格化功能
* **梯度计算**: 经向、纬向和垂直梯度计算
* **因果分析**: Liang 和 Granger 因果关系分析
* **数据可视化**: 专业的大气数据绘图工具
* **ROF 分析**: 河流流量和水文分析工具

快速示例
========

.. code-block:: python

   import skyborn

   # GRIB 转换为 NetCDF
   skyborn.grib2nc('input.grib', 'output.nc')

   # 计算梯度
   import xarray as xr
   ds = xr.open_dataset('output.nc')
   grad = skyborn.calculate_gradient(ds['temperature'], 'latitude')

获取帮助
========

* **GitHub Issues**: 在 `GitHub Issues <https://github.com/yourusername/skyborn/issues>`_ 报告错误和请求功能
* **文档**: 完整的 API 文档和示例
* **邮箱**: 联系作者 suqianye2000@gmail.com

索引和表格
==========

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
