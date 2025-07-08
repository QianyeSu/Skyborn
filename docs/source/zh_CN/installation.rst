安装指南
========

系统要求
--------

Skyborn 需要 Python 3.8 或更高版本。该库依赖于几个科学计算 Python 包：

核心依赖
~~~~~~~~

* **NumPy** (>=1.19.0) - 数值计算
* **Pandas** (>=1.3.0) - 数据处理
* **Xarray** (>=0.19.0) - N 维标记数组
* **SciPy** (>=1.7.0) - 科学计算
* **Matplotlib** (>=3.4.0) - 绘图

可选依赖
~~~~~~~~

* **eccodes** (>=1.4.0) - 用于 GRIB 到 NetCDF 转换
* **Cartopy** (>=0.20.0) - 用于地理空间绘图
* **Dask** (>=2021.6.0) - 用于并行计算

安装方法
--------

从 PyPI 安装（推荐）
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install skyborn

这将安装 Skyborn 及其核心依赖。

从 Conda 安装
~~~~~~~~~~~~~

.. code-block:: bash

   conda install -c conda-forge skyborn

从源码安装
~~~~~~~~~~

获取最新开发版本：

.. code-block:: bash

   git clone https://github.com/yourusername/skyborn.git
   cd skyborn
   pip install -e .

安装可选依赖
------------

GRIB 转换功能
~~~~~~~~~~~~~

要使用 GRIB 到 NetCDF 转换功能：

.. code-block:: bash

   # 使用 conda（推荐）
   conda install -c conda-forge eccodes

   # 使用 pip
   pip install eccodes

地理空间绘图
~~~~~~~~~~~~

.. code-block:: bash

   # 使用 conda
   conda install -c conda-forge cartopy

   # 使用 pip
   pip install cartopy

完整安装
~~~~~~~~

安装所有可选依赖：

.. code-block:: bash

   # 使用 conda
   conda install -c conda-forge skyborn eccodes cartopy dask

   # 使用 pip
   pip install skyborn[complete]

验证安装
--------

验证安装是否成功：

.. code-block:: python

   import skyborn
   print(f"Skyborn 版本: {skyborn.__version__}")

   # 测试基本功能
   skyborn.convert_longitude_range  # 不应该报错

   # 测试 GRIB 转换（如果安装了 eccodes）
   try:
       skyborn.grib2nc
       print("GRIB 转换功能可用 ✓")
   except AttributeError:
       print("GRIB 转换功能不可用（未安装 eccodes）")

故障排除
--------

常见问题
~~~~~~~~

1. **导入错误**: 确保所有依赖已安装
2. **GRIB 转换错误**: 安装 eccodes 库
3. **绘图问题**: 安装 matplotlib 和可选的 cartopy

平台特定说明
~~~~~~~~~~~~

**Windows**:
  - 推荐使用 Anaconda/Miniconda
  - 某些包可能需要 Visual Studio Build Tools

**macOS**:
  - 可能需要安装 Xcode 命令行工具
  - 使用 conda 更容易管理依赖

**Linux**:
  - 通常开箱即用
  - 某些包可能需要开发头文件

开发者安装
----------

对于开发者：

.. code-block:: bash

   git clone https://github.com/yourusername/skyborn.git
   cd skyborn
   pip install -e .[dev]

   # 安装 pre-commit 钩子
   pre-commit install

这将安装额外的开发依赖，包括测试和文档工具。
