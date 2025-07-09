# Skyborn Calc Module - Code Reorganization

## 📁 **新的代码结构 (New Code Structure)**

```
src/skyborn/
├── __init__.py                 # 主包入口
├── calc/                       # 🆕 新的计算模块文件夹
│   ├── __init__.py            # 计算模块入口
│   ├── calculations.py        # 统计计算和线性回归
│   └── emergent_constraints.py # 涌现约束方法
├── conversion/                 # GRIB转换模块
├── gradients.py               # 梯度计算
├── causality.py               # 因果关系分析
├── plot/                      # 绘图模块
├── interp/                    # 插值模块
└── ROF/                       # ROF模块
```

## 🔄 **移动的文件 (Moved Files)**

| 原位置 | 新位置 | 描述 |
|--------|--------|------|
| `src/skyborn/calculations.py` | `src/skyborn/calc/calculations.py` | 统计计算函数 |
| `emergent_constraint.py` | `src/skyborn/calc/emergent_constraints.py` | 涌现约束方法 |

## 📦 **可用函数 (Available Functions)**

### **统计计算 (Statistical Calculations)**
```python
from skyborn.calc import (
    linear_regression,           # 线性回归
    convert_longitude_range,     # 经度范围转换
    pearson_correlation,         # 皮尔逊相关系数
    spearman_correlation,        # 斯皮尔曼相关系数
    kendall_correlation,         # 肯德尔相关系数
    calculate_potential_temperature  # 位温计算
)
```

### **涌现约束方法 (Emergent Constraint Methods)**

#### **新函数名 (New Function Names - Recommended)**
```python
from skyborn.calc import (
    gaussian_pdf,                    # 高斯PDF计算 (NEW!)
    emergent_constraint_posterior,   # 涌现约束后验PDF (NEW!)
    emergent_constraint_prior,       # 涌现约束先验PDF (NEW!)
)
```

#### **传统函数名 (Legacy Function Names - For Compatibility)**
```python
from skyborn.calc import (
    calc_GAUSSIAN_PDF,          # 高斯PDF计算 (legacy)
    calc_PDF_EC,                # 涌现约束PDF计算 (legacy)
    find_std_from_PDF,          # 从PDF计算标准差 (legacy)
    calc_PDF_EC_PRIOR           # 先验概率计算 (legacy)
)
```

### **🔗 参考来源 (References)**
涌现约束方法的实现参考了以下来源：
- **GitHub**: https://github.com/blackcata/Emergent_Constraints/tree/master
- **论文**: Cox, P. M., et al. (2013). Nature, 494(7437), 341-344.

## 🎯 **使用示例 (Usage Examples)**

### **方式1：使用新的改进函数名（推荐）**
```python
import skyborn as skb
import numpy as np

# 使用新的函数名 - 更清晰明确
x_values = np.linspace(-3, 3, 100)
pdf = skb.gaussian_pdf(mu=0, sigma=1, x=x_values)

# 应用涌现约束（新函数名更明确表达用途）
posterior_pdf, sigma, mean = skb.emergent_constraint_posterior(
    constraint_data, target_data, x_grid, y_grid, obs_pdf
)
```

### **方式2：从calc子模块导入（推荐用于开发）**
```python
from skyborn.calc import (
    linear_regression,
    gaussian_pdf,  # 新名称：清晰表达功能
    emergent_constraint_posterior  # 新名称：明确是后验约束
)

# 执行线性回归
coeff, p_values = linear_regression(data_3d, predictor)

# 应用涌现约束 - 函数名现在清晰表达了用途
constrained_pdf, sigma, mean = emergent_constraint_posterior(
    tmp_x, tmp_y, x, y, obs_pdf
)
```

### **方式3：传统函数名（向后兼容）**
```python
from skyborn.calc import calc_PDF_EC  # 旧名称仍然可用

# 旧的函数名仍然工作，但建议使用新名称
pdf, sigma, mean = calc_PDF_EC(tmp_x, tmp_y, x, y, obs_pdf)
```

## 📄 **示例文件 (Example Files)**

- **`examples/calc_example.py`** - Python脚本示例
- **`examples/emergent_constraints_demo.ipynb`** - Jupyter Notebook示例
- **`DATA/`** - 示例数据文件夹

## ✅ **优势 (Benefits)**

1. **🗂️ 更好的组织** - 计算相关的函数集中管理
2. **📚 清晰的模块化** - 每个模块职责明确
3. **🔍 易于发现** - 用户更容易找到需要的函数
4. **🛠️ 便于维护** - 开发者更容易维护和扩展
5. **📖 完整的文档** - 改进的函数文档字符串

## 🚀 **下一步 (Next Steps)**

1. 测试新的模块结构
2. 在Jupyter Notebook中运行示例
3. 使用DATA文件夹中的真实数据
4. 探索其他Skyborn模块的集成使用
