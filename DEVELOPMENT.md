# Skyborn 开发指南

## 🚀 快速开始开发

### 1. 克隆仓库
```bash
git clone https://github.com/yourusername/skyborn.git
cd skyborn
```

### 2. 自动设置开发环境
```bash
python setup_dev.py
```

这会自动：
- 安装所有开发依赖
- 设置 pre-commit 钩子
- 创建 VS Code 配置

### 3. 手动设置（可选）
```bash
# 安装开发依赖
pip install -e .
pip install pre-commit black isort flake8 mypy pytest

# 设置 pre-commit
pre-commit install
```

## 🔧 开发工具

### Pre-commit 钩子
每次提交时自动运行：
- **Black**: 代码格式化
- **isort**: 导入排序
- **flake8**: 代码质量检查
- **mypy**: 类型检查
- **bandit**: 安全检查

### 手动运行检查
```bash
# 运行所有检查
pre-commit run --all-files

# 单独运行工具
black src/                  # 格式化代码
isort src/                  # 排序导入
flake8 src/                 # 质量检查
mypy src/                   # 类型检查
bandit -r src/              # 安全检查
```

### 测试
```bash
# 运行所有测试
pytest

# 运行测试并生成覆盖率报告
pytest --cov=src/skyborn --cov-report=html

# 运行特定测试
pytest tests/test_conversion.py
```

## 📝 代码规范

### 代码风格
- 使用 **Black** 格式化（行长度 88）
- 使用 **isort** 排序导入
- 遵循 **PEP 8** 规范

### 文档风格
- 使用 **NumPy 风格** 的 docstring
- 所有公共函数必须有文档
- 包含参数、返回值和示例

### 示例 docstring：
```python
def example_function(param1: str, param2: int = 10) -> bool:
    """
    Example function with proper documentation.

    Parameters
    ----------
    param1 : str
        Description of param1.
    param2 : int, default 10
        Description of param2.

    Returns
    -------
    bool
        Description of return value.

    Examples
    --------
    >>> example_function("test", 5)
    True
    """
    return True
```

## 🧪 CI/CD 流程

### GitHub Actions 工作流

1. **test-cross-platform.yml**
   - 测试 Python 3.8-3.12
   - 测试 Windows/macOS/Linux
   - 代码质量检查

2. **pre-commit.yml**
   - 运行 pre-commit 检查
   - 自动修复格式问题

3. **test-coverage.yml**
   - 代码覆盖率分析
   - 性能测试
   - 文档构建测试

4. **publish.yml**
   - 自动发布到 PyPI
   - 创建 GitHub Release

## 📦 发布流程

### 发布新版本
1. 更新版本号：`src/skyborn/__init__.py`
2. 更新 CHANGELOG
3. 创建标签并推送：
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```
4. GitHub Actions 自动构建和发布

## 🔐 Secrets 设置

在 GitHub 仓库设置中添加以下 secrets：

### 必需的 Secrets
- `PYPI_API_TOKEN`: PyPI 发布令牌

### 可选的 Secrets  
- `TEST_PYPI_API_TOKEN`: 测试 PyPI 令牌
- `CODECOV_TOKEN`: Codecov 上传令牌

### 获取 Codecov Token
1. 访问 [codecov.io](https://codecov.io)
2. 用 GitHub 账号登录
3. 添加您的仓库
4. 复制 token 并添加到 GitHub Secrets

## 🆘 常见问题

### Pre-commit 失败
```bash
# 跳过 pre-commit（不推荐）
git commit -m "message" --no-verify

# 修复问题后重新提交
pre-commit run --all-files
git add .
git commit -m "fix: code formatting"
```

### 测试失败
```bash
# 查看详细错误
pytest -v --tb=long

# 运行特定测试
pytest tests/test_conversion.py::TestGribToNetCDF::test_basic -v
```

### 构建失败
```bash
# 检查包构建
python -m build
twine check dist/*
```

## 📚 更多资源

- [Pre-commit 文档](https://pre-commit.com/)
- [Black 文档](https://black.readthedocs.io/)
- [pytest 文档](https://docs.pytest.org/)
- [Codecov 文档](https://docs.codecov.com/)
