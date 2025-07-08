# Pre-commit Configuration Simplification

## 🔧 Changes Made

### Before (复杂配置)
- 多个工具：black, isort, flake8, mypy, bandit
- 复杂的 GitHub Actions workflow (85 行)
- 自动修复机制
- 缓存和错误处理

### After (简化配置)
- 基本工具：pre-commit-hooks, black
- 简单的 GitHub Actions workflow (16 行)
- 直接使用官方 pre-commit actions

## ✅ New Configuration

### `.pre-commit-config.yaml`
```yaml
# See https://pre-commit.com for more information
# Set the default stage to commit, only checking during submission
default_stages:
  - pre-commit

# Set the default language version
default_language_version:
  python: python3.10

repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v5.0.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-symlinks
      - id: check-merge-conflict
      - id: check-json
      - id: check-toml
      - id: name-tests-test

  - repo: https://github.com/psf/black
    rev: 24.10.0
    hooks:
      - id: black
        types: [python]
        exclude: ^docs/
```

### `.github/workflows/pre-commit-ci.yml`
```yaml
name: Workflow for pre-commit

on:
  pull_request:
  push:
    branches: [main, develop]

jobs:
  main:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - uses: pre-commit/action@v3.0.1
      - uses: pre-commit-ci/lite-action@v1.1.0
        if: always()
```

## 🎯 Benefits of Simplification

1. **更少的依赖**: 只有最基本的检查
2. **更快的运行**: 减少了工具链复杂度
3. **更容易维护**: 配置文件简洁明了
4. **官方支持**: 使用 pre-commit 官方 actions
5. **自动修复**: pre-commit-ci/lite-action 提供基本自动修复

## 🚀 Usage

### Local Development
```bash
# 安装 pre-commit
pip install pre-commit

# 安装钩子
pre-commit install

# 手动运行检查
pre-commit run --all-files
```

### What Gets Checked
- ✅ 尾随空格清理
- ✅ 文件末尾换行
- ✅ YAML/JSON/TOML 语法检查
- ✅ 合并冲突检查
- ✅ 测试文件命名规范
- ✅ Python 代码格式化 (black)

## 📋 Removed Features
- ❌ isort (import sorting) - 可手动运行
- ❌ flake8 (linting) - 减少噪音
- ❌ mypy (type checking) - 可单独配置
- ❌ bandit (security) - 有专门的 security workflow
- ❌ 复杂的缓存机制
- ❌ 自定义自动修复逻辑

这样的配置更加简洁和实用，减少了不必要的复杂性，同时保持了核心的代码质量检查。
