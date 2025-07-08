# 文档部署故障排除指南

## 🔧 常见的 GitHub Pages 部署问题

### 1. GitHub Pages 设置检查

首先确保您的 GitHub 仓库已正确配置 Pages：

1. 进入 GitHub 仓库设置
2. 点击 "Pages" 选项
3. 在 "Source" 下选择 "GitHub Actions"（不是传统的分支方式）

### 2. 权限问题

确保工作流程有正确的权限：

```yaml
permissions:
  contents: read
  pages: write
  id-token: write
```

### 3. 工作流程触发器

检查工作流程是否被正确触发：

```yaml
on:
  push:
    branches: [ main ]  # 确保分支名称正确
    paths:
      - 'docs/**'
      - 'src/**'
  workflow_dispatch:    # 允许手动触发
```

### 4. 依赖安装问题

最常见的问题是依赖安装失败：

```yaml
- name: Install dependencies
  run: |
    python -m pip install --upgrade pip
    pip install -r requirements.txt
    pip install -r requirements-docs.txt
    pip install -e .
```

### 5. 文档构建失败

检查 Sphinx 构建是否成功：

```yaml
- name: Build documentation
  run: |
    cd docs
    python build_docs.py --lang en --clean
    python build_docs.py --lang zh_CN --clean
```

## 🚀 修复步骤

### 步骤 1: 使用新的工作流程

我已经创建了 `docs-deploy-new.yml`，它包含：

- ✅ 最新的 Actions 版本
- ✅ 改进的错误检查
- ✅ 调试输出
- ✅ 正确的权限设置

### 步骤 2: 启用新工作流程

1. 删除或重命名旧的 `docs-deploy.yml`
2. 将 `docs-deploy-new.yml` 重命名为 `docs-deploy.yml`
3. 推送到 main 分支

### 步骤 3: 手动触发测试

1. 在 GitHub 仓库页面
2. 点击 "Actions" 标签
3. 选择 "Deploy Documentation"
4. 点击 "Run workflow" 手动触发

### 步骤 4: 检查日志

如果失败，检查以下几个关键步骤的日志：

1. **Install dependencies** - 依赖安装
2. **Build English documentation** - 英语文档构建
3. **Build Chinese documentation** - 中文文档构建
4. **Prepare GitHub Pages structure** - 页面结构准备
5. **Deploy to GitHub Pages** - 实际部署

## 🔍 调试技巧

### 本地测试

在推送之前，在本地测试文档构建：

```bash
cd docs
python build_docs.py --clean
```

检查输出：
- `docs/build/en/html/index.html` 应该存在
- `docs/build/zh_CN/html/index.html` 应该存在

### 检查 requirements-docs.txt

确保包含所有必要的包：

```txt
sphinx>=5.0.0
sphinx-rtd-theme>=1.0.0
myst-parser>=0.18.0
sphinx-autodoc-typehints>=1.0.0
nbsphinx>=0.8.0
```

### 检查 Python 版本

确保 Actions 使用的 Python 版本与本地一致：

```yaml
- name: Set up Python
  uses: actions/setup-python@v5
  with:
    python-version: '3.10'
```

## 📋 快速修复清单

- [ ] GitHub Pages 设置为 "GitHub Actions"
- [ ] 工作流程文件语法正确
- [ ] 权限设置正确
- [ ] 依赖文件存在且完整
- [ ] 本地文档构建成功
- [ ] 分支名称匹配（main vs master）
- [ ] 文档构建脚本可执行

## 🆘 如果还是不工作

1. **检查 Actions 日志**：找到具体的错误消息
2. **简化工作流程**：先尝试只构建一种语言
3. **使用官方模板**：参考 GitHub Pages 官方文档示例
4. **检查网络问题**：有时是临时的 GitHub 服务问题

## 🎯 推荐的新工作流程

新的 `docs-deploy-new.yml` 包含：

- 最新的 Actions 版本
- 更好的错误处理
- 详细的调试输出
- 环境变量设置
- 正确的 GitHub Pages 配置

请按照上述步骤操作，应该能解决您的文档部署问题。
