#!/usr/bin/env python3
"""
设置脚本：安装开发环境和 pre-commit 钩子

运行此脚本来设置完整的开发环境，包括 pre-commit 钩子。
"""

import subprocess
import sys
import os
from pathlib import Path

def run_command(cmd, check=True, shell=False):
    """运行命令并打印输出"""
    print(f"🔧 Running: {cmd}")
    try:
        if shell:
            result = subprocess.run(cmd, shell=True, check=check, text=True)
        else:
            result = subprocess.run(cmd.split(), check=check, text=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {e}")
        return False

def check_python_version():
    """检查 Python 版本"""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print("❌ Python 3.8+ is required")
        return False
    print(f"✅ Python {version.major}.{version.minor}.{version.micro}")
    return True

def install_dev_dependencies():
    """安装开发依赖"""
    print("\n📦 Installing development dependencies...")
    
    # 安装包本身（可编辑模式）
    if not run_command("pip install -e ."):
        return False
    
    # 安装开发工具
    dev_packages = [
        "pre-commit",
        "black",
        "isort",
        "flake8",
        "mypy",
        "pytest",
        "pytest-cov",
        "bandit",
        "safety"
    ]
    
    for package in dev_packages:
        if not run_command(f"pip install {package}"):
            print(f"⚠️ Failed to install {package}, continuing...")
    
    return True

def setup_pre_commit():
    """设置 pre-commit 钩子"""
    print("\n🪝 Setting up pre-commit hooks...")
    
    # 安装 pre-commit 钩子
    if not run_command("pre-commit install"):
        return False
    
    # 运行一次以确保工作正常
    print("🧪 Running pre-commit on all files (first time setup)...")
    run_command("pre-commit run --all-files", check=False)
    
    return True

def create_vscode_settings():
    """创建 VS Code 设置"""
    print("\n🔧 Creating VS Code settings...")
    
    vscode_dir = Path(".vscode")
    vscode_dir.mkdir(exist_ok=True)
    
    settings = {
        "python.defaultInterpreterPath": "./venv/bin/python",
        "python.formatting.provider": "black",
        "python.linting.enabled": True,
        "python.linting.flake8Enabled": True,
        "python.linting.mypyEnabled": True,
        "python.sortImports.args": ["--profile", "black"],
        "editor.formatOnSave": True,
        "editor.codeActionsOnSave": {
            "source.organizeImports": True
        },
        "[python]": {
            "editor.rulers": [88],
            "editor.tabSize": 4
        }
    }
    
    import json
    with open(vscode_dir / "settings.json", "w") as f:
        json.dump(settings, f, indent=2)
    
    print("✅ VS Code settings created")

def main():
    """主函数"""
    print("🚀 Setting up Skyborn development environment")
    print("=" * 50)
    
    # 检查 Python 版本
    if not check_python_version():
        sys.exit(1)
    
    # 安装依赖
    if not install_dev_dependencies():
        print("❌ Failed to install dependencies")
        sys.exit(1)
    
    # 设置 pre-commit
    if not setup_pre_commit():
        print("❌ Failed to setup pre-commit")
        sys.exit(1)
    
    # 创建 VS Code 设置
    create_vscode_settings()
    
    print("\n" + "=" * 50)
    print("🎉 Development environment setup complete!")
    print("\n📝 Next steps:")
    print("1. 开始开发代码")
    print("2. pre-commit 会在每次提交时自动运行")
    print("3. 运行 'pre-commit run --all-files' 手动检查所有文件")
    print("4. 运行 'pytest' 执行测试")
    print("\n💡 Useful commands:")
    print("  pre-commit run --all-files  # 运行所有检查")
    print("  black src/                  # 格式化代码")
    print("  isort src/                  # 排序导入")
    print("  flake8 src/                 # 代码质量检查")
    print("  pytest tests/               # 运行测试")

if __name__ == "__main__":
    main()
