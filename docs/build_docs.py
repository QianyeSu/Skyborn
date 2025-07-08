#!/usr/bin/env python3
"""
Multi-language documentation builder for Skyborn.

This script builds both English and Chinese versions of the documentation
using Sphinx with proper configuration for each language.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

# Configuration
DOCS_DIR = Path(__file__).parent
SOURCE_DIR = DOCS_DIR / "source"
BUILD_DIR = DOCS_DIR / "build"
LANGUAGES = {
    'en': 'English',
    'zh_CN': '中文 (简体)'
}

def run_command(cmd, cwd=None):
    """Run a shell command and return the result."""
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            cwd=cwd,
            capture_output=True, 
            text=True,
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error output: {e.stderr}")
        return None

def check_dependencies():
    """Check if required packages are installed."""
    required_packages = [
        'sphinx',
        'sphinx-rtd-theme', 
        'myst-parser',
        'sphinx-autodoc-typehints'
    ]
    
    missing = []
    for package in required_packages:
        try:
            __import__(package.replace('-', '_'))
        except ImportError:
            missing.append(package)
    
    if missing:
        print("Missing required packages:")
        for pkg in missing:
            print(f"  - {pkg}")
        print("\nInstall with:")
        print(f"pip install {' '.join(missing)}")
        return False
    
    return True

def clean_build():
    """Clean the build directory."""
    if BUILD_DIR.exists():
        print("Cleaning build directory...")
        shutil.rmtree(BUILD_DIR)
    BUILD_DIR.mkdir(exist_ok=True)

def build_language(lang_code):
    """Build documentation for a specific language."""
    lang_name = LANGUAGES[lang_code]
    source_path = SOURCE_DIR / lang_code
    output_path = BUILD_DIR / "html" / lang_code
    
    print(f"\nBuilding {lang_name} documentation...")
    print(f"Source: {source_path}")
    print(f"Output: {output_path}")
    
    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Build command
    cmd = f"sphinx-build -b html {source_path} {output_path}"
    
    # Run sphinx-build
    result = run_command(cmd, cwd=DOCS_DIR)
    
    if result is not None:
        print(f"✅ {lang_name} documentation built successfully!")
        print(f"   Open: {output_path}/index.html")
        return True
    else:
        print(f"❌ Failed to build {lang_name} documentation")
        return False

def create_index_page():
    """Create a main index page with language selection."""
    index_content = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Skyborn Documentation</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 800px;
            margin: 50px auto;
            padding: 20px;
            text-align: center;
        }
        .language-selector {
            display: flex;
            justify-content: center;
            gap: 30px;
            margin-top: 50px;
        }
        .language-option {
            background: #f8f9fa;
            border: 2px solid #e9ecef;
            border-radius: 10px;
            padding: 30px;
            text-decoration: none;
            color: #333;
            transition: all 0.3s ease;
            min-width: 200px;
        }
        .language-option:hover {
            background: #e9ecef;
            border-color: #007bff;
            transform: translateY(-5px);
        }
        .language-title {
            font-size: 24px;
            font-weight: bold;
            margin-bottom: 10px;
        }
        .language-desc {
            color: #666;
        }
    </style>
</head>
<body>
    <h1>Skyborn Documentation</h1>
    <p>Select your preferred language / 选择您的语言</p>
    
    <div class="language-selector">
        <a href="html/en/index.html" class="language-option">
            <div class="language-title">English</div>
            <div class="language-desc">English Documentation</div>
        </a>
        
        <a href="html/zh_CN/index.html" class="language-option">
            <div class="language-title">中文</div>
            <div class="language-desc">中文文档</div>
        </a>
    </div>
    
    <footer style="margin-top: 50px; color: #666;">
        <p>Skyborn - Atmospheric Data Processing Library</p>
    </footer>
</body>
</html>"""
    
    index_file = BUILD_DIR / "index.html"
    with open(index_file, 'w', encoding='utf-8') as f:
        f.write(index_content)
    
    print(f"✅ Main index page created: {index_file}")

def main():
    """Main build function."""
    print("🚀 Building Skyborn Multi-language Documentation")
    print("=" * 50)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Clean build directory
    clean_build()
    
    # Build each language
    success_count = 0
    for lang_code in LANGUAGES:
        if build_language(lang_code):
            success_count += 1
    
    # Create main index page
    create_index_page()
    
    # Summary
    print("\n" + "=" * 50)
    print(f"📊 Build Summary:")
    print(f"   Languages built: {success_count}/{len(LANGUAGES)}")
    
    if success_count == len(LANGUAGES):
        print("🎉 All documentation built successfully!")
        print(f"\n📖 View documentation:")
        print(f"   Main page: file://{BUILD_DIR.absolute()}/index.html")
        print(f"   English:   file://{BUILD_DIR.absolute()}/html/en/index.html")
        print(f"   Chinese:   file://{BUILD_DIR.absolute()}/html/zh_CN/index.html")
    else:
        print("⚠️  Some builds failed. Check the output above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()
