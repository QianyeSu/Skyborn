#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-language documentation builder for Skyborn.

This script builds both English and Chinese versions of the documentation
using Sphinx with proper configuration for each language.
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path

# Set proper encoding for Windows
os.environ['PYTHONIOENCODING'] = 'utf-8'

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
        # Set proper encoding for Windows
        result = subprocess.run(
            cmd, 
            shell=True, 
            cwd=cwd,
            capture_output=True, 
            text=True,
            encoding='utf-8',
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        if e.stderr:
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
    output_path = BUILD_DIR / lang_code / "html"
    
    print(f"\nBuilding {lang_name} documentation...")
    print(f"Source: {source_path}")
    print(f"Output: {output_path}")
    
    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Build command with explicit encoding
    cmd = f"sphinx-build -b html -E -a {source_path} {output_path}"
    
    # Set environment variables for proper encoding
    env = os.environ.copy()
    env['PYTHONIOENCODING'] = 'utf-8'
    env['LANG'] = 'en_US.UTF-8'
    env['LC_ALL'] = 'en_US.UTF-8'
    
    try:
        # Run sphinx-build with proper encoding
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=DOCS_DIR,
            capture_output=True,
            text=True,
            encoding='utf-8',
            env=env,
            check=True
        )
        print(f"✅ {lang_name} documentation built successfully!")
        print(f"   Open: {output_path}/index.html")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to build {lang_name} documentation")
        print(f"Command: {cmd}")
        if e.stdout:
            print(f"Output: {e.stdout}")
        if e.stderr:
            print(f"Error: {e.stderr}")
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
        <a href="en/html/index.html" class="language-option">
            <div class="language-title">English</div>
            <div class="language-desc">English Documentation</div>
        </a>
        
        <a href="zh_CN/html/index.html" class="language-option">
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
    parser = argparse.ArgumentParser(description='Build Skyborn documentation')
    parser.add_argument('--lang', choices=['en', 'zh_CN'], 
                       help='Build specific language only')
    parser.add_argument('--clean', action='store_true',
                       help='Clean build directory before building')
    
    args = parser.parse_args()
    
    print("🚀 Building Skyborn Multi-language Documentation")
    print("=" * 50)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Clean build directory if requested
    if args.clean:
        clean_build()
    
    # Build specific language or all languages
    if args.lang:
        languages_to_build = [args.lang]
    else:
        languages_to_build = list(LANGUAGES.keys())
    
    # Build each language
    success_count = 0
    for lang_code in languages_to_build:
        if build_language(lang_code):
            success_count += 1
    
    # Create main index page if building all languages
    if not args.lang:
        create_index_page()
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 Build Summary:")
    print(f"   Languages built: {success_count}/{len(languages_to_build)}")
    
    if success_count == len(languages_to_build):
        print("🎉 All documentation built successfully!")
        if not args.lang:
            print("\n📖 View documentation:")
            print(f"   Main page: file://{BUILD_DIR.absolute()}/index.html")
        print(f"   English:   file://{BUILD_DIR.absolute()}/en/html/index.html")
        print(f"   Chinese:   file://{BUILD_DIR.absolute()}/zh_CN/html/index.html")
    else:
        print("⚠️  Some builds failed. Check the output above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()
