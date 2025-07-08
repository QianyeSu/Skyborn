@echo off
REM Build script for Windows
echo Building Skyborn documentation...

echo.
echo Building English documentation...
sphinx-build -b html source/en docs/build/html/en

echo.
echo Building Chinese documentation...
sphinx-build -b html source/zh_CN docs/build/html/zh_CN

echo.
echo Documentation built successfully!
echo English docs: docs/build/html/en/index.html
echo Chinese docs: docs/build/html/zh_CN/index.html

pause
