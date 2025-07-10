@echo off
echo Building Skyborn documentation...
python "%~dp0build_docs.py" %*
echo.
echo Documentation build process finished.
pause
