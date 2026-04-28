@echo off
SET PYTHON_EXECUTABLE=pyw
cd /d "%~dp0"
cd
start %PYTHON_EXECUTABLE% .\code\main.py
