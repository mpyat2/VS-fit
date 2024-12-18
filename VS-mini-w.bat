@echo off
SET PYTHON_EXECUTABLE=python
start "" "%COMSPEC%" /c %PYTHON_EXECUTABLE% .\code\main.py ^& pause
