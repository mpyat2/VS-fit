import os
import subprocess
import sys

def open_pdf(path):
    if sys.platform.startswith("win"):
        os.startfile(path)
    elif sys.platform.startswith("darwin"):
        subprocess.run(["open", path])
    else:
        subprocess.run(["xdg-open", path])

def open_folder(path):
    if sys.platform.startswith("win"):
        os.startfile(path)
    elif sys.platform.startswith("darwin"):
        subprocess.run(["open", path])
    else:
        subprocess.run(["xdg-open", path])

def safe_eval(s):
  return eval(s, {"__builtins__": {}}, {})
