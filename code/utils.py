def safe_eval(s):
  return eval(s, {"__builtins__": {}}, {})
