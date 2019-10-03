import re


# Config
patterns = {'parameters': r"\$:?([a-zA-Z_]*)\$?",
            'units': r"([0-9]*?\.?[0-9]+(?:[eE][-+]?[0-9]*?)?)\ *?\[([-+.*/()a-zA-Z0-9]*)\]",
            'units_b': r"([a-zA-Z]*)",
            'symbolic': r"\`([-+.*/() 0-9eE]*)\`",
            'sanitize': r"[a-z-[e]A-Z-[E]]",
            'strip_trailing': r"\.?0+(?=e)",
            'strip_trailing_b': r"e\+00|\+0?|(?<=-)0"}
symbolic_format = '%1.6e'


def SymbolicMathRegexHandler(match):
  k = match.group(1)
  if k:
    # Sanitize the input
    sanitized = re.sub(patterns['sanitize'], '', k).strip()
    value = eval(sanitized, {'__builtins__': None})

    # Format the string, removing any trailing zeros and decimals
    str_value = re.sub(patterns['strip_trailing'], '', symbolic_format % (value))
    # Strip + and trailing zero exponents
    str_value = re.sub(patterns['strip_trailing_b'], '', str_value)
    return str_value
  else:
    return


class DictRegexHandler():
  def __init__(self):
    self.target = {}

  def __call__(self, match):
    k = match.group(1)
    if k:
      if (k not in self.target.keys()):
          raise Exception('Error: Target (%s) is not defined in the regex handler' % k)
      value = self.target[k]
      return str(value)
    else:
      return




