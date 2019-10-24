import re


# Define regex patterns used throughout the module:
#
# Pattern         |  Example targets             | Notes
# ------------------------------------------------------------------------------------
# parameters      | $Parameter, $Parameter$      | Matches entire parameter string
# units           | 9.81[m**2/s],  1.0 [bbl/day] | Matches entire unit string
# units_b         | m, bbl, day                  | Matches unit names
# symbolic        | `1 + 2.34e5*2`               | Matches the entire symbolic string
# sanitize        |                              | Removes any residual characters before
#                 |                              |       evaluating symbolic expressions
# strip_trailing  | 3.0000, 5.150050             | Removes unnecessary float strings
# strip_trailing_b| 3.0000e0, 1.23e0             | Removes unnecessary float strings
#
patterns = {'parameters': r"\$:?([a-zA-Z_]*)\$?",
            'units': r"([0-9]*?\.?[0-9]+(?:[eE][-+]?[0-9]*?)?)\ *?\[([-+.*/()a-zA-Z0-9]*)\]",
            'units_b': r"([a-zA-Z]*)",
            'symbolic': r"\`([-+.*/() 0-9eE]*)\`",
            'sanitize': r"[a-z-[e]A-Z-[E]]",
            'strip_trailing': r"\.?0+(?=e)",
            'strip_trailing_b': r"e\+00|\+0?|(?<=-)0"}

# String formatting for symbolic expressions
symbolic_format = '%1.6e'


# This function is used to evaluate symbolic expressions that are identified
# using the patterns['symbolic'] regex
def SymbolicMathRegexHandler(match):
  k = match.group(1)
  if k:
    # Sanitize the input
    sanitized = re.sub(patterns['sanitize'], '', k).strip()
    value = eval(sanitized, {'__builtins__': None})

    # Format the string, removing any trailing zeros, decimals, etc.
    str_value = re.sub(patterns['strip_trailing'], '', symbolic_format % (value))
    str_value = re.sub(patterns['strip_trailing_b'], '', str_value)
    return str_value
  else:
    return


# This class is used to substitute matched values with those stored in a dict
class DictRegexHandler():
  def __init__(self):
    # Substitution definitions
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


