import re

class RegexConfig():
  def __init__(self):
    self.parameters = r"\$:?([a-zA-Z_]*)\$?"
    self.units = r"([0-9]*?\.?[0-9]+(?:[eE][-+]?[0-9]*?)?)\ *?\[([-+.*/()a-zA-Z0-9]*)\]"
    self.units_b = r"([a-zA-Z]*)"
    self.symbolic = r"\{([-+.*/() 0-9eE]*)\}"
    self.sanitize = r"[a-z-[e]A-Z-[E]]"


regexConfig = RegexConfig()


def symbolicMathRegexHandler(match):
  k = match.group(1)
  if k:
    # Sanitize the input
    sanitized = re.sub(regexConfig.sanitize, '', k).strip()
    value = eval(sanitized, {'__builtins__': None})
    return str(value)
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


parameterHandler = DictRegexHandler()



