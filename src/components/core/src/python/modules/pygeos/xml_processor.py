from xml.etree import ElementTree
import re
import sys
import os
from numpy import *
from . import DictRegexHandler, UnitManager

def MergeIncludedXMLFiles(root, fname, includeCount, maxInclude=100):
  # Check to see if the code has fallen into a loop
  includeCount += 1
  if (includeCount > maxInclude):
    raise Exception('Reached maximum recursive includes...  Is there an include loop?')

  # Load target xml
  includeTree = ElementTree.parse(fname)
  includeRoot = includeTree.getroot()

  # Recursively add the includes:
  for includeNode in includeRoot.findall('Included'):
    for f in includeNode.findall('File'):
      MergeIncludedXMLFiles(root, f.get('name'), includeCount)

  # Merge the results into the xml tree
  for topLevelNode in list(includeRoot):
    rootMatchingNodes = root.findall(topLevelNode.tag)
    if (rootMatchingNodes):
      for secondLevelNode in list(topLevelNode):
        rootMatchingNodes[0].insert(-1, secondLevelNode)
    else:
      root.insert(-1, topLevelNode)


def generateRandomName(suffix=''):
  from hashlib import md5
  from time import time

  return '%s_%s' % (md5(str(time())).hexdigest(), suffix)


def symbolicMathRegexHandler(match):
  k = match.group(1)
  print 'symbolic math: %s' % k
  if k:
    # Sanitize the input
    sanitized = re.sub(r"[a-zA-Z]", '', k).strip()
    value = eval(sanitized, {'__builtins__':None})
    return str(value)
  else:
    return


def PreprocessGEOSXML(inputFile, schema='/g/g17/sherman/GEOS/core/src/schema/gpac.xsd'):
  print('\nReading input xml parameters and parsing symbolic math...')

  # Load the xml files and merge includes.
  tree = ElementTree.parse(inputFile)
  root = tree.getroot()
  includeCount = 0
  for includeNode in root.findall('Included'):
    for f in includeNode.findall('File'):
      MergeIncludedXMLFiles(root, f.get('name'), includeCount)

  # Build the parameter map, convert function
  Pmap = {}
  for parameters in root.findall('Parameters'):
    for p in parameters.findall('Parameter'):
      Pmap[p.get('name')] = p.get('value')
  tmp_fname_a = generateRandomName(suffix='int.xml')
  tree.write(tmp_fname_a)
  parameterHandler = DictRegexHandler()
  parameterHandler.target = Pmap

  # Parse the raw xml file:  
  tmp_fname_b = generateRandomName(suffix='prep.xml')
  unitManager = UnitManager()
  with open(tmp_fname_a, 'r') as ifile, open(tmp_fname_b, 'w') as ofile:
    for line in ifile:
      # Fill in any paramters (format:  $Parameter or $:Parameter)
      if ('$' in line):
        line = re.sub(r"\$:?([a-zA-Z_]*\$?)", parameterHandler, line)

      # Parse any units (format: 9.81[m**2/s] or 1.0 [bbl/day])
      if ('[' in line):
        line = re.sub(r"([0-9]*\.?[0-9]*?[eE]?[-+]?[0-9]*?)\ *?\[([-+.*/()a-zA-Z0-9]*)\]", unitManager.regexHandler, line)

      # Evaluate symbolic math (format: `1 + 2.34e5*2 * ...`)
      if ('`' in line):
        line = re.sub(r"`([-+.*/() 0-9eE]*)`", symbolicMathRegexHandler, line)
      
      ofile.write(line)

  # Validate against the schema
  print('Validating the xml against the schema..')
  err = os.system('xmllint --noout --schema %s %s' % (schema, tmp_fname_b))
  if err:
    print('\nWarning: XML includes invalid elements!\n')

  print('Preprocessed xml file stored in %s' % (tmp_fname_b))
  return tmp_fname_b

if (__name__ == "__main__"):
  ifile = sys.argv[1]
  PreprocessGEOSXML(ifile)