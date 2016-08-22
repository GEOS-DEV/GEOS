from lxml import etree as ElementTree
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
  parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
  includeTree = ElementTree.parse(fname, parser)
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


def generateRandomName(prefix='', suffix='.xml'):
  from hashlib import md5
  from time import time

  return '%s%s%s' % (prefix, md5(str(time())).hexdigest(), suffix)


def symbolicMathRegexHandler(match):
  k = match.group(1)
  if k:
    # Sanitize the input
    sanitized = re.sub(r"[a-z-[e]A-Z-[E]]", '', k).strip()
    value = eval(sanitized, {'__builtins__':None})
    return str(value)
  else:
    return


def PreprocessGEOSXML(inputFile, schema='/g/g17/sherman/GEOS/geosx/src/components/core/src/schema/gpac_new.xsd'):
  print('\nReading input xml parameters and parsing symbolic math...')

  # Load the xml files and merge includes
  parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
  tree = ElementTree.parse(inputFile, parser=parser)
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
  tmp_fname_a = generateRandomName(prefix='int_')
  tree.write(tmp_fname_a, pretty_print=True)
  parameterHandler = DictRegexHandler()
  parameterHandler.target = Pmap

  # Parse the raw xml file:  
  tmp_fname_b = generateRandomName(prefix='prep_')
  unitManager = UnitManager()
  with open(tmp_fname_a, 'r') as ifile, open(tmp_fname_b, 'w') as ofile:
    for line in ifile:
      # Fill in any paramters (format:  $Parameter or $:Parameter)
      if ('$' in line):
        line = re.sub(r"\$:?([a-zA-Z_]*\$?)", parameterHandler, line)

      # Parse any units (format: 9.81[m**2/s] or 1.0 [bbl/day])
      if ('[' in line):
        line = re.sub(r"([0-9]*\.?[0-9]*?[eE]?[-+]?[0-9]*?)\ *?\[([-+.*/()a-zA-Z0-9]*)\]", unitManager.regexHandler, line)

      # Evaluate symbolic math (format: {1 + 2.34e5*2 * ...})
      if ('{' in line):
        line = re.sub(r"\{([-+.*/() 0-9eE]*)\}", symbolicMathRegexHandler, line)
      
      ofile.write(line)

  # Check for un-matched special characters
  with open(tmp_fname_b, 'r') as ofile:
    for line in ofile:
      if any([sc in line for sc in ['$', '[', ']', '{', '}']]):
        raise Exception('Found un-matched special characters in the pre-processed input file on line:\n%s\n Check your input xml for errors!' % (line))

  # Validate against the schema 
  print('Validating the xml against the schema...')
  try:
    ofile = ElementTree.parse(tmp_fname_b) 
    sfile = ElementTree.XMLSchema(ElementTree.parse(schema))
    sfile.assertValid(ofile)
  except ElementTree.DocumentInvalid as err:
    print('\nWarning: input XML contains potentially invalid input parameters:')
    print('-'*20+'\n')
    print sfile.error_log
    print('\n'+'-'*20)
    print('(Total schema warnings: %i)\n' % (len(sfile.error_log)))

  os.remove(tmp_fname_a)
  print('Preprocessed xml file stored in %s' % (tmp_fname_b))
  return tmp_fname_b

if (__name__ == "__main__"):
  ifile = sys.argv[1]
  PreprocessGEOSXML(ifile)