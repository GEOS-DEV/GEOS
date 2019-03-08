from lxml import etree as ElementTree
from lxml.etree import XMLSyntaxError
import re
import os
from . import unitManager, parameterHandler, symbolicMathRegexHandler, regexConfig


def mergeXMLNodes(existingNode, targetNode, level):
  # Copy attributes
  for tk in targetNode.attrib.keys():
    existingNode.set(tk, targetNode.get(tk))

  # Copy target children into the xml structure
  currentTag = ''
  matchingSubNodes = []

  for target in targetNode.getchildren():
    insertCurrentLevel = True

    if (currentTag != target.tag):
      currentTag = target.tag
      matchingSubNodes = existingNode.findall(target.tag)

    if (matchingSubNodes):
      targetName = target.get('name')

      if (level == 0):
        insertCurrentLevel = False
        mergeXMLNodes(matchingSubNodes[0], target, level + 1)

      elif (targetName and (currentTag not in ['Nodeset'])):
        for match in matchingSubNodes:
          if (match.get('name') == targetName):
            insertCurrentLevel = False
            mergeXMLNodes(match, target, level + 1)

    if (insertCurrentLevel):
      existingNode.insert(-1, target)


def mergeIncludedXMLFiles(root, fname, includeCount, maxInclude=100):
  # Expand the input path
  pwd = os.getcwd()
  includePath, fname = os.path.split(os.path.abspath(os.path.expanduser(fname)))
  os.chdir(includePath)

  # Check to see if the code has fallen into a loop
  includeCount += 1
  if (includeCount > maxInclude):
    raise Exception('Reached maximum recursive includes...  Is there an include loop?')

  # Check to make sure the file exists
  if (not os.path.isfile(fname)):
    print('Included file does not exist: %s' % (fname))
    raise Exception('Check included file path!')

  # Load target xml
  try:
    parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
    includeTree = ElementTree.parse(fname, parser)
    includeRoot = includeTree.getroot()
  except XMLSyntaxError as err:
    print('\nCould not load included file: %s' % (fname))
    print(err.msg)
    raise Exception('\nCheck included file!')

  # Recursively add the includes:
  for includeNode in includeRoot.findall('Included'):
    for f in includeNode.findall('File'):
      mergeIncludedXMLFiles(root, f.get('name'), includeCount)

  # Merge the results into the xml tree
  mergeXMLNodes(root, includeRoot, 0)
  os.chdir(pwd)


def applyRegexToNode(node):
  for k in node.attrib.keys():
    value = node.get(k)

    # Parameter format:  $Parameter or $:Parameter
    ii = 0
    while ('$' in value):
      value = re.sub(regexConfig.parameters, parameterHandler, value)
      ii += 1
      if (ii > 100):
        raise Exception('Reached maximum parameter expands (Node=%s, value=%s)' % (node.tag, value))

    # Unit format:       9.81[m**2/s] or 1.0 [bbl/day]
    if ('[' in value):
      value = re.sub(regexConfig.units, unitManager.regexHandler, value)

    # Symbolic format:   {1 + 2.34e5*2 * ...}
    ii = 0
    while ('{' in value):
      value = re.sub(regexConfig.symbolic, symbolicMathRegexHandler, value)
      ii += 1
      if (ii > 100):
        raise Exception('Reached maximum symbolic expands (Node=%s, value=%s)' % (node.tag, value))

    node.set(k, value)

  for subNode in node.getchildren():
    applyRegexToNode(subNode)


def generateRandomName(prefix='', suffix='.xml'):
  from hashlib import md5
  from time import time
  from os import getpid

  return '%s%s%s' % (prefix, md5(str(time())+str(getpid())).hexdigest(), suffix)


def preprocessGEOSXML(inputFile, schema='/g/g17/sherman/GEOS/geosx/src/components/core/src/schema/gpac_new.xsd', verbose=1):

  if (verbose > 0):
    print('\nReading input xml parameters and parsing symbolic math...')

  # Expand the input path
  pwd = os.getcwd()
  rootPath, inputFile = os.path.split(os.path.abspath(os.path.expanduser(inputFile)))
  os.chdir(rootPath)

  # Load the xml files and merge includes
  try:
    parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
    tree = ElementTree.parse(inputFile, parser=parser)
    root = tree.getroot()
  except XMLSyntaxError as err:
    print('\nCould not load input file: %s' % (inputFile))
    print(err.msg)
    raise Exception('\nCheck input file!')

  # Add the included files to the xml structure
  includeCount = 0
  for includeNode in root.findall('Included'):
    for f in includeNode.findall('File'):
      mergeIncludedXMLFiles(root, f.get('name'), includeCount)
  os.chdir(pwd)

  # Build the parameter map
  Pmap = {}
  for parameters in root.findall('Parameters'):
    for p in parameters.findall('Parameter'):
      Pmap[p.get('name')] = p.get('value')
  parameterHandler.target = Pmap

  # Apply regexes and write the new file
  applyRegexToNode(root)
  recordName = generateRandomName(prefix='prep_')
  tree.write(recordName, pretty_print=True)

  # Check for un-matched special characters
  with open(recordName, 'r') as ofile:
    for line in ofile:
      if any([sc in line for sc in ['$', '[', ']', '{', '}']]):
        raise Exception('Found un-matched special characters in the pre-processed input file on line:\n%s\n Check your input xml for errors!' % (line))

  if (verbose > 0):
    print('Preprocessed xml file stored in %s' % (recordName))
    validateXML(recordName, schema)

  return recordName


def validateXML(fname, schema):
  print('Validating the xml against the schema...')
  try:
    ofile = ElementTree.parse(fname)
    sfile = ElementTree.XMLSchema(ElementTree.parse(os.path.expanduser(schema)))
    sfile.assertValid(ofile)
    print('Done!')
  except ElementTree.DocumentInvalid as err:
    print(err)
    print('\nWarning: input XML contains potentially invalid input parameters:')
    print('-'*20+'\n')
    print(sfile.error_log)
    print('\n'+'-'*20)
    print('(Total schema warnings: %i)\n' % (len(sfile.error_log)))


