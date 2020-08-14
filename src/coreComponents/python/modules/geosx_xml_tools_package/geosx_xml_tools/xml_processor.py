"""Tools for processing xml files in GEOSX"""

from lxml import etree as ElementTree
from lxml.etree import XMLSyntaxError
import re
import os
from geosx_xml_tools import regex_tools, unit_manager
from geosx_xml_tools import xml_formatter


# Create an instance of the unit, parameter regex handlers
unitManager = unit_manager.UnitManager()
parameterHandler = regex_tools.DictRegexHandler()


def merge_xml_nodes(existingNode, targetNode, level):
  """Merge nodes in an included file into the current structure level by level.

     @param existingNode The current node in the base xml structure.
     @param targetNode The node to insert.
     @param level The xml file depth.
  """

  # Copy attributes on the current level
  for tk in targetNode.attrib.keys():
    existingNode.set(tk, targetNode.get(tk))

  # Copy target children into the xml structure
  currentTag = ''
  matchingSubNodes = []

  for target in targetNode.getchildren():
    insertCurrentLevel = True

    # Check to see if a node with the appropriate type
    # exists at this level
    if (currentTag != target.tag):
      currentTag = target.tag
      matchingSubNodes = existingNode.findall(target.tag)

    if (matchingSubNodes):
      targetName = target.get('name')

      # Special case for the root Problem node (which may be unnamed)
      if (level == 0):
        insertCurrentLevel = False
        merge_xml_nodes(matchingSubNodes[0], target, level + 1)

      # Handle named xml nodes
      elif (targetName and (currentTag not in ['Nodeset'])):
        for match in matchingSubNodes:
          if (match.get('name') == targetName):
            insertCurrentLevel = False
            merge_xml_nodes(match, target, level + 1)

    # Insert any unnamed nodes or named nodes that aren't present
    # in the current xml structure
    if (insertCurrentLevel):
      existingNode.insert(-1, target)


def merge_included_xml_files(root, fname, includeCount, maxInclude=100):
  """Recursively merge included files into the current structure.

     @param root The root node of the base xml structure.
     @param fname The name of the target xml file to merge.
     @param includeCount The current recursion depth.
     @param maxInclude The maximum number of xml files to include (default = 100)
  """

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
      merge_included_xml_files(root, f.get('name'), includeCount)

  # Merge the results into the xml tree
  merge_xml_nodes(root, includeRoot, 0)
  os.chdir(pwd)


def apply_regex_to_node(node):
  """Apply regexes that handle parameters, units, and symbolic math to each
     xml attribute in the structure.

     @param node The target node in the xml structure.
  """

  for k in node.attrib.keys():
    value = node.get(k)

    # Parameter format:  $Parameter or $:Parameter
    ii = 0
    while ('$' in value):
      value = re.sub(regex_tools.patterns['parameters'],
                     parameterHandler,
                     value)
      ii += 1
      if (ii > 100):
        raise Exception('Reached maximum parameter expands (Node=%s, value=%s)' % (node.tag, value))

    # Unit format:       9.81[m**2/s] or 1.0 [bbl/day]
    if ('[' in value):
      value = re.sub(regex_tools.patterns['units'],
                     unitManager.regexHandler,
                     value)

    # Symbolic format:   `1 + 2.34e5*2 * ...`
    ii = 0
    while ('`' in value):
      value = re.sub(regex_tools.patterns['symbolic'],
                     regex_tools.SymbolicMathRegexHandler,
                     value)
      ii += 1
      if (ii > 100):
        raise Exception('Reached maximum symbolic expands (Node=%s, value=%s)' % (node.tag, value))

    node.set(k, value)

  for subNode in node.getchildren():
    apply_regex_to_node(subNode)


def generate_random_name(prefix='', suffix='.xml'):
  """If the target name is not specified, generate a random name for the compiled xml

     @param prefix The file prefix (default = '').
     @param suffix The file suffix (default = '.xml')
  """
  from hashlib import md5
  from time import time
  from os import getpid

  tmp = str(time()) + str(getpid())
  return '%s%s%s' % (prefix, md5(tmp.encode('utf-8')).hexdigest(), suffix)


def process(inputFile, outputFile='', schema='', verbose=0, keep_parameters=True, keep_includes=True):
  """Process an xml file

     @param inputFile Input file name.
     @param outputFile Output file name (if not specified, then generate randomly).
     @param schema Schema file name to validate the final xml (if not specified, then do not validate).
     @param verbose Verbosity level.
     @param keep_parameters If True, then keep parameters in the compiled file (default = True)
     @param keep_includes If True, then keep includes in the compiled file (default = True)
  """
  if verbose:
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
  # Note: doing this first assumes that parameters aren't used in Included block
  includeCount = 0
  for includeNode in root.findall('Included'):
    for f in includeNode.findall('File'):
      merge_included_xml_files(root, f.get('name'), includeCount)
  os.chdir(pwd)

  # Build the parameter map
  Pmap = {}
  for parameters in root.findall('Parameters'):
    for p in parameters.findall('Parameter'):
      Pmap[p.get('name')] = p.get('value')
  parameterHandler.target = Pmap

  # Process any parameters, units, and symbolic math in the xml
  apply_regex_to_node(root)

  # Comment out or remove the Parameter, Included nodes
  for includeNode in root.findall('Included'):
    if keep_includes:
      root.insert(-1, ElementTree.Comment(ElementTree.tostring(includeNode)))
    root.remove(includeNode)
  for parameterNode in root.findall('Parameters'):
    if keep_parameters:
      root.insert(-1, ElementTree.Comment(ElementTree.tostring(parameterNode)))
    root.remove(parameterNode)

  # Generate a random output name if not specified
  if not outputFile:
    outputFile = generate_random_name(prefix='prep_')

  # Write the output file
  tree.write(outputFile, pretty_print=True)

  # Check for un-matched special characters
  with open(outputFile, 'r') as ofile:
    for line in ofile:
      if any([sc in line for sc in ['$', '[', ']', '`']]):
        raise Exception('Found un-matched special characters in the pre-processed input file on line:\n%s\n Check your input xml for errors!' % (line))

  # Apply formatting to the file
  xml_formatter.format_file(outputFile)

  if verbose:
    print('Preprocessed xml file stored in %s' % (outputFile))

  if schema:
    validate_xml(outputFile, schema, verbose)

  return outputFile


def validate_xml(fname, schema, verbose):
  """Validate an xml file, and parse the warnings.

     @param fname Target xml file name.
     @param schema Schema file name.
     @param verbose Verbosity level.
  """
  if verbose:
    print('Validating the xml against the schema...')
  try:
    ofile = ElementTree.parse(fname)
    sfile = ElementTree.XMLSchema(ElementTree.parse(os.path.expanduser(schema)))
    sfile.assertValid(ofile)
  except ElementTree.DocumentInvalid as err:
    print(err)
    print('\nWarning: input XML contains potentially invalid input parameters:')
    print('-'*20+'\n')
    print(sfile.error_log)
    print('\n'+'-'*20)
    print('(Total schema warnings: %i)\n' % (len(sfile.error_log)))

  if verbose:
    print('Done!')


