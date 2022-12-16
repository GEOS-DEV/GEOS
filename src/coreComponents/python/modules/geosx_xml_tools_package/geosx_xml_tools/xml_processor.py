"""Tools for processing xml files in GEOSX"""

from lxml import etree as ElementTree    # type: ignore[import]
from lxml.etree import XMLSyntaxError    # type: ignore[import]
import re
import os
from geosx_xml_tools import regex_tools, unit_manager
from geosx_xml_tools import xml_formatter
from typing import Iterable, Tuple, List

# Create an instance of the unit, parameter regex handlers
unitManager = unit_manager.UnitManager()
parameterHandler = regex_tools.DictRegexHandler()


def merge_xml_nodes(existingNode: ElementTree.Element, targetNode: ElementTree.Element, level: int) -> None:
    """Merge nodes in an included file into the current structure level by level.

    Args:
        existingNode (lxml.etree.Element): The current node in the base xml structure.
        targetNode (lxml.etree.Element): The node to insert.
        level (int): The xml file depth.
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


def merge_included_xml_files(root: ElementTree.Element, fname: str, includeCount: int, maxInclude: int = 100) -> None:
    """Recursively merge included files into the current structure.

    Args:
        root (lxml.etree.Element): The root node of the base xml structure.
        fname (str): The name of the target xml file to merge.
        includeCount (int): The current recursion depth.
        maxInclude (int): The maximum number of xml files to include (default = 100)
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


def apply_regex_to_node(node: ElementTree.Element) -> None:
    """Apply regexes that handle parameters, units, and symbolic math to each
    xml attribute in the structure.

    Args:
        node (lxml.etree.Element): The target node in the xml structure.
    """

    for k in node.attrib.keys():
        value = node.get(k)

        # Parameter format:  $Parameter or $:Parameter
        ii = 0
        while ('$' in value):
            value = re.sub(regex_tools.patterns['parameters'], parameterHandler, value)
            ii += 1
            if (ii > 100):
                raise Exception('Reached maximum parameter expands (Node=%s, value=%s)' % (node.tag, value))

        # Unit format:       9.81[m**2/s] or 1.0 [bbl/day]
        if ('[' in value):
            value = re.sub(regex_tools.patterns['units'], unitManager.regexHandler, value)

        # Symbolic format:   `1 + 2.34e5*2 * ...`
        ii = 0
        while ('`' in value):
            value = re.sub(regex_tools.patterns['symbolic'], regex_tools.SymbolicMathRegexHandler, value)
            ii += 1
            if (ii > 100):
                raise Exception('Reached maximum symbolic expands (Node=%s, value=%s)' % (node.tag, value))

        node.set(k, value)

    for subNode in node.getchildren():
        apply_regex_to_node(subNode)


def generate_random_name(prefix: str = '', suffix: str = '.xml') -> str:
    """If the target name is not specified, generate a random name for the compiled xml

    Args:
        prefix (str): The file prefix (default = '').
        suffix (str): The file suffix (default = '.xml')

    Returns:
        str: Random file name
    """
    from hashlib import md5
    from time import time
    from os import getpid

    tmp = str(time()) + str(getpid())
    return '%s%s%s' % (prefix, md5(tmp.encode('utf-8')).hexdigest(), suffix)


def process(inputFiles: Iterable[str],
            outputFile: str = '',
            schema: str = '',
            verbose: int = 0,
            parameter_override: List[Tuple[str, str]] = [],
            keep_parameters: bool = True,
            keep_includes: bool = True) -> str:
    """Process an xml file

    Args:
        inputFiles (list): Input file names.
        outputFile (str): Output file name (if not specified, then generate randomly).
        schema (str): Schema file name to validate the final xml (if not specified, then do not validate).
        verbose (int): Verbosity level.
        parameter_override (list): Parameter value overrides
        keep_parameters (bool): If True, then keep parameters in the compiled file (default = True)
        keep_includes (bool): If True, then keep includes in the compiled file (default = True)

    Returns:
        str: Output file name
    """
    if verbose:
        print('\nReading input xml parameters and parsing symbolic math...')

    # Check the type of inputFiles
    if isinstance(inputFiles, str):
        inputFiles = [inputFiles]

    # Expand the input path
    pwd = os.getcwd()
    expanded_files = [os.path.abspath(os.path.expanduser(f)) for f in inputFiles]
    single_path, single_input = os.path.split(expanded_files[0])
    os.chdir(single_path)

    # Handle single vs. multiple command line inputs
    root = ElementTree.Element("Problem")
    tree = ElementTree.ElementTree()
    if (len(expanded_files) == 1):
        # Load single files directly
        try:
            parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
            tree = ElementTree.parse(single_input, parser=parser)
            root = tree.getroot()
        except XMLSyntaxError as err:
            print('\nCould not load input file: %s' % (single_input))
            print(err.msg)
            raise Exception('\nCheck input file!')

    else:
        # For multiple inputs, create a simple xml structure to hold
        # the included files.  These will be saved as comments in the compiled file
        root = ElementTree.Element('Problem')
        tree = ElementTree.ElementTree(root)
        included_node = ElementTree.Element("Included")
        root.append(included_node)
        for f in expanded_files:
            included_file = ElementTree.Element("File")
            included_file.set('name', f)
            included_node.append(included_file)

    # Add the included files to the xml structure
    # Note: doing this first assumes that parameters aren't used in Included block
    includeCount = 0
    for includeNode in root.findall('Included'):
        for f in includeNode.findall('File'):
            merge_included_xml_files(root, f.get('name'), includeCount)    # type: ignore[attr-defined]
    os.chdir(pwd)

    # Build the parameter map
    Pmap = {}
    for parameters in root.findall('Parameters'):
        for p in parameters.findall('Parameter'):
            Pmap[p.get('name')] = p.get('value')

    # Apply any parameter overrides
    if len(parameter_override):
        # Save overriden values to a new xml element
        command_override_node = ElementTree.Element("CommandLineOverride")
        root.append(command_override_node)
        for ii in range(len(parameter_override)):
            pname = parameter_override[ii][0]
            pval = ' '.join(parameter_override[ii][1:])
            Pmap[pname] = pval
            override_parameter = ElementTree.Element("Parameter")
            override_parameter.set('name', pname)
            override_parameter.set('value', pval)
            command_override_node.append(override_parameter)

    # Add the parameter map to the handler
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
    for overrideNode in root.findall('CommandLineOverride'):
        if keep_parameters:
            root.insert(-1, ElementTree.Comment(ElementTree.tostring(overrideNode)))
        root.remove(overrideNode)

    # Generate a random output name if not specified
    if not outputFile:
        outputFile = generate_random_name(prefix='prep_')

    # Write the output file
    tree.write(outputFile, pretty_print=True)

    # Check for un-matched special characters
    with open(outputFile, 'r') as ofile:
        for line in ofile:
            if any([sc in line for sc in ['$', '[', ']', '`']]):
                raise Exception(
                    'Found un-matched special characters in the pre-processed input file on line:\n%s\n Check your input xml for errors!'
                    % (line))

    # Apply formatting to the file
    xml_formatter.format_file(outputFile)

    if verbose:
        print('Preprocessed xml file stored in %s' % (outputFile))

    if schema:
        validate_xml(outputFile, schema, verbose)

    return outputFile


def validate_xml(fname: str, schema: str, verbose: int) -> None:
    """Validate an xml file, and parse the warnings.

    Args:
        fname (str): Target xml file name.
        schema (str): Schema file name.
        verbose (int): Verbosity level.
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
        print('-' * 20 + '\n')
        print(sfile.error_log)
        print('\n' + '-' * 20)
        print('(Total schema warnings: %i)\n' % (len(sfile.error_log)))

    if verbose:
        print('Done!')
