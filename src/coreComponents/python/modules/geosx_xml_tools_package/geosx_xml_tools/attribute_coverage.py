from lxml import etree as ElementTree
import os
from pathlib import Path
import argparse


def parse_schema_element(root,
                         node,
                         xsd='{http://www.w3.org/2001/XMLSchema}',
                         recursive_types=['PeriodicEvent', 'SoloEvent', 'HaltEvent'],
                         folders=['src', 'examples']):
    """Parse the xml schema at the current level

    Args:
        root (lxml.etree.Element): the root schema node
        node (lxml.etree.Element): current schema node
        xsd (str): the file namespace
        recursive_types (list): node tags that allow recursive nesting
        folders (list): folders to sort xml attribute usage into

    Returns:
        dict: Dictionary of attributes and children for the current node
    """

    element_type = node.get('type')
    element_name = node.get('name')
    element_def = root.find("%scomplexType[@name='%s']" % (xsd, element_type))
    local_types = {'attributes': {}, 'children': {}}

    # Parse attributes
    for attribute in element_def.findall('%sattribute' % (xsd)):
        attribute_name = attribute.get('name')
        local_types['attributes'][attribute_name] = {ka: [] for ka in folders}
        if ('default' in attribute.attrib):
            local_types['attributes'][attribute_name]['default'] = attribute.get('default')

    # Parse children
    choice_node = element_def.findall('%schoice' % (xsd))
    if choice_node:
        for child in choice_node[0].findall('%selement' % (xsd)):
            child_name = child.get('name')
            if not ((child_name in recursive_types) and (element_name in recursive_types)):
                local_types['children'][child_name] = parse_schema_element(root, child)

    return local_types


def parse_schema(fname):
    """Parse the schema file into the xml attribute usage dict

    Args:
        fname (str): schema name

    Returns:
        dict: Dictionary of attributes and children for the entire schema
    """
    xml_tree = ElementTree.parse(fname)
    xml_root = xml_tree.getroot()
    problem_node = xml_root.find("{http://www.w3.org/2001/XMLSchema}element")
    return {'Problem': parse_schema_element(xml_root, problem_node)}


def collect_xml_attributes_level(local_types, node, folder):
    """Collect xml attribute usage at the current level

    Args:
        local_types (dict): dictionary containing attribute usage
        node (lxml.etree.Element): current xml node
        folder (str): the source folder for the current file
    """
    for ka in node.attrib.keys():
        local_types['attributes'][ka][folder].append(node.get(ka))

    for child in node:
        if child.tag in local_types['children']:
            collect_xml_attributes_level(local_types['children'][child.tag], child, folder)


def collect_xml_attributes(xml_types, fname, folder):
    """Collect xml attribute usage in a file

    Args:
        xml_types (dict): dictionary containing attribute usage
        fname (str): name of the target file
        folder (str): the source folder for the current file
    """
    parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
    xml_tree = ElementTree.parse(fname, parser=parser)
    xml_root = xml_tree.getroot()

    collect_xml_attributes_level(xml_types['Problem'], xml_root, folder)


def write_attribute_usage_xml_level(local_types, node, folders=['src', 'examples']):
    """Write xml attribute usage file at a given level

    Args:
        local_types (dict): dict containing attribute usage at the current level
        node (lxml.etree.Element): current xml node
    """

    # Write attributes
    for ka in local_types['attributes'].keys():
        attribute_node = ElementTree.Element(ka)
        node.append(attribute_node)

        if ('default' in local_types['attributes'][ka]):
            attribute_node.set('default', local_types['attributes'][ka]['default'])

        unique_values = []
        for f in folders:
            sub_values = list(set(local_types['attributes'][ka][f]))
            unique_values.extend(sub_values)
            attribute_node.set(f, ' | '.join(sub_values))

        unique_length = len(set(unique_values))
        attribute_node.set('unique_values', str(unique_length))

    # Write children
    for ka in sorted(local_types['children']):
        child = ElementTree.Element(ka)
        node.append(child)
        write_attribute_usage_xml_level(local_types['children'][ka], child)


def write_attribute_usage_xml(xml_types, fname):
    """Write xml attribute usage file

    Args:
        xml_types (dict): dictionary containing attribute usage by xml type
        fname (str): output file name
    """
    xml_root = ElementTree.Element('Problem')
    xml_tree = ElementTree.ElementTree(xml_root)

    write_attribute_usage_xml_level(xml_types['Problem'], xml_root)
    xml_tree.write(fname, pretty_print=True)


def process_xml_files(geosx_root, output_name):
    """Test for xml attribute usage

    Args:
        geosx_root (str): GEOSX root directory
        output_name (str): output file name
    """

    # Parse the schema
    geosx_root = os.path.expanduser(geosx_root)
    schema = '%ssrc/coreComponents/schema/schema.xsd' % (geosx_root)
    xml_types = parse_schema(schema)

    # Find all xml files, collect their attributes
    for folder in ['src', 'examples']:
        print(folder)
        xml_files = Path(os.path.join(geosx_root, folder)).rglob('*.xml')
        for f in xml_files:
            print('  %s' % (str(f)))
            collect_xml_attributes(xml_types, str(f), folder)

    # Consolidate attributes
    write_attribute_usage_xml(xml_types, output_name)


def main():
    """Entry point for the xml attribute usage test script

    Args:
        -r/--root (str): GEOSX root directory
        -o/--output (str): output file name
    """

    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--root', type=str, help='GEOSX root', default='')
    parser.add_argument('-o', '--output', type=str, help='Output file name', default='attribute_test.xml')
    args = parser.parse_args()

    # Parse the xml files
    process_xml_files(args.root, args.output)


if __name__ == "__main__":
    main()
