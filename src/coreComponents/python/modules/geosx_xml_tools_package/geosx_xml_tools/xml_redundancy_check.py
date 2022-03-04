
from geosx_xml_tools.attribute_coverage import parse_schema
from geosx_xml_tools.xml_formatter import format_file
from lxml import etree as ElementTree
import os
from pathlib import Path
import argparse


def check_redundancy_level(local_schema, node, whitelist=['component']):
    """Check xml redundancy at the current level

    @arg local_schema dict containing schema definitions
    @arg node current xml node
    """
    node_is_required = 0
    for ka in node.attrib.keys():
        if (ka in whitelist):
            node_is_required += 1
        elif (ka not in local_schema['attributes']):
            node_is_required += 1
        elif ('default' not in local_schema['attributes'][ka]):
            node_is_required += 1
        elif (node.get(ka) != local_schema['attributes'][ka]['default']):
            node_is_required += 1
        else:
            node.attrib.pop(ka)

    for child in node:
        # Comments will not appear in the schema
        if child.tag in local_schema['children']:
            child_is_required = check_redundancy_level(local_schema['children'][child.tag],
                                                       child)
            node_is_required += child_is_required
            if not child_is_required:
                node.remove(child)

    return node_is_required


def check_xml_redundancy(schema, fname):
    """Check redundancy in an xml file

    @arg schema dict containing schema definitions
    @arg fname name of the target file
    """
    xml_tree = ElementTree.parse(fname)
    xml_root = xml_tree.getroot()
    check_redundancy_level(schema['Problem'], xml_root)
    xml_tree.write(fname)
    format_file(fname)


def process_xml_files(geosx_root):
    """Test for xml redundancy

    @arg geosx_root GEOSX root directory
    """

    # Parse the schema
    geosx_root = os.path.expanduser(geosx_root)
    schema_fname = '%ssrc/coreComponents/schema/schema.xsd' % (geosx_root)
    schema = parse_schema(schema_fname)

    # Find all xml files, collect their attributes
    for folder in ['src', 'examples']:
        print(folder)
        xml_files = Path(os.path.join(geosx_root, folder)).rglob('*.xml')
        for f in xml_files:
            print('  %s' % (str(f)))
            check_xml_redundancy(schema, str(f))


def main():
    """Entry point for the xml attribute usage test script

    @arg -r/--root GEOSX root directory
    """

    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--root', type=str, help='GEOSX root', default='')
    args = parser.parse_args()

    # Parse the xml files
    process_xml_files(args.root)


if __name__ == "__main__":
    main()
