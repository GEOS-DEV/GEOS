from lxml import etree as ElementTree
import os
import argparse
import re


def parse_schema_element(
    root,
    node,
    path,
    xsd="{http://www.w3.org/2001/XMLSchema}",
    recursive_types=["PeriodicEvent", "SoloEvent", "HaltEvent"],
):
    """
  @brief Check element/attribute names at the current schema level
  @param root the root schema node
  @param node the current node
  @param path the path to the current level
  @param xsd a string cointaining namespace information
  @param recursive_types a list of elements which are recursive
  """
    # Attributes should be in camelCase, Elements should be in PascalCase
    camel_case_regex = re.compile(r"[a-z][a-zA-Z]*")
    pascal_case_regex = re.compile(r"[A-Z][a-zA-Z]*")

    element_type = node.get("type")
    element_name = node.get("name")
    element_def = root.find("%scomplexType[@name='%s']" % (xsd, element_type))

    # Parse attributes
    for attribute in element_def.findall("%sattribute" % (xsd)):
        attribute_name = attribute.get("name")
        if not camel_case_regex.match(attribute_name):
            print("attribute is not camelCase: %s/%s" % (path, attribute_name))

    # Parse children
    choice_node = element_def.findall("%schoice" % (xsd))
    if choice_node:
        for child in choice_node[0].findall("%selement" % (xsd)):
            child_name = child.get("name")
            if not pascal_case_regex.match(child_name):
                print("Element is not PascalCase: %s/%s" % (path, child_name))

            if not ((child_name in recursive_types) and (element_name in recursive_types)):
                sub_path = "%s/%s" % (path, child_name)
                parse_schema_element(root, child, sub_path)


def parse_schema(fname):
    """
  @brief Check the element/attribute names in the schema
  @param fname the schema file name
  """
    xml_tree = ElementTree.parse(fname)
    xml_root = xml_tree.getroot()
    problem_node = xml_root.find("{http://www.w3.org/2001/XMLSchema}element")
    parse_schema_element(xml_root, problem_node, "Problem")


def main():
    """
  @brief Entry point for the name checking script
  """
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--root", type=str, help="GEOSX root", default=".")
    args = parser.parse_args()

    # Parse the xml files
    geosx_root = os.path.expanduser(args.root)
    schema = "%s/src/coreComponents/fileIO/schema/schema.xsd" % (geosx_root)
    parse_schema(schema)


if __name__ == "__main__":
    main()
