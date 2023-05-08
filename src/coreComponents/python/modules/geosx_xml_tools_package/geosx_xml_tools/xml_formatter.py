import os
from lxml import etree as ElementTree    # type: ignore[import]
import re
from typing import List, Any, TextIO
from geosx_xml_tools import command_line_parsers


def format_attribute(attribute_indent: str, ka: str, attribute_value: str) -> str:
    """Format xml attribute strings

    Args:
        attribute_indent (str): Attribute indent string
        ka (str): Attribute name
        attribute_value (str): Attribute value

    Returns:
        str: Formatted attribute value
    """
    # Make sure that a space follows commas
    attribute_value = re.sub(r",\s*", ", ", attribute_value)

    # Handle external brackets
    attribute_value = re.sub(r"{\s*", "{ ", attribute_value)
    attribute_value = re.sub(r"\s*}", " }", attribute_value)

    # Consolidate whitespace
    attribute_value = re.sub(r"\s+", " ", attribute_value)

    # Identify and split multi-line attributes
    if re.match(r"\s*{\s*({[-+.,0-9a-zA-Z\s]*},?\s*)*\s*}", attribute_value):
        split_positions: List[Any] = [match.end() for match in re.finditer(r"}\s*,", attribute_value)]
        newline_indent = '\n%s' % (' ' * (len(attribute_indent) + len(ka) + 4))
        new_values = []
        for a, b in zip([None] + split_positions, split_positions + [None]):
            new_values.append(attribute_value[a:b].strip())
        if new_values:
            attribute_value = newline_indent.join(new_values)

    return attribute_value


def format_xml_level(output: TextIO,
                     node: ElementTree.Element,
                     level: int,
                     indent: str = ' ' * 2,
                     block_separation_max_depth: int = 2,
                     modify_attribute_indent: bool = False,
                     sort_attributes: bool = False,
                     close_tag_newline: bool = False,
                     include_namespace: bool = False) -> None:
    """Iteratively format the xml file

    Args:
        output (file): the output text file handle
        node (lxml.etree.Element): the current xml element
        level (int): the xml depth
        indent (str): the xml indent style
        block_separation_max_depth (int): the maximum depth to separate adjacent elements
        modify_attribute_indent (bool): option to have flexible attribute indentation
        sort_attributes (bool): option to sort attributes alphabetically
        close_tag_newline (bool): option to place close tag on a separate line
        include_namespace (bool): option to include the xml namespace in the output
    """

    # Handle comments
    if node.tag is ElementTree.Comment:
        output.write('\n%s<!--%s-->' % (indent * level, node.text))

    else:
        # Write opening line
        opening_line = '\n%s<%s' % (indent * level, node.tag)
        output.write(opening_line)

        # Write attributes
        if (len(node.attrib) > 0):
            # Choose indentation
            attribute_indent = '%s' % (indent * (level + 1))
            if modify_attribute_indent:
                attribute_indent = ' ' * (len(opening_line))

            # Get a copy of the attributes
            attribute_dict = {}
            if ((level == 0) & include_namespace):
                # Handle the optional namespace information at the root level
                # Note: preferably, this would point to a schema we host online
                attribute_dict['xmlns:xsi'] = 'http://www.w3.org/2001/XMLSchema-instance'
                attribute_dict['xsi:noNamespaceSchemaLocation'] = '/usr/gapps/GEOS/schema/schema.xsd'
            elif (level > 0):
                attribute_dict = node.attrib

            # Sort attribute names
            akeys = list(attribute_dict.keys())
            if sort_attributes:
                akeys = sorted(akeys)

            # Format attributes
            for ka in akeys:
                # Avoid formatting object paths (adds spaces in object sublists which break parsing)
                if node.tag == "FieldSpecification" and ka == "objectPath":
                    continue
                # Avoid formatting mathpresso expressions
                if node.tag in ["SymbolicFunction", "CompositeFunction"] and ka == "expression":
                    continue
                attribute_dict[ka] = format_attribute(attribute_indent, ka, attribute_dict[ka])

            for ii in range(0, len(akeys)):
                k = akeys[ii]
                if ((ii == 0) & modify_attribute_indent):
                    output.write(' %s=\"%s\"' % (k, attribute_dict[k]))
                else:
                    output.write('\n%s%s=\"%s\"' % (attribute_indent, k, attribute_dict[k]))

        # Write children
        if len(node):
            output.write('>')
            Nc = len(node)
            for ii, child in zip(range(Nc), node):
                format_xml_level(output, child, level + 1, indent, block_separation_max_depth, modify_attribute_indent,
                                 sort_attributes, close_tag_newline, include_namespace)

                # Add space between blocks
                if ((level < block_separation_max_depth) & (ii < Nc - 1) & (child.tag is not ElementTree.Comment)):
                    output.write('\n')

            # Write the end tag
            output.write('\n%s</%s>' % (indent * level, node.tag))
        else:
            if close_tag_newline:
                output.write('\n%s/>' % (indent * level))
            else:
                output.write('/>')


def format_file(input_fname: str,
                indent_size: int = 2,
                indent_style: bool = False,
                block_separation_max_depth: int = 2,
                alphebitize_attributes: bool = False,
                close_style: bool = False,
                namespace: bool = False) -> None:
    """Script to format xml files

    Args:
        input_fname (str): Input file name
        indent_size (int): Indent size
        indent_style (bool): Style of indentation (0=fixed, 1=hanging)
        block_separation_max_depth (int): Max depth to separate xml blocks
        alphebitize_attributes (bool): Alphebitize attributes
        close_style (bool): Style of close tag (0=same line, 1=new line)
        namespace (bool): Insert this namespace in the xml description
    """
    fname = os.path.expanduser(input_fname)
    try:
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        prologue_comments = [tmp.text for tmp in root.itersiblings(preceding=True)]
        epilog_comments = [tmp.text for tmp in root.itersiblings()]

        with open(fname, 'w') as f:
            f.write('<?xml version=\"1.0\" ?>\n')

            for comment in reversed(prologue_comments):
                f.write('\n<!--%s-->' % (comment))

            format_xml_level(f,
                             root,
                             0,
                             indent=' ' * indent_size,
                             block_separation_max_depth=block_separation_max_depth,
                             modify_attribute_indent=indent_style,
                             sort_attributes=alphebitize_attributes,
                             close_tag_newline=close_style,
                             include_namespace=namespace)

            for comment in epilog_comments:
                f.write('\n<!--%s-->' % (comment))
            f.write('\n')

    except ElementTree.ParseError as err:
        print('\nCould not load file: %s' % (fname))
        print(err.msg)
        raise Exception('\nCheck input file!')


def main() -> None:
    """Script to format xml files

    Args:
        input (str): Input file name
        -i/--indent (int): Indent size
        -s/--style (int): Indent style
        -d/--depth (int): Block separation depth
        -a/--alphebitize (int): Alphebitize attributes
        -c/--close (int): Close tag style
        -n/--namespace (int): Include namespace
    """
    parser = command_line_parsers.build_xml_formatter_input_parser()
    args = parser.parse_args()
    format_file(args.input,
                indent_size=args.indent,
                indent_style=args.style,
                block_separation_max_depth=args.depth,
                alphebitize_attributes=args.alphebitize,
                close_style=args.close,
                namespace=args.namespace)


if __name__ == "__main__":
    main()
