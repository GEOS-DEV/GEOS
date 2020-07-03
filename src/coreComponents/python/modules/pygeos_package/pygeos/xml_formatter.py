import argparse
import os
from lxml import etree as ElementTree


def format_xml_level(output,
                     node,
                     level,
                     indent=' '*2,
                     block_separation_max_depth=2,
                     modify_attribute_indent=False,
                     sort_attributes=False,
                     close_tag_newline=False,
                     include_namespace=False):
  """Iteratively format the xml file

     @param output the output filename
     @param node the current xml element
     @param level the xml depth
     @param indent the xml indent style
     @param block_separation_max_depth the maximum depth to separate adjacent elements
     @param modify_attribute_indent option to have flexible attribute indentation
     @param sort_attributes option to sort attributes alphabetically
     @param close_tag_newline option to place close tag on a separate line
     @param include_namespace option to include the xml namespace in the output
  """

  # Handle comments
  if node.tag is ElementTree.Comment:
    output.write('\n%s<!--%s-->' % (indent*level, node.text))

  else:
    # Write opening line
    opening_line = '\n%s<%s' % (indent*level, node.tag)
    output.write(opening_line)

    # Write attributes
    if (len(node.attrib) > 0):
      # Choose indentation
      attribute_indent = '%s' % (indent*(level+1))
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
        format_xml_level(output, child, level+1, indent,
                         block_separation_max_depth,
                         modify_attribute_indent,
                         sort_attributes,
                         close_tag_newline,
                         include_namespace)

        # Add space between blocks
        if ((level < block_separation_max_depth) & (ii < Nc-1) & (child.tag is not ElementTree.Comment)):
          output.write('\n')

      # Write the end tag
      output.write('\n%s</%s>' % (indent*level, node.tag))
    else:
      if close_tag_newline:
        output.write('\n%s/>' % (indent*level))
      else:
        output.write('/>')


def main():
  """Script to format xml files

  @arg input Input file name
  """

  # Parse the user arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('input', type=str, help='Input file name')
  parser.add_argument('-i', '--indent', type=int, help='Indent size', default=2)
  parser.add_argument('-s', '--style', type=int, help='Indent style', default=0)
  parser.add_argument('-d', '--depth', type=int, help='Block separation depth', default=2)
  parser.add_argument('-a', '--alphebitize', type=int, help='Alphebetize attributes', default=0)
  parser.add_argument('-c', '--close', type=int, help='Close tag style', default=0)
  parser.add_argument('-n', '--namespace', type=int, help='Include namespace', default=0)
  args = parser.parse_args()

  # Process the xml file
  fname = os.path.expanduser(args.input)
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
                       indent=' '*args.indent,
                       block_separation_max_depth=args.depth,
                       modify_attribute_indent=args.style,
                       sort_attributes=args.alphebitize,
                       close_tag_newline=args.close,
                       include_namespace=args.namespace)

      for comment in epilog_comments:
        f.write('\n<!--%s-->' % (comment))
      f.write('\n')

  except ElementTree.ParseError as err:
    print('\nCould not load file: %s' % (fname))
    print(err.msg)
    raise Exception('\nCheck input file!')


if __name__ == "__main__":
  main()


