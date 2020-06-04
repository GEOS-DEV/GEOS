import argparse
import os
from xml.etree import ElementTree


def format_xml_level(output,
                     node,
                     level,
                     indent=' '*2,
                     block_separation_max_depth=2,
                     modify_attribute_indent=False,
                     sort_attributes=False,
                     close_tag_newline=False):
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

      # Sort attribute names
      akeys = list(node.attrib.keys())
      if sort_attributes:
        akeys = sorted(akeys)

      for ii in range(0, len(akeys)):
        k = akeys[ii]
        if ((ii == 0) & modify_attribute_indent):
          output.write(' %s=\"%s\"' % (k, node.attrib[k]))
        else:
          output.write('\n%s%s=\"%s\"' % (attribute_indent, k, node.attrib[k]))

    # Write children
    if len(node):
      output.write('>')
      Nc = len(node)
      for ii, child in zip(range(Nc), node):
        format_xml_level(output, child, level+1)

        # Add space between blocks
        if ((level < block_separation_max_depth) & (ii < Nc-1)):
          output.write('\n')

      # Write the end tag
      output.write('\n%s</%s>' % (indent*level, node.tag))
    else:
      if close_tag_newline:
        output.write('\n%s/>' % (indent*level))
      else:
        output.write('/>')


# Class to handle commented xml structure
class CommentedTreeBuilder(ElementTree.TreeBuilder):
    def comment(self, data):
        self.start(ElementTree.Comment, {})
        self.data(data)
        self.end(ElementTree.Comment)


def main():
  """Script to format xml files

  @arg input Input file name
  """

  # Parse the user arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('input', type=str, help='Input file name')
  args = parser.parse_args()

  # Process the xml file
  fname = os.path.expanduser(args.input)
  try:
    xml_parser = parser = ElementTree.XMLParser(target=CommentedTreeBuilder())
    tree = ElementTree.parse(fname, xml_parser)
    root = tree.getroot()

    with open(fname, 'w') as f:
      f.write('<?xml version=\"1.0\" ?>\n')
      format_xml_level(f, root, 0)
      f.write('\n')

  except ElementTree.XMLSyntaxError as err:
    print('\nCould not load file: %s' % (fname))
    print(err.msg)
    raise Exception('\nCheck input file!')


if __name__ == "__main__":
  main()


