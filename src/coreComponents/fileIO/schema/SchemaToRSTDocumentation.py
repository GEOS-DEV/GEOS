
import os
import re
import numpy as np
from xml.etree import ElementTree as etree


class TreeBuilderWithComments(etree.TreeBuilder):
  def comment(self, data):
    self.start(etree.Comment, {})
    self.data(data)
    self.end(etree.Comment)


def writeTableRST(file_name, values):
  L = [[len(x) for x in row] for row in values]
  M = tuple(np.amax(np.array(L), axis=0))

  # Build column formats (the final column is allowed to use multiple lines)
  single_row_format = ('{:<%is} ' * len(M) + '\n') % M

  # The final column is allowed to span multiple lines
  description_index = sum(M[:-1]) + len(M) - 1
  multiple_row_format_a = ('{:<%is} ' * (len(M) - 1) + '| {:<%is} \n') % M
  multiple_row_format_b = ' ' * description_index + '| {:<%is} \n' % (M[-1])

  # Build headers
  boundary = single_row_format.format(*['=' * x for x in M])

  # Format lines
  formatted_lines = []
  for ii in range(0, len(values)):
    if ('\\n' not in values[ii][-1]):
      # Single line entry
      formatted_lines.append(single_row_format.format(*values[ii]))

    else:
      # Multi-line entry
      description = values[ii][-1].split('\\n')
      values[ii][-1] = description[0]
      formatted_lines.append(multiple_row_format_a.format(*values[ii]))

      for d in description[1:]:
        formatted_lines.append(multiple_row_format_b.format(d))

  # Build table
  with open(file_name, 'w') as f:
    f.write('\n\n')
    f.write(boundary)
    f.write(formatted_lines[0])
    f.write(boundary)
    for line in formatted_lines[1:]:
      f.write(line)
    f.write(boundary)
    f.write('\n\n')


def format_type_data(x):
  # Make sure that there are spaces within angled brackets
  # This allows for easier table line-breaks
  x = re.sub(r"<([a-zA-Z])", r"< \1", x)
  x = re.sub(r"([a-zA-Z])>", r"\1 >", x)

  # Remove any source tags
  x = re.sub(r"[a-zA-Z]*::", "", x)

  return(x)


def buildAttributeMap(root_node, xsd='{http://www.w3.org/2001/XMLSchema}'):
  attribute_map = {}

  # Build forward maps on the first pass (name, type, default, use, description, registered_by)
  for child_node in root_node.findall(xsd + 'complexType'):
    type_name = child_node.get('name')[:-4]
    attribute_map[type_name] = {}

    # Read the values supported by the xsd format (name, type, default, use)
    for attribute_node in child_node.findall(xsd + 'attribute'):
      att_name = attribute_node.get('name')
      att_type = format_type_data(attribute_node.get('type', default=' '))
      att_default = attribute_node.get('default', default=' ')

      # Handle the special case for required values
      useValue = attribute_node.get('use')
      if useValue:
        att_default = useValue

      attribute_map[type_name][att_name] = {'Type': att_type, 'Default': att_default}

    # Parse nodes
    for choice_node in child_node.findall(xsd + 'choice'):
      for sub_node in choice_node.findall(xsd + 'element'):
        sub_name = sub_node.get('name')
        sub_required = sub_node.get('minOccurs')
        sub_unique = sub_node.get('maxOccurs')
        node_use = ''
        if sub_required:
          if sub_unique:
            node_use = 'unique, required'
          else:
            node_use = 'required'
        elif (sub_unique):
          node_use = 'unique'

        attribute_map[type_name][sub_name] = {'Type': 'node',
                                              'Default': node_use}

    # Read the remaining values stored in comments (description, registered_by)
    for comment_node in child_node.iter(etree.Comment):
      tmp = comment_node.text.split(' => ')
      att_name = tmp[0]
      att_description = tmp[1].replace('\\\\', '\\').replace('\n', '\\n')
      attribute_map[type_name][att_name]['Description'] = att_description
      if (len(tmp) > 2):
       attribute_map[type_name][att_name]['Registered By'] = [x for x in tmp[2].split(', ')]

  # Build the reverse maps on the second pass (registered_on)
  for ka in attribute_map.keys():
    for kb in attribute_map[ka].keys():
      if ('Registered By' in attribute_map[ka][kb]):
        for kc in attribute_map[ka][kb]['Registered By']:
          attribute_map[kc][kb] = attribute_map[ka][kb]
          attribute_map[kc][kb]['Registered On'] = [kc]

  return attribute_map


def buildTableValues(type_map,
                     link_string='XML'):
  table_headers = ['Name', 'Type', 'Default', 'Registered On', 'Registered By', 'Description']
  optional_headers = ['Default', 'Registered On', 'Registered By']
  akeys = type_map[type_name].keys()

  # Remove any unused headers
  for ii in range(0, len(optional_headers)):
    header_count = 0
    for k in akeys:
      if (optional_headers[ii] in type_map[type_name][k]):
        header_count += 1
    if (header_count == 0):
      table_headers.remove(optional_headers[ii])

  # Setup the empty table
  N_rows = len(akeys)
  N_cols = len(table_headers)
  table_values = [[' ' for jj in N_cols] for ii in range(0, N_rows+1)]
  table_values[0] = table_headers

  # Add values to the table
  for ii in range(0, N_rows):
    # Get the row, set the name
    att_name = akeys[ii]
    table_row = table_values[ii+1]
    table_row[0] = att_name

    # Set the other parameters
    for jj in range(1, len(table_headers)):
      k = table_headers[jj]
      if k in type_map[type_name][att_name]:
        table_row[jj] = type_map[type_name][att_name][k]

        # Format any registration entries as links
        if ('Registered' in k):
          table_row[jj] = ", ".join([':ref:`%s_%s`' % (link_string, x) for x in table_row[jj]])

    # Set the link if the target is a node
    if (table_row[1] == 'node'):
      table_row[-1] = ':ref:`%s_%s`' % (link_string, table_row[0])

  return table_values


# Config
schema_name = 'schema.xsd'
additional_documentation_name = 'schema.xsd.other'
complete_output = '../../../docs/sphinx/CompleteXMLSchema'
output_folder = 'docs'
sphinx_path = '../../coreComponents/fileIO/schema/docs'
xsd = '{http://www.w3.org/2001/XMLSchema}'


# Parse the input/non-input schemas
parser = etree.XMLParser(target=TreeBuilderWithComments())
include_tree = etree.parse(schema_name, parser=parser)
include_root = include_tree.getroot()
input_attribute_map = buildAttributeMap(include_root)

parser = etree.XMLParser(target=TreeBuilderWithComments())
include_tree = etree.parse(additional_documentation_name, parser=parser)
include_root = include_tree.getroot()
other_attribute_map = buildAttributeMap(include_root)


# TODO: handle non-unique (and case-sensitive) links


# Build documentation tables
os.system('mkdir -p %s' % (output_folder))
with open('%s.rst' % (complete_output), 'w') as output_handle:
  # Write the file header
  output_handle.write('======================\n')
  output_handle.write('GEOSX Data Structure\n')
  output_handle.write('======================\n\n')

  # Parse the input schema definitions
  output_handle.write('********************************\n')
  output_handle.write('Input Schema Definitions\n')
  output_handle.write('********************************\n\n')

  for type_name in input_attribute_map.keys():
    # Write the individual tables
    table_values = buildTableValues(input_attribute_map[type_name])
    writeTableRST('%s/%s.rst' % (output_folder, type_name), table_values)

    # Write to the master list
    element_header = 'Element: %s' % (type_name)
    output_handle.write('\n.. _XML_%s:\n\n' % (type_name))
    output_handle.write('%s\n' % (element_header))
    output_handle.write('='*len(element_header) + '\n')
    output_handle.write('.. include:: %s/%s.rst\n\n' % (sphinx_path, type_name))

  # Parse the non-input schema definitions
  output_handle.write('********************************\n')
  output_handle.write('Datastructure Definitions\n')
  output_handle.write('********************************\n\n')

  for type_name in other_attribute_map.keys():
    # Write the individual tables
    table_values = buildTableValues(other_attribute_map[type_name], link_string='DATASTRUCTURE')
    writeTableRST('%s/%s_other.rst' % (output_folder, type_name), table_values)

    # Write to the master list
    element_header = 'Datastructure: %s' % (type_name)
    output_handle.write('\n.. _DATASTRUCTURE_%s:\n\n' % (type_name))
    output_handle.write('%s\n' % (element_header))
    output_handle.write('='*len(element_header) + '\n')
    output_handle.write('.. include:: %s/%s_other.rst\n\n' % (sphinx_path, type_name))


