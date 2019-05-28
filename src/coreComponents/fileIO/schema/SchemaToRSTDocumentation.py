
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


def format_value(x):
  # Make sure that there are spaces within angled brackets
  # This allows for easier table line-breaks
  x = re.sub(r"<([a-zA-Z])", r"< \1", x)
  x = re.sub(r"([a-zA-Z])>", r"\1 >", x)

  # Remove any source tags
  x = re.sub(r"[a-zA-Z]*::", "", x)

  return(x)


def parseSchemaNode(node,
                    link_string='XML',
                    include_defaults=True,
                    include_registrant=False):
  type_name = child_node.get('name')[:-4]
  table_headers = ['Name', 'Type', 'Default', 'Description', 'Registered By']
  table_dict = {k: [] for k in table_headers}

  # Parse comments
  attribute_comments = {}
  attribute_registrants = {}
  for comment_node in child_node.iter(etree.Comment):
    tmp = comment_node.text.split(' => ')
    attribute_comments[tmp[0]] = tmp[1].replace('\\\\', '\\').replace('\n', '\\n')
    if (len(tmp) > 2):
      attribute_registrants[tmp[0]] = ", ".join([':ref:`%s_%s`' % (link_string, x) for x in tmp[2].split(', ')])

  # Parse attributes
  for attribute_node in child_node.findall(xsd + 'attribute'):
    # Read the name, type, default values
    table_dict['Name'].append(format_value(attribute_node.get('name', default=' ')))
    table_dict['Type'].append(format_value(attribute_node.get('type', default=' ')))
    table_dict['Default'].append(format_value(attribute_node.get('default', default=' ')))

    # Handle the special case for required values
    useValue = attribute_node.get('use')
    if useValue:
      table_dict['Default'][-1] = useValue

    # Add any available descriptions, registrant information
    ka = table_dict['Name'][-1]
    if ka in attribute_comments:
      table_dict['Description'].append(attribute_comments[ka])
    else:
      table_dict['Description'].append('')

    if ka in attribute_registrants:
      table_dict['Registered By'].append(attribute_registrants[ka])
    else:
      table_dict['Registered By'].append('')

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

      table_dict['Name'].append(sub_name)
      table_dict['Type'].append('node')
      table_dict['Default'].append(node_use)
      table_dict['Description'].append(':ref:`%s_%s`' % (link_string, sub_name))
      table_dict['Registered By'].append('')

  # Handle empty tables
  if (len(table_dict['Name']) == 0):
    table_dict['Name'].append('')
    table_dict['Type'].append('')
    table_dict['Default'].append('')
    table_dict['Description'].append('(no documentation available)')
    table_dict['Registered By'].append('')

  # Remove unused columns
  if not include_defaults:
    table_headers.remove('Default')
  if (not any(table_dict['Registered By'])) | (not include_registrant):
    table_headers.remove('Registered By')

  # Move values into a table
  table_values = [[k for k in table_headers]]
  for ii in range(0, len(table_dict['Name'])):
    table_values.append([table_dict[k][ii] for k in table_headers])

  return type_name, table_values


# Config
schema_name = 'schema.xsd'
additional_documentation_name = 'schema.xsd.other'
complete_output = '../../../docs/sphinx/CompleteXMLSchema'
output_folder = 'docs'
sphinx_path = '../../coreComponents/fileIO/schema/docs'
xsd = '{http://www.w3.org/2001/XMLSchema}'


# Parse the schema, build documentation tables
os.system('mkdir -p %s' % (output_folder))

with open('%s.rst' % (complete_output), 'w') as output_handle:
  # Write the file header
  output_handle.write('======================\n')
  output_handle.write('GEOSX Data Structure\n')
  output_handle.write('======================\n\n')

  # Parse the schema definitions
  output_handle.write('********************************\n')
  output_handle.write('Input Schema Definitions\n')
  output_handle.write('********************************\n\n')

  parser = etree.XMLParser(target=TreeBuilderWithComments())
  include_tree = etree.parse(schema_name, parser=parser)
  include_root = include_tree.getroot()

  for child_node in include_root.findall(xsd + 'complexType'):
    # Parse schema into a list of key parameters
    type_name, table_values = parseSchemaNode(child_node)

    # Write table
    writeTableRST('%s/%s.rst' % (output_folder, type_name), table_values)

    element_header = 'Element: %s' % (type_name)
    output_handle.write('\n.. _XML_%s:\n\n' % (type_name))
    output_handle.write('%s\n' % (element_header))
    output_handle.write('='*len(element_header) + '\n')
    output_handle.write('.. include:: %s/%s.rst\n\n' % (sphinx_path, type_name))

  # Parse the non-schema definitions
  output_handle.write('********************************\n')
  output_handle.write('Datastructure Definitions\n')
  output_handle.write('********************************\n\n')

  parser = etree.XMLParser(target=TreeBuilderWithComments())
  include_tree = etree.parse(additional_documentation_name, parser=parser)
  include_root = include_tree.getroot()
  all_links = []

  for child_node in include_root.findall(xsd + 'complexType'):
    # The additional documentation uses the same format as the schema
    type_name, table_values = parseSchemaNode(child_node, link_string='DATASTRUCTURE', include_defaults=False, include_registrant=True)
    type_name_lower = type_name.lower()

    # Write table
    writeTableRST('%s/%s_other.rst' % (output_folder, type_name), table_values)

    element_header = 'Datastructure: %s' % (type_name)

    # There are some names in the GEOSX datastructure that only vary by case.
    # This can break the documentation, so skip link tags if necessary
    if (type_name_lower not in all_links):
      output_handle.write('\n.. _DATASTRUCTURE_%s:\n\n' % (type_name))
      all_links.append(type_name.lower())

    output_handle.write('%s\n' % (element_header))
    output_handle.write('='*len(element_header) + '\n')
    output_handle.write('.. include:: %s/%s_other.rst\n\n' % (sphinx_path, type_name))


