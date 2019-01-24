
import os
import numpy as np
from lxml import etree
import string


def writeTableRST(file_name, element_name, values):
  L = [[len(x) for x in row] for row in values]
  M = np.amax(np.array(L), axis=0)

  # Build headers
  boundary = ''
  for x in M:
    boundary += '='*x + ' '
  boundary += '\n'

  # Format lines
  formatted_lines = []
  for ii in range(0, len(values)):
    formatted_lines.append('')
    for jj in range(0, len(values[ii])):
      formatted_lines[ii] += string.ljust(values[ii][jj], M[jj]) + ' '
    formatted_lines[ii] += '\n'

  # Build table
  with open(file_name, 'w') as f:
    tmp = 'Element: %s' % (element_name)
    f.write('\n%s\n' % (tmp))
    f.write('='*len(tmp) + '\n\n')
    f.write(boundary)
    f.write(formatted_lines[0])
    f.write(boundary)
    for line in formatted_lines[1:]:
      f.write(line)
    f.write(boundary)
    f.write('\n\n')


# Config
input_name = 'schema.xsd'
complete_output = 'xml_schema'
output_folder = 'docs'
xsd = '{http://www.w3.org/2001/XMLSchema}'


# Parse schema, build docs
os.system('mkdir -p %s' % (output_folder))

include_tree = etree.parse(input_name)
include_root = include_tree.getroot()

with open('%s/%s.rst' % (output_folder, complete_output), 'w') as output_handle:
  output_handle.write('=============\n')
  output_handle.write('XML Structure\n')
  output_handle.write('=============\n')

  for child_node in include_root.findall(xsd + 'complexType'):
    type_name = child_node.get('name')[:-4]
    table_values = [['Name', 'Type', 'Default', 'Use', 'Description']]

    # Parse comments
    attribute_comments = {}
    for comment_node in child_node.iterchildren(etree.Comment):
      tmp = str(comment_node)[4:-3].split(' = ')
      attribute_comments[tmp[0]] = tmp[1]

    # Parse attributes
    for attribute_node in child_node.findall(xsd + 'attribute'):
      table_values.append([attribute_node.get(v, default=' ') for v in ['name', 'type', 'default', 'use']])
      k = table_values[-1][0]
      if k in attribute_comments:
        table_values[-1].append(attribute_comments[k])
      else:
        table_values[-1].append(' ')

    # Parse nodes
    for choice_node in child_node.findall(xsd + 'choice'):
      for sub_node in choice_node.findall(xsd + 'element'):
        sub_name = sub_node.get('name')
        sub_required = sub_node.get('minOccurs')
        sub_unique = sub_node.get('maxOccurs')
        use = ''
        if sub_required:
          if sub_unique:
            use = 'unique, required'
          else:
            use = 'required'
        elif (sub_unique):
          use = 'unique'

        table_values.append([sub_name, 'node', '', use, '`Element: %s`_' % (sub_name)])

    # Write table
    writeTableRST('%s/%s.rst' % (output_folder, type_name), type_name, table_values)
    output_handle.write('\n.. include:: ./%s.rst\n' % (type_name))

