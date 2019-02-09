
from lxml import etree as ElementTree


def format_xml_file(input_file, output_file=''):

  active_blocks = ['Solvers',
                   'Outputs',
                   'Events',
                   'Mesh',
                   'Partition',
                   'NumericalMethods',
                   'ElementRegions',
                   'Constitutive',
                   'InitialConditions',
                   'FieldSpecifications',
                   'Functions',
                   'Geometry',
                   'Fracture']

  documentation_blocks = ['Included',
                          'Parameters']


  # Open new files
  parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
  tree_old = ElementTree.parse(input_file, parser=parser)
  root_old = tree_old.getroot()

  tree_new = ElementTree.parse(input_file, parser=parser)
  root_new = tree_new.getroot()
  for child in root_new:
    root_new.remove(child)

  # Insert elements in the desired order
  inserted_nodes = []
  comment = ElementTree.Comment('         End            ')
  root_new.insert(-1, comment)

  # Active blocks
  comment = ElementTree.Comment('                        ')
  root_new.insert(-1, comment)
  comment = ElementTree.Comment('     Active Blocks      ')
  root_new.insert(-1, comment)
  for block in active_blocks:
    for child in root_old.findall(block):
      comment = ElementTree.Comment('                        ')
      root_new.insert(-1, comment)
      root_new.insert(-1, child)
      inserted_nodes.append(block)

  # Documentation blocks
  comment = ElementTree.Comment('                        ')
  root_new.insert(-1, comment)
  comment = ElementTree.Comment('  Documentation Blocks  ')
  root_new.insert(-1, comment)
  for block in documentation_blocks:
    for child in root_old.findall(block):
      comment = ElementTree.Comment('                        ')
      root_new.insert(-1, comment)
      root_new.insert(-1, child)
      inserted_nodes.append(block)

  # Check for unsorted blocks
  unsorted_blocks = []
  for child in root_old:
    if (child.tag not in inserted_nodes):
      unsorted_blocks.append(child.tag)

  if unsorted_blocks:
    comment = ElementTree.Comment('                        ')
    root_new.insert(-1, comment)
    comment = ElementTree.Comment('    Unsorted Blocks     ')
    root_new.insert(-1, comment)
    for block in unsorted_blocks:
      for child in root_old.findall(block):
        comment = ElementTree.Comment('                        ')
        root_new.insert(-1, comment)
        root_new.insert(-1, child)
        inserted_nodes.append(block)

  # Write output
  if output_file:
    tree_new.write(output_file, pretty_print=True)
  else:
    tree_new.write(input_file, pretty_print=True)

