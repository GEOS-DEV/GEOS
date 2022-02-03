
"""
visitVTKConversion.py

VisIt macros designed to convert GEOSX vtk format files containing multiple
regions into separate files, which can then be normally read by the tool.
To install the macros, do the following:
  1) Select Controls/Command from the top menu in VisIt
  2) In the newly opened Commands window, select the Macros tab
  3) Copy the following scripts to the panel and click Update Macros
  4) The macro will appear under the name 'Convert pvd' in the macros window

If the macros window is not visible, select Controls/Macros from the top menu in VisIt.
When you select 'Convert pvd', visit will attempt to convert any pvd files in the current
directory (where the visit command was ran).  The macro will create a new set of .vtm files
in the vtk output directory (named '[Region]_*.vtm').

For multi-region problems, you can load each set of .vtm files individual.
After adding plots to the interface, you may be prompted to 'Correlate databases'.
If you select 'Yes', VisIt will add a new time slider 'Correlation*' to the interface, which can
be used to change the visualization time, while keeping the datasets  syncronized.

Note: the VisIt python interpreter does not allow empty lines within functions or multi-line arguments
"""


def convert_vtm_file(vtm_file, output_dir, root_path):
  from xml.etree import ElementTree
  import os
  vtm_dir = os.path.split(vtm_file)[0]
  vtm_header = os.path.split(vtm_file)[1].split('.')[0]
  vtm_tree = ElementTree.parse(vtm_file)
  vtm_root = vtm_tree.getroot()
  multiblock = vtm_root[0]
  for block in multiblock:
    for child_block in block:
      # Setup a new tree and copy parent elements
      visit_root = ElementTree.Element(vtm_root.tag, attrib=vtm_root.attrib)
      visit_multiblock = ElementTree.Element(multiblock.tag, attrib=multiblock.attrib)
      visit_root.append(visit_multiblock)
      visit_block = ElementTree.Element(block.tag, attrib=block.attrib)
      visit_multiblock.append(visit_block)
      for dataset in child_block:
        # Copy the dataset, then correct the relative path
        visit_block.append(dataset)
        dataset_path = os.path.relpath(os.path.join(root_path, vtm_dir, dataset.get('file')), start=output_dir)
        dataset.set('file', dataset_path)
      visit_tree = ElementTree.ElementTree(element=visit_root)
      visit_tree.write(os.path.join(output_dir, '%s_%s.vtm' % (child_block.get('name'), vtm_header)))

def user_macro_convert_pvd_files():
  import glob
  import os
  from xml.etree import ElementTree
  print('Converting pvd files in current directory:')
  for file in glob.glob('*.pvd'):
    print('  ' + file)
    pvd_dir = os.path.split(file)[0]
    output_dir = file[:file.rfind('.')]
    pvd_tree = ElementTree.parse(file)
    pvd_root = pvd_tree.getroot()
    collection_root = pvd_root[0]
    for dataset in collection_root:
      convert_vtm_file(dataset.get('file'), output_dir, pvd_dir)
  print('Done!')
RegisterMacro("Convert pvd", user_macro_convert_pvd_files)
