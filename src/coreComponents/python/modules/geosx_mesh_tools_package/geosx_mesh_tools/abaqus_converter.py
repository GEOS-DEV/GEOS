
import meshio
from meshio._mesh import CellBlock
import numpy as np
import argparse


def convert_abaqus_to_gmsh(input_mesh, output_mesh):
  """
  @brief Convert an abaqus mesh to gmsh 2 format, preserving nodeset information
  @param input_mesh path of the input abaqus file
  @param output_mesh path of the output gmsh file
  """
  # Load the mesh
  print('Reading abaqus mesh...')
  mesh = meshio.Mesh.read(input_mesh, "abaqus")

  # Add a default region tag
  print('Converting region tags...')
  region_list = list(mesh.cell_sets.keys())
  Nregion = len(region_list)
  cell_ids = np.zeros(len(mesh.cells[0][1]), dtype=int)
  for ii, ka in zip(range(Nregion), region_list):
    mesh.field_data[ka] = [ii + 1, 3]
    region_list.append(ka)
    cell_ids[mesh.cell_sets[ka][0]] = ii + 1
  mesh.cell_data['gmsh:physical'] = [cell_ids]
  mesh.cell_data['gmsh:geometrical'] = [cell_ids]

  # Build the face elements
  new_tris = []
  tri_nodeset = []
  tri_region = []
  new_quads = []
  quad_nodeset = []
  quad_region = []
  print('Converting nodesets to face elements...')
  for ii, ka in zip(range(len(mesh.point_sets)), mesh.point_sets.keys()):
    print('  %s' % (ka))
    nodeset = mesh.point_sets[ka]

    # Search the elements
    for jj in range(len(mesh.cells[0][1])):
      # Find any matching nodes
      matching_nodes = []
      for kk in mesh.cells[0][1][jj]:
        if kk in nodeset:
          matching_nodes.append(kk)

      # Find the region
      region_id = -1
      if (len(matching_nodes) >= 3):
        for kb in region_list:
          if (jj in mesh.cell_sets[kb][0]):
            region_id = mesh.field_data[kb][0]

      # Test to see if they match a quad or triangle
      if (len(matching_nodes) == 3):
        new_tris.append(matching_nodes)
        tri_nodeset.append(ii + Nregion + 1)
        tri_region.append(region_id)

      elif (len(matching_nodes) == 4):
        new_quads.append(matching_nodes)
        quad_nodeset.append(ii + Nregion + 1)
        quad_region.append(region_id)

  # Add the nodeset tags
  for ii, ka in zip(range(len(mesh.point_sets)), mesh.point_sets.keys()):
    mesh.field_data[ka] = [ii + Nregion + 1, 2]

  # Add new tris
  if len(new_tris):
    print('  Adding %i new triangles...' % (len(new_tris)))
    mesh.cells.append(CellBlock('triangle', np.array(new_tris)))
    mesh.cell_data['gmsh:geometrical'].append(np.array(tri_region))
    mesh.cell_data['gmsh:physical'].append(np.array(tri_nodeset))

  # Add  new quads
  if len(new_quads):
    print('  Adding %i new quads...' % (len(new_quads)))
    mesh.cells.append(CellBlock('quad', np.array(new_quads)))
    mesh.cell_data['gmsh:geometrical'].append(np.array(quad_region))
    mesh.cell_data['gmsh:physical'].append(np.array(quad_nodeset))

  # Write the final mesh
  print('Writing gmsh mesh...')
  meshio.write(output_mesh, mesh, file_format="gmsh22", binary=False)

  print('Done!')


def main():
  """
  @brief Entry point for the abaqus convertor console script
  @arg input_mesh Input abaqus file name
  @arg output_mesh Output gmsh file name
  """

  # Parse the user arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('input', type=str, help='Input abaqus mesh file name')
  parser.add_argument('output', type=str, help='Output gmsh mesh file name')
  args = parser.parse_args()

  convert_abaqus_to_gmsh(args.input_mesh, args.output_mesh)


