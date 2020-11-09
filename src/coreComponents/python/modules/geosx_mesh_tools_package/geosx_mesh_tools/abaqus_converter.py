
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

    # Convert the element regions to tags
    print('Converting region tags...')
    region_list = list(mesh.cell_sets.keys())
    n_regions = len(region_list)
    cell_ids = np.zeros(len(mesh.cells[0][1]), dtype=int) - 1
    for ii, region in enumerate(region_list):
        mesh.field_data[region] = [ii + 1, 3]
        cell_ids[mesh.cell_sets[region][0]] = ii + 1

    # Add to the meshio datastructure
    mesh.cell_data['gmsh:physical'] = [cell_ids]
    mesh.cell_data['gmsh:geometrical'] = [cell_ids]
    if (-1 in cell_ids):
        raise Exception('Element regions did not convert correctly to tags!')

    # Build the face elements
    print('Converting nodesets to face elements, tags...')
    new_tris, tri_nodeset, tri_region = [], [], []
    new_quads, quad_nodeset, quad_region = [], [], []

    for nodeset_id, nodeset_name in enumerate(mesh.point_sets):
        print('  %s' % (nodeset_name))
        mesh.field_data[nodeset_name] = [nodeset_id + n_regions + 1, 2]
        nodeset = mesh.point_sets[nodeset_name]

        # Search the elements
        for element_id in range(len(mesh.cells[0][1])):
            # Find any matching nodes
            matching_nodes = []
            for node_id in mesh.cells[0][1][element_id]:
                if node_id in nodeset:
                    matching_nodes.append(node_id)

            # Find the region
            region_id = -1
            if (len(matching_nodes) >= 3):
                for region in region_list:
                    if (element_id in mesh.cell_sets[region][0]):
                        region_id = mesh.field_data[region][0]

            # Test to see if they match a quad or triangle
            tag_id = mesh.field_data[nodeset_name][0]
            n_matching = len(matching_nodes)
            if (n_matching == 3):
                new_tris.append(matching_nodes)
                tri_nodeset.append(tag_id)
                tri_region.append(region_id)

            elif (n_matching == 4):
                new_quads.append(matching_nodes)
                quad_nodeset.append(tag_id)
                quad_region.append(region_id)

            elif (n_matching > 4):
                raise Exception('Unexpected number of nodes (%i) for element %i in set: %s' % (n_matching, element_id, nodeset_name))

    # Add new tris
    if new_tris:
        print('  Adding %i new triangles...' % (len(new_tris)))
        if (-1 in tri_region):
            raise Exception('Traingles with empty region information found!')
        mesh.cells.append(CellBlock('triangle', np.array(new_tris)))
        mesh.cell_data['gmsh:geometrical'].append(np.array(tri_region))
        mesh.cell_data['gmsh:physical'].append(np.array(tri_nodeset))

    # Add  new quads
    if new_quads:
        print('  Adding %i new quads...' % (len(new_quads)))
        if (-1 in quad_region):
            raise Exception('Quads with empty region information found!')
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

    convert_abaqus_to_gmsh(args.input, args.output)


