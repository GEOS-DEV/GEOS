
import meshio
from meshio._mesh import CellBlock
import numpy as np
import argparse
import logging
import sys


def expand_nodesets_gmsh(input_mesh, output_mesh, logger=None, nodesets=[], depth=1):
    """
    @brief Expand nodesets in a gmsh file to include neighbors
    @param input_mesh path of the input abaqus file
    @param output_mesh path of the output gmsh file
    @param logger an instance of logging.Logger
    @param nodesets A list of nodeset names to expand
    @param depth The maximum depth
    """
    # Initialize the logger if it is empty
    if not logger:
        logging.basicConfig(level=logging.WARNING)
        logger = logging.getLogger(__name__)

    # Load the mesh
    n_warnings = 0
    logger.info('Reading the initial mesh...')
    mesh = meshio.read(input_mesh)

    region_list = list(mesh.field_data.keys())
    tri_nodeset_ids = mesh.cell_data['gmsh:physical'][1]
    all_tris = mesh.cells[1].data
    all_tetra = mesh.cells[0].data
    tetra_geometrical_ids = mesh.cell_data['gmsh:geometrical'][0]

    logger.info('Expanding nodesets...')
    new_tris, tri_nodeset, tri_region = [], [], []

    for nodeset_id, nodeset_name in enumerate(nodesets):
        if isinstance(nodeset_name, list):
            nodeset_name = nodeset_name[0]
        if nodeset_name not in region_list:
            print('Available nodesets:')
            print(region_list)
            raise Exception('Could not find target nodeset: %s' % (nodeset_name))

        # Locate nodesets
        logger.info('  %s' % (nodeset_name))
        initial_nodeset_id = mesh.field_data[nodeset_name][0]
        new_nodesest_name = '%s_expanded' % (nodeset_name)
        new_nodeset_id = len(region_list) + nodeset_id + 1
        mesh.field_data[new_nodesest_name] = [new_nodeset_id, 2]
        Ia = np.where(tri_nodeset_ids == initial_nodeset_id)[0]

        # Search for adjacent tets
        nodes = np.unique(np.reshape(all_tris[Ia], (-1)))
        adjacent_tetra = []
        for level in range(depth):
            logger.info('    (searching for adjacent tets, level=%i)' % (level))
            for ii in nodes:
                row, col = np.where(all_tetra == ii)
                adjacent_tetra.extend(row)

            if (level < depth - 1):
                for ii in adjacent_tetra:
                    nodes = np.append(nodes, all_tetra[ii, :])
                nodes = np.unique(nodes)

        # Build matching tris
        adjacent_tetra = np.unique(adjacent_tetra)
        adjacent_tris = []
        adjacent_tri_geom_ids = []
        logger.info('    (building trimesh)')
        for ii in adjacent_tetra:
            tmp = sorted(all_tetra[ii, :])
            adjacent_tris.append(tmp[:-1])
            adjacent_tris.append(tmp[1:])
            adjacent_tri_geom_ids.append(tetra_geometrical_ids[ii])
            adjacent_tri_geom_ids.append(tetra_geometrical_ids[ii])

        adjacent_tris, Ia = np.unique(adjacent_tris, axis=0, return_index=True)
        adjacent_tri_geom_ids = np.array(adjacent_tri_geom_ids)[Ia]
        adjacent_tri_ids = np.zeros(len(adjacent_tris)) + new_nodeset_id

        new_tris.append(adjacent_tris)
        tri_nodeset.append(adjacent_tri_geom_ids)
        tri_region.append(adjacent_tri_ids)

    # Add new tris
    if new_tris:
        new_tris = np.concatenate(new_tris, axis=0)
        tri_nodeset = np.concatenate(tri_nodeset, axis=0)
        tri_region = np.concatenate(tri_region, axis=0)

        logger.info('  Adding %i new triangles...' % (len(new_tris)))
        if (-1 in tri_region):
            logger.warning('Triangles with empty region information found!')
            logger.warning('Note: These will be indicated by a -1 in the output file.')
            n_warnings += 1
        mesh.cells.append(CellBlock('triangle', np.array(new_tris)))
        mesh.cell_data['gmsh:geometrical'].append(np.array(tri_region))
        mesh.cell_data['gmsh:physical'].append(np.array(tri_nodeset))

    # Write the final mesh
    logger.info('Writing gmsh mesh...')
    meshio.write(output_mesh, mesh, file_format="gmsh22", binary=False)
    logger.info('Done!')

    return (n_warnings > 0)


def main():
    """
    @brief Expand a nodeset to contain neighboring elements
    @arg input Input gmsh file name
    @arg output Output gmsh file name
    @arg nodeset Name of the target nodeset
    @arg depth Include nodes within this number of connections from the target
    """

    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input gmsh mesh file name')
    parser.add_argument('output', type=str, help='Output gmsh mesh file name')
    parser.add_argument('-n', '--nodeset', help='Target nodeset', nargs='+', action='append', default=[])
    parser.add_argument('-d', '--depth', type=int, help='Maximum connection depth', default=1)
    parser.add_argument('-v', '--verbose', help='Increase verbosity level', action="store_true")
    args = parser.parse_args()

    # Set up a logger
    logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger(__name__)
    if args.verbose:
        logger.setLevel(logging.INFO)

    # Call the converter
    err = expand_nodesets_gmsh(args.input, args.output, nodesets=args.nodeset, depth=args.depth, logger=logger)
    if err:
        sys.exit('Warnings detected: check the output file for potential errors!')

