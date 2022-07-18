import meshio
from meshio._mesh import CellBlock
import numpy as np
import argparse
import logging
import sys


def convert_abaqus_to_gmsh(input_mesh, output_mesh, logger=None):
    """
    @brief Convert an abaqus mesh to gmsh 2 format, preserving nodeset information.
    @details If the code encounters any issues with region/element indices,
             the conversion will attempt to continue, with errors
             indicated by -1 values in the output file.
    @param input_mesh path of the input abaqus file
    @param output_mesh path of the output gmsh file
    @param logger an instance of logging.Logger
    """
    # Initialize the logger if it is empty
    if not logger:
        logging.basicConfig(level=logging.WARNING)
        logger = logging.getLogger(__name__)

    # Keep track of the number of warnings
    n_warnings = 0

    # Load the mesh
    logger.info('Reading abaqus mesh...')
    mesh = meshio.read(input_mesh, file_format="abaqus")

    # Convert the element regions to tags
    logger.info('Converting region tags...')
    region_list = list(mesh.cell_sets.keys())
    n_regions = len(region_list)
    cell_ids = []
    for block_id, block in enumerate(mesh.cells):
        cell_ids.append(np.zeros(len(block[1]), dtype=int) - 1)
        for region_id, region in enumerate(region_list):
            mesh.field_data[region] = [region_id + 1, 3]
            cell_ids[block_id][mesh.cell_sets[region][block_id]] = region_id + 1

        # Check for bad element region conversions
        if (-1 in cell_ids[-1]):
            logger.warning('Some element regions in block %i did not convert correctly to tags!' % (block_id))
            logger.warning('Note: These will be indicated by a -1 in the output file.')
            n_warnings += 1

    # Add to the meshio datastructure
    # Note: the copy here is required, so that later appends
    #       do not break these dicts
    mesh.cell_data['gmsh:physical'] = cell_ids.copy()
    mesh.cell_data['gmsh:geometrical'] = cell_ids.copy()

    # Build the face elements
    logger.info('Converting nodesets to face elements, tags...')
    new_tris, tri_nodeset, tri_region = [], [], []
    new_quads, quad_nodeset, quad_region = [], [], []

    for nodeset_id, nodeset_name in enumerate(mesh.point_sets):
        logger.info('  %s' % (nodeset_name))
        mesh.field_data[nodeset_name] = [nodeset_id + n_regions + 1, 2]
        nodeset = mesh.point_sets[nodeset_name]

        # Search by block, then element
        for block_id, block in enumerate(mesh.cells):
            for element_id, element in enumerate(block[1]):
                # Find any matching nodes
                matching_nodes = [x for x in element if x in nodeset]

                # Add a new face element if there are enough nodes
                n_matching = len(matching_nodes)
                if (n_matching >= 3):
                    # Find the region
                    region_id = -1
                    for region in region_list:
                        if (element_id in mesh.cell_sets[region][block_id]):
                            region_id = mesh.field_data[region][block_id]

                    # Test to see if the element is a quad or triangle
                    tag_id = mesh.field_data[nodeset_name][0]
                    if (n_matching == 3):
                        new_tris.append(matching_nodes)
                        tri_nodeset.append(tag_id)
                        tri_region.append(region_id)

                    elif (n_matching == 4):
                        new_quads.append(matching_nodes)
                        quad_nodeset.append(tag_id)
                        quad_region.append(region_id)

                    else:
                        logger.warning('  Discarding an element with an unexpected number of nodes')
                        logger.warning('    n_nodes=%i, element=%i, set=%s' % (n_matching, element_id, nodeset_name))
                        n_warnings += 1

    # Add new tris
    if new_tris:
        logger.info('  Adding %i new triangles...' % (len(new_tris)))
        if (-1 in tri_region):
            logger.warning('Triangles with empty region information found!')
            logger.warning('Note: These will be indicated by a -1 in the output file.')
            n_warnings += 1
        mesh.cells.append(CellBlock('triangle', np.array(new_tris)))
        mesh.cell_data['gmsh:geometrical'].append(np.array(tri_region))
        mesh.cell_data['gmsh:physical'].append(np.array(tri_nodeset))

    # Add new quads
    if new_quads:
        logger.info('  Adding %i new quads...' % (len(new_quads)))
        if (-1 in quad_region):
            logger.warning('Quads with empty region information found!')
            logger.warning('Note: These will be indicated by a -1 in the output file.')
            n_warnings += 1
        mesh.cells.append(CellBlock('quad', np.array(new_quads)))
        mesh.cell_data['gmsh:geometrical'].append(np.array(quad_region))
        mesh.cell_data['gmsh:physical'].append(np.array(quad_nodeset))

    # Write the final mesh
    logger.info('Writing gmsh mesh...')
    meshio.write(output_mesh, mesh, file_format="gmsh22", binary=False)
    logger.info('Done!')

    return (n_warnings > 0)


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
    parser.add_argument('-v', '--verbose', help='Increase verbosity level', action="store_true")
    args = parser.parse_args()

    # Set up a logger
    logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger(__name__)
    if args.verbose:
        logger.setLevel(logging.INFO)

    # Call the converter
    err = convert_abaqus_to_gmsh(args.input, args.output, logger)
    if err:
        sys.exit('Warnings detected: check the output file for potential errors!')
