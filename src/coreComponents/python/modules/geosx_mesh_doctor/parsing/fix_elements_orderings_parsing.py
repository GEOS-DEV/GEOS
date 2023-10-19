import logging
import random

from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAGONAL_PRISM,
    VTK_HEXAHEDRON,
    VTK_PENTAGONAL_PRISM,
    VTK_PYRAMID,
    VTK_TETRA,
    VTK_VOXEL,
    VTK_WEDGE,
)

from checks.fix_elements_orderings import Options, Result

from . import vtk_output_parsing, FIX_ELEMENTS_ORDERINGS


__CELL_TYPE_MAPPING = {
    "Hexahedron": VTK_HEXAHEDRON,
    "Prism5": VTK_PENTAGONAL_PRISM,
    "Prism6": VTK_HEXAGONAL_PRISM,
    "Pyramid": VTK_PYRAMID,
    "Tetrahedron": VTK_TETRA,
    "Voxel": VTK_VOXEL,
    "Wedge": VTK_WEDGE,
}

__CELL_TYPE_SUPPORT_SIZE = {
    VTK_HEXAHEDRON: 8,
    VTK_PENTAGONAL_PRISM: 10,
    VTK_HEXAGONAL_PRISM: 12,
    VTK_PYRAMID: 5,
    VTK_TETRA: 4,
    VTK_VOXEL: 8,
    VTK_WEDGE: 6,
}


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(FIX_ELEMENTS_ORDERINGS,
                              help="Reorders the support nodes for the given cell types.")
    for key, vtk_key in __CELL_TYPE_MAPPING.items():
        tmp = list(range(__CELL_TYPE_SUPPORT_SIZE[vtk_key]))
        random.Random(4).shuffle(tmp)
        p.add_argument('--' + key,
                       type=str,
                       metavar=",".join(map(str, tmp)),
                       default=None,
                       required=False,
                       help=f"[list of integers]: node permutation for \"{key}\".")
    vtk_output_parsing.fill_vtk_output_subparser(p)


def convert(parsed_options) -> Options:
    """
    From the parsed cli options, return the converted options for self intersecting elements check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    cell_type_to_ordering = {}
    for key, vtk_key in __CELL_TYPE_MAPPING.items():
        raw_mapping = parsed_options[key]
        if raw_mapping:
            tmp = tuple(map(int, raw_mapping.split(",")))
            if not set(tmp) == set(range(__CELL_TYPE_SUPPORT_SIZE[vtk_key])):
                err_msg = f"Permutation {raw_mapping} for type {key} is not valid."
                logging.error(err_msg)
                raise ValueError(err_msg)
            cell_type_to_ordering[vtk_key] = tmp
    vtk_output = vtk_output_parsing.convert(parsed_options)
    return Options(vtk_output=vtk_output,
                   cell_type_to_ordering=cell_type_to_ordering)


def display_results(options: Options, result: Result):
    if result.output:
        logging.info(f"New mesh was written to file '{result.output}'")
        if result.unchanged_cell_types:
            logging.info(f"Those vtk types were not reordered: [{', '.join(map(str, result.unchanged_cell_types))}].")
        else:
            logging.info("All the cells of the mesh were reordered.")
    else:
        logging.info("No output file was written.")
