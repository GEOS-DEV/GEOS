import logging
import textwrap
from typing import Dict

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

from . import cli_parsing, vtk_output_parsing, FIX_ELEMENTS_ORDERINGS


__TETRAHEDRON_NAME = "Tetrahedron"

__CELL_TYPE_MAPPING = {
    "Hexahedron": VTK_HEXAHEDRON,
    "Prism5": VTK_PENTAGONAL_PRISM,
    "Prism6": VTK_HEXAGONAL_PRISM,
    "Pyramid": VTK_PYRAMID,
    __TETRAHEDRON_NAME: VTK_TETRA,
    "Voxel": VTK_VOXEL,
    "Wedge": VTK_WEDGE,
}

__ALL_KEYWORDS = {
    *vtk_output_parsing.get_vtk_output_keywords(),
    *__CELL_TYPE_MAPPING.keys()
}


def get_help():
    msg = f"""\
    Reorders the support nodes for the given cell types.
    
    Supported cell types are '{", ".join(__CELL_TYPE_MAPPING.keys())}'.
    For example, use the '{__TETRAHEDRON_NAME}=1,2,3,0' option to change the orders of all the tetrahedra of the input mesh.
    {vtk_output_parsing.get_vtk_output_help()}
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for self intersecting elements check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    cell_type_support_size = {
        VTK_HEXAHEDRON: 8,
        VTK_PENTAGONAL_PRISM: 10,
        VTK_HEXAGONAL_PRISM: 12,
        VTK_PYRAMID: 5,
        VTK_TETRA: 4,
        VTK_VOXEL: 8,
        VTK_WEDGE: 6,
    }
    options: Dict[str, str] = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(FIX_ELEMENTS_ORDERINGS, __ALL_KEYWORDS, options)
    cell_type_to_ordering = {}
    for key, vtk_key in __CELL_TYPE_MAPPING.items():
        raw_mapping = options.get(key)
        if raw_mapping:
            tmp = tuple(map(int, raw_mapping.split(",")))
            if not set(tmp) == set(range(cell_type_support_size[vtk_key])):
                err_msg = f"Permutation {raw_mapping} for type {key} is not valid."
                logging.error(err_msg)
                raise ValueError(err_msg)
            cell_type_to_ordering[vtk_key] = tmp
    vtk_output = vtk_output_parsing.parse_cli_options(options)
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
