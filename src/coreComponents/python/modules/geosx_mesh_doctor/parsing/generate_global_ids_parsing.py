from dataclasses import dataclass
import logging

from checks.generate_global_ids import Options, Result

from . import vtk_output_parsing, GENERATE_GLOBAL_IDS


__CELLS, __POINTS = "cells", "points"


@dataclass(frozen=True)
class GlobalIdsInfo:
    cells: bool
    points: bool


def convert_global_ids(parsed_options) -> GlobalIdsInfo:
    return GlobalIdsInfo(cells=parsed_options[__CELLS],
                         points=parsed_options[__POINTS])


def convert(parsed_options) -> Options:
    gids: GlobalIdsInfo = convert_global_ids(parsed_options)
    return Options(vtk_output=vtk_output_parsing.convert(parsed_options),
                   generate_cells_global_ids=gids.cells,
                   generate_points_global_ids=gids.points)


def fill_generate_global_ids_subparser(p):
    p.add_argument('--' + __CELLS,
                   action="store_true",
                   help=f"[bool]: Generate global ids for cells. Defaults to true.")
    p.add_argument('--no-' + __CELLS,
                   action="store_false",
                   dest=__CELLS,
                   help=f"[bool]: Don't generate global ids for cells.")
    p.set_defaults(**{__CELLS: True})
    p.add_argument('--' + __POINTS,
                   action="store_true",
                   help=f"[bool]: Generate global ids for points. Defaults to true.")
    p.add_argument('--no-' + __POINTS,
                   action="store_false",
                   dest=__POINTS,
                   help=f"[bool]: Don't generate global ids for points.")
    p.set_defaults(**{__POINTS: True})


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(GENERATE_GLOBAL_IDS,
                              help="Adds globals ids for points and cells.")
    fill_generate_global_ids_subparser(p)
    vtk_output_parsing.fill_vtk_output_subparser(p)


def display_results(options: Options, result: Result):
    logging.info(result.info)
