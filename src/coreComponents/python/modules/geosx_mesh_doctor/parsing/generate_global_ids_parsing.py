import logging

from checks.generate_global_ids import Options, Result

from . import vtk_output_parsing, GENERATE_GLOBAL_IDS


def convert(parsed_options) -> Options:
    return Options(vtk_output=vtk_output_parsing.convert(parsed_options))


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(GENERATE_GLOBAL_IDS,
                              help="Adds globals ids for points and cells.")
    vtk_output_parsing.fill_vtk_output_subparser(p)


def display_results(options: Options, result: Result):
    logging.info(result.info)
