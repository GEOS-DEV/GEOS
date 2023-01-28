import logging
import textwrap

from checks.generate_global_ids import Options, Result

from . import cli_parsing, vtk_output_parsing, GENERATE_GLOBAL_IDS

__ALL_KEYWORDS = {
    *vtk_output_parsing.get_vtk_output_keywords(),
}


def get_help():
    msg = f"""\
    Adds globals ids for points and cells.
    
    {vtk_output_parsing.get_vtk_output_help()}
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for collocated nodes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(GENERATE_GLOBAL_IDS, __ALL_KEYWORDS, options)
    vtk_output = vtk_output_parsing.parse_cli_options(options)
    return Options(vtk_output=vtk_output)


def display_results(options: Options, result: Result):
    logging.info(result.info)
