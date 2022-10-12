import logging
import textwrap

from checks import GENERATE_GLOBAL_IDS
from checks.generate_global_ids import Options, Result

from . import cli_parsing


__OUTPUT_FILE = "output"

__ALL_KEYWORDS = {__OUTPUT_FILE}


def get_help():
    msg = f"""\
    Adds globals ids for points and cells.
    
    {__OUTPUT_FILE} [string]: The vtk output destination.
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
    output = options.get(__OUTPUT_FILE, "")
    return Options(output=output)


def display_results(options: Options, result: Result):
    logging.info(result.info)
