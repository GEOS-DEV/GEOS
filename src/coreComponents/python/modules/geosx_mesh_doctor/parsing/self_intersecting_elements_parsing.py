import logging
import textwrap
from typing import Dict

import numpy

from checks.self_intersecting_elements import Options, Result

from . import cli_parsing, SELF_INTERSECTING_ELEMENTS

__TOLERANCE = "min"
__TOLERANCE_DEFAULT = numpy.finfo(float).eps


def get_help():
    msg = f"""\
    Checks if the faces of the elements are self intersecting.
    
    {__TOLERANCE} [float]: The tolerance in the computation. Defaults to machine precision {__TOLERANCE_DEFAULT}.
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for self intersecting elements check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options: Dict[str, str] = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(SELF_INTERSECTING_ELEMENTS, {__TOLERANCE}, options)
    tolerance = float(options.get(__TOLERANCE, __TOLERANCE_DEFAULT))
    if tolerance == 0:
        logging.warning("Having tolerance set to 0 can induce lots of false positive results (adjacent faces may be considered intersecting).")
    elif tolerance < 0:
        raise ValueError(f"Negative tolerance ({tolerance}) in the {SELF_INTERSECTING_ELEMENTS} check is not allowed.")
    return Options(tolerance=tolerance)


def display_results(options: Options, result: Result):
    logging.error(f"You have {len(result.intersecting_faces_elements)} elements with self intersecting faces.")
    if result.intersecting_faces_elements:
        logging.error("The elements indices are:\n" + ", ".join(map(str, result.intersecting_faces_elements)))
