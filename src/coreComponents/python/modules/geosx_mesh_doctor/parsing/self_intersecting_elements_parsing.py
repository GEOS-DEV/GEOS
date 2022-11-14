import logging
import textwrap
from typing import Dict

from checks.self_intersecting_elements import Options, Result

from . import cli_parsing, SELF_INTERSECTING_ELEMENTS

__TOLERANCE = "min"
__TOLERANCE_DEFAULT = 0.


def get_help():
    msg = f"""\
    Checks if the faces of the elements are self intersecting.
        {__TOLERANCE} [float]: The tolerance in the computation. Defaults to {__TOLERANCE_DEFAULT}.
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
    return Options(float(options.get(__TOLERANCE, __TOLERANCE_DEFAULT)))


def display_results(options: Options, result: Result):
    logging.error(f"You have {len(result.jumbled_elements)} elements with self intersecting faces.")
    if result.jumbled_elements:
        logging.error("The elements indices are:\n" + "\n".join(map(str, result.jumbled_elements)))
