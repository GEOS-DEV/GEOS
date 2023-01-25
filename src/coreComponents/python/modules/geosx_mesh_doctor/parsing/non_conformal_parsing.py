import logging
import textwrap

from checks.non_conformal import Options, Result

from . import cli_parsing, NON_CONFORMAL

__ANGLE_TOLERANCE = "angle_tolerance"
__POINT_TOLERANCE = "point_tolerance"
__FACE_TOLERANCE = "face_tolerance"

__ANGLE_TOLERANCE_DEFAULT = 10.

__ALL_KEYWORDS = {__ANGLE_TOLERANCE, __POINT_TOLERANCE, __FACE_TOLERANCE}


def get_help():
    msg = f"""\
    Detects non conformal elements. [EXPERIMENTAL]
    
    {__ANGLE_TOLERANCE} [float]: angle tolerance in degrees. Defaults to {__ANGLE_TOLERANCE_DEFAULT}.
    {__POINT_TOLERANCE} [float]: tolerance for two points to be considered collocated.
    {__FACE_TOLERANCE} [float]: tolerance for two faces to be considered "touching".
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for collocated nodes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(NON_CONFORMAL, __ALL_KEYWORDS, options)
    angle_tolerance = float(options.get(__ANGLE_TOLERANCE, __ANGLE_TOLERANCE_DEFAULT))
    point_tolerance = float(options[__POINT_TOLERANCE])
    face_tolerance = float(options[__FACE_TOLERANCE])
    return Options(angle_tolerance=angle_tolerance,
                   point_tolerance=point_tolerance,
                   face_tolerance=face_tolerance)


def display_results(options: Options, result: Result):
    non_conformal_cells = []
    for i, j in result.non_conformal_cells:
        non_conformal_cells += i, j
    non_conformal_cells = set(non_conformal_cells)
    logging.error(f"You have {len(non_conformal_cells)} non conformal cells.\n{', '.join(map(str, sorted(non_conformal_cells)))}")
