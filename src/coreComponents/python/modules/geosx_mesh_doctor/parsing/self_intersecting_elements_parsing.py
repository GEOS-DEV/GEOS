import logging

import numpy

from checks.self_intersecting_elements import Options, Result

from . import SELF_INTERSECTING_ELEMENTS

__TOLERANCE = "min"
__TOLERANCE_DEFAULT = numpy.finfo(float).eps


def convert(parsed_options) -> Options:
    tolerance = parsed_options[__TOLERANCE]
    if tolerance == 0:
        logging.warning("Having tolerance set to 0 can induce lots of false positive results (adjacent faces may be considered intersecting).")
    elif tolerance < 0:
        raise ValueError(f"Negative tolerance ({tolerance}) in the {SELF_INTERSECTING_ELEMENTS} check is not allowed.")
    return Options(tolerance=tolerance)


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(SELF_INTERSECTING_ELEMENTS,
                              help="Checks if the faces of the elements are self intersecting.")
    p.add_argument('--' + __TOLERANCE,
                   type=float,
                   required=False,
                   metavar=__TOLERANCE_DEFAULT,
                   default=__TOLERANCE_DEFAULT,
                   help=f"[float]: The tolerance in the computation. Defaults to your machine precision {__TOLERANCE_DEFAULT}.")


def display_results(options: Options, result: Result):
    logging.error(f"You have {len(result.intersecting_faces_elements)} elements with self intersecting faces.")
    if result.intersecting_faces_elements:
        logging.error("The elements indices are:\n" + ", ".join(map(str, result.intersecting_faces_elements)))
