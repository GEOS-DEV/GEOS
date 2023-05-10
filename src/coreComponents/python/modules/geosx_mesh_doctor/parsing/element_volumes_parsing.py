import logging

from checks.element_volumes import Options, Result

from . import ELEMENT_VOLUMES

__MIN_VOLUME = "min"
__MIN_VOLUME_DEFAULT = 0.


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(ELEMENT_VOLUMES,
                              help=f"Checks if the volumes of the elements are greater than \"{__MIN_VOLUME}\".")
    p.add_argument('--' + __MIN_VOLUME,
                   type=float,
                   metavar=__MIN_VOLUME_DEFAULT,
                   default=__MIN_VOLUME_DEFAULT,
                   required=True,
                   help=f"[float]: The minimum acceptable volume. Defaults to {__MIN_VOLUME_DEFAULT}.")


def convert(parsed_options) -> Options:
    """
    From the parsed cli options, return the converted options for elements volumes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    return Options(min_volume=parsed_options[__MIN_VOLUME])


def display_results(options: Options, result: Result):
    logging.error(f"You have {len(result.element_volumes)} elements with volumes smaller than {options.min_volume}.")
    if result.element_volumes:
        logging.error("The elements indices and their volumes are:\n" + "\n".join(map(str, result.element_volumes)))
