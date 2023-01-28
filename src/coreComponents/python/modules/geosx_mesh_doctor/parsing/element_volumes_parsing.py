import logging
import textwrap
from typing import Dict

from checks.element_volumes import Options, Result

from . import cli_parsing, ELEMENT_VOLUMES

__MIN_VOLUME = "min"
__MIN_VOLUME_DEFAULT = 0.


def get_help():
    msg = f"""\
    Checks if the volumes of the elements are greater than {__MIN_VOLUME}.
    
    {__MIN_VOLUME} [float]: The minimum acceptable volume. Defaults to {__MIN_VOLUME_DEFAULT}.
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for elements volumes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options: Dict[str, str] = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(ELEMENT_VOLUMES, {__MIN_VOLUME}, options)
    return Options(float(options.get(__MIN_VOLUME, __MIN_VOLUME_DEFAULT)))


def display_results(options: Options, result: Result):
    logging.error(f"You have {len(result.element_volumes)} elements with volumes smaller than {options.min_volume}.")
    if result.element_volumes:
        logging.error("The elements indices and their volumes are:\n" + "\n".join(map(str, result.element_volumes)))
