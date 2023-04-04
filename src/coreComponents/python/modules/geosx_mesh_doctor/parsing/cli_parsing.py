import argparse
from collections import OrderedDict
from dataclasses import dataclass
import logging
from typing import List, Dict, Set


__OPTIONS_SEP = ":"
__KV_SEP = "="

__VERBOSE_KEY = "verbose"
__QUIET_KEY = "quiet"


@dataclass(frozen=True)
class Arguments:
    vtk_input_file: str
    checks: OrderedDict()


def parse_cli_option(s: str) -> Dict[str, str]:
    """
    Taking a CLI option string, and converts into a dict.
    :param s: A string of options like "foo=bar:fizz=buzz".
    :return: The corresponding {foo:bar, fizz:buzz} dict.
    """
    result = OrderedDict()
    for opt in s.split(__OPTIONS_SEP):
        k, v = opt.split(__KV_SEP)
        result[k] = v
    return result


def validate_cli_options(check_name: str, valid_keys: Set[str], options: Dict[str, str]):
    """
    Checks the user input options and logs if an options is not recognized.
    :param check_name: The key/name of the check.
    :param valid_keys: The option keys that are valid for the considered check.
    :param options: The user options.
    :return: None
    :todo: Deal with default options.
    """
    invalid_keys = set(options.keys()) - valid_keys
    if invalid_keys:
        if len(invalid_keys) == 1:
            logging.warning(f"Key \"{invalid_keys.pop()}\" is not a valid option of \"{check_name}\". Ignoring it.")
        else:
            logging.warning(
                f"Keys \"{', '.join(invalid_keys)}\" are not valid options of \"{check_name}\". Ignoring them.")


def parse_and_set_verbosity(cli_args: List[str]) -> None:
    """
    Parse the verbosity flag only. And sets the logger's level accordingly.
    :param cli_args: The list of arguments (as strings)
    :return: None
    """
    dummy_verbosity_parser = argparse.ArgumentParser(add_help=None)
    dummy_verbosity_parser.add_argument('-v', '--verbose', action='count', default=2, dest=__VERBOSE_KEY)
    dummy_verbosity_parser.add_argument('-q', '--quiet', action='count', default=0, dest=__QUIET_KEY)
    args = dummy_verbosity_parser.parse_known_args(cli_args[1:])[0]
    d = vars(args)
    v = d[__VERBOSE_KEY] - d[__QUIET_KEY]
    verbosity = logging.CRITICAL - (10 * v)
    if verbosity < logging.DEBUG:
        verbosity = logging.DEBUG
    elif verbosity > logging.CRITICAL:
        verbosity = logging.CRITICAL
    logging.getLogger().setLevel(verbosity)
    logging.info(f"Logger level set to \"{logging.getLevelName(verbosity)}\"")
