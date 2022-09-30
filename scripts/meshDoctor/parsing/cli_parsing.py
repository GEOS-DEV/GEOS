import argparse
from collections import OrderedDict
from dataclasses import dataclass
import logging
import textwrap
from typing import List, Dict, Set

from . import all_checks_helpers


__OPTIONS_SEP = ":"
__KV_SEP = "="

__VERBOSE_KEY = "verbose"


@dataclass(frozen=True)
class Arguments:
    vtk_input_file: str
    checks: OrderedDict()


def __get_checks_help_msg() -> str:
    """
    Gathers all the doc messages into one string.
    :return: A string.
    """
    tmp = []
    for check_name, check_helper in all_checks_helpers.items():
        h = check_name + ": " + check_helper.get_help()
        tmp.append(h)
    return "\n".join(tmp)


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
            logging.warning(f"Keys \"{', '.join(invalid_keys)}\" are not valid options of \"{check_name}\". Ignoring them.")


def parse_and_set_verbosity(cli_args: List[str]) -> None:
    """
    Parse the verbosity flag only. And sets the logger's level accordingly.
    :param cli_args: The list of arguments (as strings)
    :return: None
    """
    dummy_verbosity_parser = argparse.ArgumentParser(add_help=None)
    dummy_verbosity_parser.add_argument('-v',
                                        '--verbose',
                                        action='count',
                                        default=1,
                                        dest=__VERBOSE_KEY)
    args = dummy_verbosity_parser.parse_known_args(cli_args[1:])[0]
    v = vars(args)[__VERBOSE_KEY]
    verbosity = 50 - (10 * v) if v > 0 else 0
    logging.getLogger().setLevel(verbosity)
    logging.info(f"Logger level set to \"{logging.getLevelName(verbosity)}\"")


def parse(cli_args: List[str]) -> Arguments:
    """
    Parse the command line arguments and return the corresponding structure.
    :param cli_args: The list of arguments (as strings)
    :return: The struct
    """
    vtk_input_file_key = "vtk_input_file"

    epilog_msg = """\
    Note that checks are dynamically loaded.
    An option may be missing because of an unloaded module.
    Increase verbosity to get full information.
    """

    parser = argparse.ArgumentParser(description='Inspects meshes for GEOSX.',
                                     epilog=textwrap.dedent(epilog_msg),
                                     formatter_class=argparse.RawTextHelpFormatter)
    misc_grp = parser.add_argument_group('main')
    # Nothing will be done with this verbosity input.
    # It's only here for the `--help` message.
    # `parse_verbosity` does the real parsing instead.
    misc_grp.add_argument('-v',
                          '--verbose',
                          action='count',
                          default=1,
                          dest=__VERBOSE_KEY,
                          help="Use -v for warning, -vv for info, -vvv for debug. Defaults to 'ERROR'.")
    misc_grp.add_argument('-i',
                          '--vtk-input-file',
                          metavar='VTK_MESH_FILE',
                          type=str,
                          required=True,
                          dest=vtk_input_file_key)
    checks_grp = parser.add_argument_group('checks')
    for check_name, check_helper in all_checks_helpers.items():
        checks_grp.add_argument('--' + check_name,
                                metavar='opt0' + __KV_SEP + 'val0' + __OPTIONS_SEP + 'opt1' + __KV_SEP + 'val1...',
                                type=check_helper.parse_cli_options,
                                dest=check_name,
                                required=False,
                                help=check_helper.get_help())
    args = parser.parse_args(cli_args[1:])

    args = vars(args)  # converting to `dict` allows to access keys by variable.
    checks = OrderedDict()
    for check_name in all_checks_helpers.keys():
        options = args[check_name]
        checks[check_name] = options
        logging.debug(f"Check \"{check_name}\" was found with options \"{options}\".")

    return Arguments(vtk_input_file=args[vtk_input_file_key], checks=checks)
