import argparse
from collections import OrderedDict
from dataclasses import dataclass
import logging
import sys
from typing import List, Dict, Set

from . import all_checks_helpers


__OPTIONS_SEP = ":"
__KV_SEP = "="


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


def __parse_cli_options(s: List[str]) -> Dict[str, str]:
    """
    Parse a command line option into a dict:
    :param s: Ths cli string, e.g. "a=1;b=5"
    :return: The user options as a dict.
    """
    result = OrderedDict()
    for option in s:
        k, v = option.split(__KV_SEP)
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
    logging.warning(f"Key(s) [{', '.join(invalid_keys)}] is (are) not valid option(s) of {check_name}. Ignoring.")


def parse(cli_args: List[str]) -> Arguments:
    """
    Parse the command line arguments and return the corresponding structure.
    :param cli_args: The list of arguments (as strings)
    :return: The struct
    """
    parser = argparse.ArgumentParser(description='Inspects meshes for GEOSX.')
    parser.add_argument('-i',
                        '--vtk-input-file',
                        metavar='VTK_MESH_FILE',
                        type=str,
                        required=True,
                        dest="vtk_input_file")
    parser.add_argument('-c',
                        '--check',
                        metavar='CHECK' + __OPTIONS_SEP + 'opt0' + __KV_SEP + 'val0' + __OPTIONS_SEP + 'opt1' + __KV_SEP + 'val1...',
                        type=str,
                        action='extend',
                        nargs='+',
                        dest="checks",
                        required=True,
                        help=__get_checks_help_msg())
    parser.add_argument('-v',
                        '--verbose',
                        action='count',
                        default=1,
                        help="Use -v for warning, -vv for info, -vvv for debug.")
    args = parser.parse_args(cli_args[1:])

    verbosity = 50 - (10 * args.verbose) if args.verbose > 0 else 0
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s', level=verbosity)

    checks = OrderedDict()
    for c in args.checks:
        split_check_input = c.split(__OPTIONS_SEP)
        check_name, check_parameters = split_check_input[0], split_check_input[1:]
        if not check_name:
            logging.critical("Check has no name, aborting")
            sys.exit(1)
        checks[check_name] = __parse_cli_options(check_parameters)

    return Arguments(vtk_input_file=args.vtk_input_file, checks=checks)
