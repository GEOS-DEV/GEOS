import sys

try:
    min_python_version = (3, 7)
    assert sys.version_info >= min_python_version
except AssertionError as e:
    print(f"Please update python to at least version {'.'.join(map(str, min_python_version))}.")
    sys.exit(1)

import logging

from parsing import CheckHelper
from parsing.cli_parsing import parse_and_set_verbosity
import register


def main():
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s')
    parse_and_set_verbosity(sys.argv)
    main_parser, all_checks, all_checks_helpers = register.register()
    args = main_parser.parse_args(sys.argv[1:])
    logging.info(f"Checking mesh \"{args.vtk_input_file}\".")
    check_options = all_checks_helpers[args.subparsers].convert(vars(args))
    try:
        check = all_checks[args.subparsers]
    except KeyError as e:
        logging.critical(f"Check {args.subparsers} is not a valid check.")
        sys.exit(1)
    helper: CheckHelper = all_checks_helpers[args.subparsers]
    result = check(args.vtk_input_file, check_options)
    helper.display_results(check_options, result)


if __name__ == '__main__':
    main()
