import sys

try:
    min_python_version = (3, 7)
    assert sys.version_info >= min_python_version
except AssertionError as e:
    print(f"Please update python to at least version {'.'.join(map(str, min_python_version))}.")
    sys.exit(1)

import logging

from parsing import CheckHelper
from parsing.cli_parsing import parse, parse_and_set_verbosity
import register


def main():
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s')
    parse_and_set_verbosity(sys.argv)
    all_checks, all_checks_helpers = register.register()
    args = parse(all_checks_helpers, sys.argv)
    logging.info(f"Checking mesh \"{args.vtk_input_file}\".")
    for check_name, check_options in args.checks.items():
        # If there is no option, this means that the check was not requested by the user
        if not check_options:
            continue
        try:
            check = all_checks[check_name]
        except KeyError as e:
            logging.critical(f"Check {check_name} is not a valid check.")
            sys.exit(1)
        helper: CheckHelper = all_checks_helpers[check_name]
        result = check(args.vtk_input_file, check_options)
        helper.display_results(check_options, result)


if __name__ == '__main__':
    main()
