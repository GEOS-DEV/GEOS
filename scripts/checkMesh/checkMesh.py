import logging
import sys

from checks import all_checks
from parsing import all_checks_helpers
from parsing.cli_parsing import parse


def main():
    args = parse(sys.argv)
    logging.info(f"Checking mesh {args.vtk_input_file}.")
    for check_name, check_options in args.checks.items():
        try:
            check = all_checks[check_name]
        except KeyError as e:
            logging.critical(f"Check {check_name} is not a valid check.")
            sys.exit(1)
        helper = all_checks_helpers[check_name]
        options = helper.parse(check_options)
        result = check(args.vtk_input_file, options)
        helper.present(options, result)


if __name__ == '__main__':
    main()
