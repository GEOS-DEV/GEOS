import logging
import sys

from checks import all_checks
from parsing import register
from parsing import all_checks_helpers
from parsing.cli_parsing import parse, parse_and_set_verbosity


def main():
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s')
    parse_and_set_verbosity(sys.argv)
    register.register()
    args = parse(sys.argv)
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
        helper = all_checks_helpers[check_name]
        result = check(args.vtk_input_file, check_options)
        helper.display_results(check_options, result)


if __name__ == '__main__':
    main()
