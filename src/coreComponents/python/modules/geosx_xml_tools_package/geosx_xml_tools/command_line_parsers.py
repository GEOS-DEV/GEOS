
import argparse
from typing import Tuple, Iterable


def build_preprocessor_input_parser() -> argparse.ArgumentParser:
    """Build the argument parser

    Returns:
        argparse.ArgumentParser: The parser
    """
    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, action='append', help='Input file name (multiple allowed)')
    parser.add_argument('-c',
                        '--compiled-name',
                        type=str,
                        help='Compiled xml file name (otherwise, it is randomly genrated)',
                        default='')
    parser.add_argument('-s', '--schema', type=str, help='GEOSX schema to use for validation', default='')
    parser.add_argument('-v', '--verbose', type=int, help='Verbosity of outputs', default=0)
    parser.add_argument('-p',
                        '--parameters',
                        nargs='+',
                        action='append',
                        help='Parameter overrides (name value, multiple allowed)',
                        default=[])

    return parser


def parse_xml_preprocessor_arguments() -> Tuple[argparse.Namespace, Iterable[str]]:
    """Parse user arguments

    Args:
        -i/--input (str): Input file name (multiple allowed)
        -c/--compiled-name (str): Compiled xml file name
        -s/--schema (str): Path to schema to use for validation
        -v/--verbose (int): Verbosity of outputs
        -p/--parameters (str): Parameter overrides (name and value, multiple allowed)

    Returns:
        list: The remaining unparsed argument strings
    """
    parser = build_preprocessor_input_parser()
    return parser.parse_known_args()


def build_xml_formatter_input_parser() -> argparse.ArgumentParser:
    """Build the argument parser

    Returns:
        argparse.ArgumentParser: the parser instance
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file name')
    parser.add_argument('-i', '--indent', type=int, help='Indent size', default=2)
    parser.add_argument('-s', '--style', type=int, help='Indent style', default=0)
    parser.add_argument('-d', '--depth', type=int, help='Block separation depth', default=2)
    parser.add_argument('-a', '--alphebitize', type=int, help='Alphebetize attributes', default=0)
    parser.add_argument('-c', '--close', type=int, help='Close tag style', default=0)
    parser.add_argument('-n', '--namespace', type=int, help='Include namespace', default=0)
    return parser


def build_attribute_coverage_input_parser() -> argparse.ArgumentParser:
    """Build attribute coverage redundancy input parser

    Returns:
        argparse.ArgumentParser: parser instance
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--root', type=str, help='GEOSX root', default='')
    parser.add_argument('-o', '--output', type=str, help='Output file name', default='attribute_test.xml')
    return parser


def build_xml_redundancy_input_parser() -> argparse.ArgumentParser:
    """Build xml redundancy input parser

    Returns:
        argparse.ArgumentParser: parser instance
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--root', type=str, help='GEOSX root', default='')
    return parser


