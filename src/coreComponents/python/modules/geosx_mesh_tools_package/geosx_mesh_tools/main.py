import argparse
import logging
import sys


def build_abaqus_converter_input_parser() -> argparse.ArgumentParser:
    """Build the input argument parser

    Returns:
        argparse.ArgumentParser: a parser instance
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input abaqus mesh file name')
    parser.add_argument('output', type=str, help='Output gmsh/vtu mesh file name')
    parser.add_argument('-v', '--verbose', help='Increase verbosity level', action="store_true")
    return parser


def main() -> None:
    """
    Entry point for the abaqus convertor console script

    Args:
        input (str): Input abaqus mesh file name
        output (str): Output mesh file name
        -v/--verbose (flag): Increase verbosity level
    """
    from geosx_mesh_tools import abaqus_converter

    # Parse the user arguments
    parser = build_abaqus_converter_input_parser()
    args = parser.parse_args()

    # Set up a logger
    logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger(__name__)
    if args.verbose:
        logger.setLevel(logging.INFO)

    # Call the converter
    err = 0
    if ('.msh' in args.output):
        err = abaqus_converter.convert_abaqus_to_gmsh(args.input, args.output, logger)
    else:
        err = abaqus_converter.convert_abaqus_to_vtu(args.input, args.output, logger)
    if err:
        sys.exit('Warnings detected: check the output file for potential errors!')


if __name__ == '__main__':
    main()
