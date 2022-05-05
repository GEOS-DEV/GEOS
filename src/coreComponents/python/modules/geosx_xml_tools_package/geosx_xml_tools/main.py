"""Command line tools for geosx_xml_tools"""

import sys
import argparse
from geosx_xml_tools import xml_processor


def parse_arguments():
    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, action='append', help='Input file name (multiple allowed)')
    parser.add_argument('-c', '--compiled-name', type=str, help='Compiled xml file name (otherwise, it is randomly genrated)', default='')
    parser.add_argument('-s', '--schema', type=str, help='GEOSX schema to use for validation', default='')
    parser.add_argument('-v', '--verbose', type=int, help='Verbosity of outputs', default=0)
    parser.add_argument('-p', '--parameters', nargs='+', action='append', help='Parameter overrides (name value, multiple allowed)', default=[])
    return parser.parse_known_args()


def preprocess_serial():
    """
    Entry point for the geosx_xml_tools console script
    """
    # Process the xml file
    args, unknown_args = parse_arguments()
    compiled_name = xml_processor.process(args.input,
                                          outputFile=args.compiled_name,
                                          schema=args.schema,
                                          verbose=args.verbose,
                                          parameter_override=args.parameters)
    return format_geosx_arguments(compiled_name, unknown_args)


def preprocess_parallel():
    """
    MPI aware xml preprocesing
    """
    # Process the xml file
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    args, unknown_args = parse_arguments()
    compiled_name = ''
    if (rank == 0):
        compiled_name = xml_processor.process(args.input,
                                              outputFile=args.compiled_name,
                                              schema=args.schema,
                                              verbose=args.verbose,
                                              parameter_override=args.parameters)
    compiled_name = comm.bcast(compiled_name, root=0)
    return format_geosx_arguments(compiled_name, unknown_args)


def format_geosx_arguments(compiled_name, unknown_args):
    # Return GEOSX arguments
    geosx_args = [sys.argv[0], '-i', compiled_name]
    if unknown_args:
        geosx_args.extend(unknown_args)

    # Print the output name for use in bash scripts
    print(compiled_name, flush=True)
    return geosx_args


if __name__ == "__main__":
    preprocess_serial()
