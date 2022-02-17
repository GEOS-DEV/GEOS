from mpi4py import MPI
from geosx_xml_tools import xml_processor
import sys
import argparse


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def apply_xml_preprocessor():
    """
    Applies the xml preprocessor to the input file
    before handing it to GEOSX, and modifies the input
    arguments to point to the new file
    """
    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, action='append', help='Input file name')
    parser.add_argument('-f', '--file', type=str, help='Output file name (otherwise, it is randomly genrated)', default='')
    parser.add_argument('-s', '--schema', type=str, help='GEOSX schema to use for validation', default='')
    parser.add_argument('-v', '--verbose', type=int, help='Verbosity of outputs', default=0)
    parser.add_argument('-p', '--parameters', nargs='+', action='append', help='Parameter overrides', default=[])
    args, unknown_args = parser.parse_known_args()

    new_input_file = ''
    if (rank == 0):
        print('Applying xml preprocessor...', flush=True)
        new_input_file = xml_processor.process(args.input,
                                               outputFile=args.file,
                                               schema=args.schema,
                                               verbose=args.verbose,
                                               parameter_override=args.parameters)
        print('  the compiled filename is: %s' % (new_input_file), flush=True)

    # Broadcast and set the new argument list
    new_input_file = comm.bcast(new_input_file, root=0)
    new_args = [sys.argv[0], '-i', new_input_file]
    if unknown_args:
        new_args.extend(unknown_args)
    sys.argv = new_args
