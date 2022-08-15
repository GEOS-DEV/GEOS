"""Command line tools for geosx_xml_tools"""

import sys
import argparse
import os
import time
from geosx_xml_tools import xml_processor


def parse_arguments():
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
    return parser.parse_known_args()


def check_mpi_rank():
    """Check the MPI rank

    Returns:
        int: MPI rank
    """
    rank = 0
    mpi_rank_key_options = ['OMPI_COMM_WORLD_RANK', 'PMI_RANK']
    for k in mpi_rank_key_options:
        if k in os.environ:
            rank = int(os.environ[k])
    return rank


def wait_for_file_write_rank_0(target_file_argument=0, max_wait_time=100, max_startup_delay=1):
    """Constructor for a function decorator that waits for a target file to be written on rank 0

    Args:
        target_file_argument (int): Index of the filename argument in the decorated function
        max_wait_time (float): Maximum amount of time to wait (seconds)
        max_startup_delay (float): Maximum delay allowed for thread startup (seconds)

    Returns:
        Wrapped function
    """

    def wait_for_file_write_rank_0_inner(writer):
        """Intermediate constructor for the function decorator

        Args:
            writer (typing.Callable): A function that writes to a file
        """

        def wait_for_file_write_rank_0_decorator(*args, **kwargs):
            """Apply the writer on rank 0, and wait for completion on other ranks
            """
            # Check the target file status
            rank = check_mpi_rank()
            fname = ''
            if isinstance(target_file_argument, int):
                fname = args[target_file_argument]
            else:
                fname = kwargs[target_file_argument]

            target_file_exists = os.path.isfile(fname)
            target_file_edit_time = 0
            if target_file_exists:
                target_file_edit_time = os.path.getmtime(fname)

                # Variations in thread startup times may mean the file has already been processed
                # If the last edit was done within the specified time, then allow the thread to proceed
                if (abs(target_file_edit_time - time.time()) < max_startup_delay):
                    target_file_edit_time = 0

            # Go into the target process or wait for the expected file update
            if (rank == 0):
                return writer(*args, **kwargs)
            else:
                ta = time.time()
                while (time.time() - ta < max_wait_time):
                    if target_file_exists:
                        if (os.path.getmtime(fname) > target_file_edit_time):
                            break
                    else:
                        if os.path.isfile(fname):
                            break
                    time.sleep(0.1)

        return wait_for_file_write_rank_0_decorator

    return wait_for_file_write_rank_0_inner


def preprocess_serial():
    """
    Entry point for the geosx_xml_tools console script
    """
    # Process the xml file
    args, unknown_args = parse_arguments()

    # Attempt to only process the file on rank 0
    # Note: The rank here is determined by inspecting the system environment variables
    #       While this is not the preferred way of doing so, it avoids mpi environment errors
    #       If the rank detection fails, then it will preprocess the file on all ranks, which
    #       sometimes cause a (seemingly harmless) file write conflict.
    # processor = xml_processor.process
    processor = wait_for_file_write_rank_0(target_file_argument='outputFile', max_wait_time=100)(xml_processor.process)

    compiled_name = processor(args.input,
                              outputFile=args.compiled_name,
                              schema=args.schema,
                              verbose=args.verbose,
                              parameter_override=args.parameters)
    if not compiled_name:
        if args.compiled_name:
            compiled_name = args.compiled_name
        else:
            raise Exception(
                'When applying the preprocessor in parallel (outside of pygeosx), the --compiled_name argument is required'
            )

    # Note: the return value may be passed to sys.exit, and cause bash to report an error
    # return format_geosx_arguments(compiled_name, unknown_args)
    print(compiled_name)


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
    """Format GEOSX arguments

    Args:
        compiled_name (str): Name of the compiled xml file
        unknown_args (list): List of unprocessed arguments

    Returns:
        list: List of arguments to pass to GEOSX
    """
    geosx_args = [sys.argv[0], '-i', compiled_name]
    if unknown_args:
        geosx_args.extend(unknown_args)

    # Print the output name for use in bash scripts
    print(compiled_name)
    return geosx_args


if __name__ == "__main__":
    preprocess_serial()
