import logging
import multiprocessing

from checks.supported_elements import Options, Result

from . import SUPPORTED_ELEMENTS

__CHUNK_SIZE = "chunck_size"
__NUM_PROC = "nproc"


__ALL_KEYWORDS = {__CHUNK_SIZE, __NUM_PROC}

__CHUNK_SIZE_DEFAULT = 1
__NUM_PROC_DEFAULT = multiprocessing.cpu_count()


def convert(parsed_options) -> Options:
    return Options(chunk_size=parsed_options[__CHUNK_SIZE],
                   num_proc=parsed_options[__NUM_PROC])


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(SUPPORTED_ELEMENTS,
                              help="Check that all the elements of the mesh are supported by GEOSX.")
    p.add_argument('--' + __CHUNK_SIZE,
                   type=int,
                   required=False,
                   metavar=__CHUNK_SIZE_DEFAULT,
                   default=__CHUNK_SIZE_DEFAULT,
                   help=f"[int]: Defaults chunk size for parallel processing to {__CHUNK_SIZE_DEFAULT}")
    p.add_argument('--' + __NUM_PROC,
                   type=int,
                   required=False,
                   metavar=__NUM_PROC_DEFAULT,
                   default=__NUM_PROC_DEFAULT,
                   help=f"[int]: Number of threads used for parallel processing. Defaults to your CPU count {__NUM_PROC_DEFAULT}.")


def display_results(options: Options, result: Result):
    if result.unsupported_polyhedron_elements:
        logging.error(f"There is/are {len(result.unsupported_polyhedron_elements)} polyhedra that may not be converted to supported elements.")
        logging.error(f"The list of the unsupported polyhedra is\n{tuple(sorted(result.unsupported_polyhedron_elements))}.")
    else:
        logging.info("All the polyhedra (if any) can be converted to supported elements.")
    if result.unsupported_std_elements_types:
        logging.error(f"There are unsupported vtk standard element types. The list of those vtk types is {tuple(sorted(result.unsupported_std_elements_types))}.")
    else:
        logging.info("All the standard vtk element types (if any) are supported.")