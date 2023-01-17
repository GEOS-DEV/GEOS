import logging
import multiprocessing
import textwrap

from checks.supported_elements import Options, Result

from . import cli_parsing, SUPPORTED_ELEMENTS

__CHUNK_SIZE = "chunck_size"
__NUM_PROC = "num_proc"


__ALL_KEYWORDS = {__CHUNK_SIZE, __NUM_PROC}

__CHUNK_SIZE_DEFAULT = 1
__NUM_PROC_DEFAULT = multiprocessing.cpu_count()


def get_help():
    msg = f"""\
    Check that all the elements of the mesh are supported by GEOSX.
    
    {__CHUNK_SIZE} [int]: Defaults chunk size for parallel processing to {__CHUNK_SIZE_DEFAULT}.
    {__NUM_PROC} [int]: Number of threads used for parallel processinf. Defaults to {__NUM_PROC_DEFAULT}/
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for collocated nodes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(SUPPORTED_ELEMENTS, __ALL_KEYWORDS, options)
    chunk_size = int(options.get(__CHUNK_SIZE, __CHUNK_SIZE_DEFAULT))
    num_proc = int(options.get(__NUM_PROC, __NUM_PROC_DEFAULT))
    return Options(chunk_size=chunk_size, num_proc=num_proc)


def display_results(options: Options, result: Result):
    if result.unsupported_polyhedron_elements:
        logging.error(f"There is/are {len(result.unsupported_polyhedron_elements)} polyhedra that may not be converted to supported elements.")
        logging.error(f"The list unsupported polyhedra is {tuple(sorted(result.unsupported_polyhedron_elements))}.")
    if result.unsupported_std_elements_types:
        logging.error(f"The list of the unsupported types is {tuple(sorted(result.unsupported_std_elements_types))}.")
