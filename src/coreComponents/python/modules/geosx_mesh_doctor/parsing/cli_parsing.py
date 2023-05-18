import argparse
import logging
import textwrap
from typing import List

__VERBOSE_KEY = "verbose"
__QUIET_KEY = "quiet"

__VERBOSITY_FLAG = "v"
__QUIET_FLAG = "q"


def parse_and_set_verbosity(cli_args: List[str]) -> None:
    """
    Parse the verbosity flag only. And sets the logger's level accordingly.
    :param cli_args: The list of arguments (as strings)
    :return: None
    """
    dummy_verbosity_parser = argparse.ArgumentParser(add_help=False)
    dummy_verbosity_parser.add_argument('-' + __VERBOSITY_FLAG,
                                        '--' + __VERBOSE_KEY,
                                        action='count',
                                        default=2,
                                        dest=__VERBOSE_KEY)
    dummy_verbosity_parser.add_argument('-' + __QUIET_FLAG,
                                        '--' + __QUIET_KEY,
                                        action='count',
                                        default=0,
                                        dest=__QUIET_KEY)
    args = dummy_verbosity_parser.parse_known_args(cli_args[1:])[0]
    d = vars(args)
    v = d[__VERBOSE_KEY] - d[__QUIET_KEY]
    verbosity = logging.CRITICAL - (10 * v)
    if verbosity < logging.DEBUG:
        verbosity = logging.DEBUG
    elif verbosity > logging.CRITICAL:
        verbosity = logging.CRITICAL
    logging.getLogger().setLevel(verbosity)
    logging.info(f"Logger level set to \"{logging.getLevelName(verbosity)}\"")


def init_parser() -> argparse.ArgumentParser:
    vtk_input_file_key = "vtk_input_file"

    epilog_msg = f"""\
        Note that checks are dynamically loaded.
        An option may be missing because of an unloaded module.
        Increase verbosity (-{__VERBOSITY_FLAG}, -{__VERBOSITY_FLAG * 2}) to get full information.
        """
    formatter = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=8)
    parser = argparse.ArgumentParser(description='Inspects meshes for GEOSX.',
                                     epilog=textwrap.dedent(epilog_msg),
                                     formatter_class=formatter)
    # Nothing will be done with this verbosity/quiet input.
    # It's only here for the `--help` message.
    # `parse_verbosity` does the real parsing instead.
    parser.add_argument('-' + __VERBOSITY_FLAG,
                        action='count',
                        default=2,
                        dest=__VERBOSE_KEY,
                        help=f"Use -{__VERBOSITY_FLAG} 'INFO', -{__VERBOSITY_FLAG * 2} for 'DEBUG'. Defaults to 'WARNING'.")
    parser.add_argument('-' + __QUIET_FLAG,
                        action='count',
                        default=0,
                        dest=__QUIET_KEY,
                        help=f"Use -{__QUIET_FLAG} to reduce the verbosity of the output.")
    parser.add_argument('-i',
                        '--vtk-input-file',
                        metavar='VTK_MESH_FILE',
                        type=str,
                        required=True,
                        dest=vtk_input_file_key)
    return parser
