import os.path
import logging
import textwrap
from typing import Dict

from checks.vtk_utils import VtkOutput


__OUTPUT_FILE = "output"
__OUTPUT_BINARY_MODE = "data_mode"
__OUTPUT_BINARY_MODE_VALUES = "binary", "ascii"
__OUTPUT_BINARY_MODE_DEFAULT = __OUTPUT_BINARY_MODE_VALUES[0]

def get_vtk_output_keywords():
    return __OUTPUT_FILE, __OUTPUT_BINARY_MODE


def get_vtk_output_help():
    msg = \
    f"""{__OUTPUT_FILE} [string]: The vtk output file destination.
    {__OUTPUT_BINARY_MODE} [string]: For ".vtu" output format, the data mode can be {" or ".join(__OUTPUT_BINARY_MODE_VALUES)}. Defaults to {__OUTPUT_BINARY_MODE_DEFAULT}."""
    return textwrap.dedent(msg)


def fill_vtk_output_subparser(parser) -> None:
    parser.add_argument('--' + __OUTPUT_FILE,
                        type=str,
                        required=True,
                        help=f"[string]: The vtk output file destination.")
    parser.add_argument('--' + __OUTPUT_BINARY_MODE,
                        type=str,
                        metavar=", ".join(__OUTPUT_BINARY_MODE_VALUES),
                        default=__OUTPUT_BINARY_MODE_DEFAULT,
                        help=f"""[string]: For ".vtu" output format, the data mode can be {" or ".join(__OUTPUT_BINARY_MODE_VALUES)}. Defaults to {__OUTPUT_BINARY_MODE_DEFAULT}.""")


def convert(parsed_options) -> VtkOutput:
    output = parsed_options[__OUTPUT_FILE]
    if parsed_options[__OUTPUT_BINARY_MODE] and os.path.splitext(output)[-1] == ".vtk":
        logging.info("VTK data mode will be ignored for legacy file format \"vtk\".")
    is_data_mode_binary: bool = parsed_options[__OUTPUT_BINARY_MODE] == __OUTPUT_BINARY_MODE_DEFAULT
    return VtkOutput(output=output,
                     is_data_mode_binary=is_data_mode_binary)


def parse_cli_options(options: Dict[str, str]) -> VtkOutput:
    output = options[__OUTPUT_FILE]
    if options.get(__OUTPUT_BINARY_MODE) and os.path.splitext(output)[-1] == ".vtk":
        logging.info("VTK data mode will be ignored for legacy file format \"vtk\".")
    is_data_mode_binary = options.get(__OUTPUT_BINARY_MODE, __OUTPUT_BINARY_MODE_DEFAULT) == __OUTPUT_BINARY_MODE_DEFAULT
    return VtkOutput(output=output,
                     is_data_mode_binary=is_data_mode_binary)
