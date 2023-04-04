import os.path
import logging
import textwrap

from checks.vtk_utils import VtkOutput


__OUTPUT_FILE = "output"
__OUTPUT_BINARY_MODE = "data-mode"
__OUTPUT_BINARY_MODE_VALUES = "binary", "ascii"
__OUTPUT_BINARY_MODE_DEFAULT = __OUTPUT_BINARY_MODE_VALUES[0]


def get_vtk_output_help():
    msg = \
    f"""{__OUTPUT_FILE} [string]: The vtk output file destination.
    {__OUTPUT_BINARY_MODE} [string]: For ".vtu" output format, the data mode can be {" or ".join(__OUTPUT_BINARY_MODE_VALUES)}. Defaults to {__OUTPUT_BINARY_MODE_DEFAULT}."""
    return textwrap.dedent(msg)


def __build_arg(prefix, main):
    return "-".join(filter(None, (prefix, main)))


def fill_vtk_output_subparser(parser, prefix="") -> None:
    parser.add_argument('--' + __build_arg(prefix, __OUTPUT_FILE),
                        type=str,
                        required=True,
                        help=f"[string]: The vtk output file destination.")
    parser.add_argument('--' + __build_arg(prefix, __OUTPUT_BINARY_MODE),
                        type=str,
                        metavar=", ".join(__OUTPUT_BINARY_MODE_VALUES),
                        default=__OUTPUT_BINARY_MODE_DEFAULT,
                        help=f"""[string]: For ".vtu" output format, the data mode can be {" or ".join(__OUTPUT_BINARY_MODE_VALUES)}. Defaults to {__OUTPUT_BINARY_MODE_DEFAULT}.""")


def convert(parsed_options, prefix="") -> VtkOutput:
    output_key = __build_arg(prefix, __OUTPUT_FILE).replace("-", "_")
    binary_mode_key = __build_arg(prefix, __OUTPUT_BINARY_MODE).replace("-", "_")
    output = parsed_options[output_key]
    if parsed_options[binary_mode_key] and os.path.splitext(output)[-1] == ".vtk":
        logging.info("VTK data mode will be ignored for legacy file format \"vtk\".")
    is_data_mode_binary: bool = parsed_options[binary_mode_key] == __OUTPUT_BINARY_MODE_DEFAULT
    return VtkOutput(output=output,
                     is_data_mode_binary=is_data_mode_binary)
