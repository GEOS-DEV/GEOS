import logging
import textwrap

from checks import GENERATE_FRACTURES
from checks.generate_fractures import Options, Result

from . import cli_parsing

__POLICY = "policy"
__FIELD_POLICY = "field"
__POLICIES = (__FIELD_POLICY,)

__FIELD_NAME = "name"
__FIELD_VALUES = "values"

__OUTPUT_FILE = "output"
__GENERATE_FIELD_DATA = "generate_field_data"

__ALL_KEYWORDS = {__POLICY, __FIELD_NAME, __FIELD_VALUES, __FIELD_POLICY, __OUTPUT_FILE, __GENERATE_FIELD_DATA}


def get_help():  # TODO use a formatter module.
    msg = f"""\
    Splits the mesh to generate the faults and fractures. [EXPERIMENTAL]
    
    {__POLICY} [string]: The criterion to define the surfaces that will be changed into fracture zones.
    {" " * len(__POLICY)}           Possible values are "{",".join(__POLICIES)}".
    {__OUTPUT_FILE} [string]: The vtk output destination.
    {__GENERATE_FIELD_DATA} [bool]: Whether we should generate the vtk field data that currently defines the fracture in GEOSX.
    {" " * len(__GENERATE_FIELD_DATA)}         Currently default to true.
    
    If the "{__FIELD_POLICY}" {__POLICY} is selected:
        {__FIELD_NAME} [string]: which field will be considered to define the fractures.
        {__FIELD_VALUES} [list of integers]: which changes of the field will be considered as a fracture.
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for collocated nodes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(GENERATE_FRACTURES, __ALL_KEYWORDS, options)
    policy = options[__POLICY]
    assert policy in __POLICIES
    if policy == __FIELD_POLICY:
        field = options[__FIELD_NAME]
        field_values = tuple(map(int, options[__FIELD_VALUES].split(',')))
    # Do better!
    output = options.get(__OUTPUT_FILE, "")
    return Options(policy=policy, field=field, field_values=field_values, output=output)


def display_results(options: Options, result: Result):
    logging.info("Hell world!")
