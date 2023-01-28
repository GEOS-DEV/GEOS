import logging
import textwrap

from checks.generate_fractures import Options, Result

from . import cli_parsing, vtk_output_parsing, GENERATE_FRACTURES

__POLICY = "policy"
__FIELD_POLICY = "field"
__POLICIES = (__FIELD_POLICY, )

__FIELD_NAME = "name"
__FIELD_VALUES = "values"

__SPLIT_ON_DOMAIN_BOUNDARY = "split_on_domain_boundary"
__GENERATE_FIELD_DATA = "generate_field_data"

__ALL_KEYWORDS = {
    *vtk_output_parsing.get_vtk_output_keywords(),
    __POLICY, __FIELD_NAME, __FIELD_VALUES, __FIELD_POLICY, __SPLIT_ON_DOMAIN_BOUNDARY, __GENERATE_FIELD_DATA
}


def get_help():    # TODO use a formatter module.
    msg = f"""\
    Splits the mesh to generate the faults and fractures. [EXPERIMENTAL]
    
    {__GENERATE_FIELD_DATA} [bool]: Whether we should generate the vtk field data that currently defines the fracture in GEOSX.
    {" " * len(__GENERATE_FIELD_DATA)}         Currently default to true.
    {__POLICY} [string]: The criterion to define the surfaces that will be changed into fracture zones.
    {" " * len(__POLICY)}           Possible values are "{",".join(__POLICIES)}".
    {__SPLIT_ON_DOMAIN_BOUNDARY} [string]: Split policy if the fracture touches the boundary of the mesh. Defaults to true.
    If the "{__FIELD_POLICY}" {__POLICY} is selected:
        {__FIELD_NAME} [string]: which field will be considered to define the fractures.
        {__FIELD_VALUES} [list of integers]: which changes of the field will be considered as a fracture.
    {vtk_output_parsing.get_vtk_output_help()}
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
    split_on_domain_boundary = options.get(__SPLIT_ON_DOMAIN_BOUNDARY, True)
    assert policy in __POLICIES
    if policy == __FIELD_POLICY:
        field = options[__FIELD_NAME]
        field_values = tuple(map(int, options[__FIELD_VALUES].split(',')))
    vtk_output = vtk_output_parsing.parse_cli_options(options)
    return Options(policy=policy,
                   field=field,
                   field_values=field_values,
                   vtk_output=vtk_output,
                   split_on_domain_boundary=split_on_domain_boundary)


def display_results(options: Options, result: Result):
    logging.info("Hell world!")
