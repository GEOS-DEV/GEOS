import argparse
import logging

from checks.generate_fractures import Options, Result

from . import vtk_output_parsing, GENERATE_FRACTURES

__POLICY = "policy"
__FIELD_POLICY = "field"
__POLICIES = (__FIELD_POLICY, )

__FIELD_NAME = "name"
__FIELD_VALUES = "values"

__SPLIT_ON_DOMAIN_BOUNDARY = "split_on_domain_boundary"
__GENERATE_FIELD_DATA = "generate_field_data"


def fill_subparser(subparsers) -> argparse.ArgumentParser:
    p: argparse.ArgumentParser = subparsers.add_parser(GENERATE_FRACTURES,
                                                       help="Splits the mesh to generate the faults and fractures. [EXPERIMENTAL]")
    p.add_argument('--' + __GENERATE_FIELD_DATA,
                   type=bool,
                   metavar="True",
                   default="True",
                   help=f"[bool]: Whether we should generate the vtk field data that currently defines the fracture in GEOSX. Defaults to \"True\".")
    p.add_argument('--' + __POLICY,
                   type=str,
                   metavar=",".join(__POLICIES),
                   required=True,
                   help=f"[string]: The criterion to define the surfaces that will be changed into fracture zones. Possible values are {', '.join(__POLICIES)}")
    p.add_argument('--' + __SPLIT_ON_DOMAIN_BOUNDARY,
                   type=bool,
                   metavar="True",
                   default=True,
                   help=f"[bool]: Split policy if the fracture touches the boundary of the mesh. Defaults to true.")
    p.add_argument('--' + __FIELD_NAME,
                   type=str,
                   help=f"[string]: If the \"{__FIELD_POLICY}\" {__POLICY} is selected, defines which field will be considered to define the fractures.")
    p.add_argument('--' + __FIELD_VALUES,
                   type=str,
                   nargs="+",
                   help=f"[string]: If the \"{__FIELD_POLICY}\" {__POLICY} is selected, which changes of the field will be considered as a fracture.")
    vtk_output_parsing.fill_vtk_output_subparser(p)
    return p


def convert(parsed_options) -> Options:
    policy = parsed_options[__POLICY]
    split_on_domain_boundary = parsed_options[__SPLIT_ON_DOMAIN_BOUNDARY]
    assert policy in __POLICIES
    if policy == __FIELD_POLICY:
        field = parsed_options[__FIELD_NAME]
        field_values = tuple(parsed_options[__FIELD_VALUES])
    vtk_output = vtk_output_parsing.convert(parsed_options)
    return Options(policy=policy,
                   field=field,
                   field_values=field_values,
                   vtk_output=vtk_output,
                   split_on_domain_boundary=split_on_domain_boundary)


def display_results(options: Options, result: Result):
    logging.info("Hell world!")
