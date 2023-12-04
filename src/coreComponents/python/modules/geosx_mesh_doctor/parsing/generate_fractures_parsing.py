import logging

from checks.generate_fractures import Options, Result

from . import vtk_output_parsing, GENERATE_FRACTURES

__POLICY = "policy"
__FIELD_POLICY = "field"
__INTERNAL_SURFACES_POLICY = "internal_surfaces"
__POLICIES = (__FIELD_POLICY, __INTERNAL_SURFACES_POLICY )

__FIELD_NAME = "name"
__FIELD_TYPE = "type"
__FIELD_VALUES = "values"

__FIELD_POINT_TYPE = "point"  # Move this to the check part...
__FIELD_CELL_TYPE = "cell"
__FIELD_CELL_TO_FACES_TYPE = "cell_to_faces"
__FIELD_DEFAULT_TYPE = __FIELD_CELL_TYPE
__FIELD_TYPES = (__FIELD_POINT_TYPE, __FIELD_CELL_TYPE, __FIELD_CELL_TO_FACES_TYPE)

__SPLIT_ON_DOMAIN_BOUNDARY = "split_on_domain_boundary"

__FRACTURE_PREFIX = "fracture"


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(GENERATE_FRACTURES,
                              help="Splits the mesh to generate the faults and fractures. [EXPERIMENTAL]")
    p.add_argument('--' + __POLICY,
                   type=str,
                   metavar=", ".join(__POLICIES),
                   required=True,
                   help=f"[string]: The criterion to define the surfaces that will be changed into fracture zones. Possible values are \"{', '.join(__POLICIES)}\"")
    p.add_argument('--' + __SPLIT_ON_DOMAIN_BOUNDARY,
                   type=bool,
                   metavar="True",
                   default=True,
                   help=f"[bool]: Split policy if the fracture touches the boundary of the mesh. Defaults to true.")
    p.add_argument('--' + __FIELD_NAME,
                   type=str,
                   help=f"[string]: If the \"{__FIELD_POLICY}\" {__POLICY} is selected, defines which field will be considered to define the fractures. If the \"{__INTERNAL_SURFACES_POLICY}\" {__POLICY}  is selected, defines the name of the attribute will be considered to identify the fractures. ")
    p.add_argument('--' + __FIELD_TYPE,
                   type=str,
                   default=__FIELD_DEFAULT_TYPE,
                   help=f"[string]: The type of field we're going to split on. Options are \"{', '.join(__FIELD_TYPES)}\". Defaults to \"{__FIELD_DEFAULT_TYPE}\".")
    p.add_argument('--' + __FIELD_VALUES,
                   type=str,
                   help=f"[list of comma separated integers]: If the \"{__FIELD_POLICY}\" {__POLICY} is selected, which changes of the field will be considered as a fracture. If the \"{__INTERNAL_SURFACES_POLICY}\" {__POLICY} is selected, list of the fracture attributes.")
    vtk_output_parsing.fill_vtk_output_subparser(p)
    vtk_output_parsing.fill_vtk_output_subparser(p, prefix=__FRACTURE_PREFIX)


def convert(parsed_options) -> Options:
    policy = parsed_options[__POLICY]
    split_on_domain_boundary = parsed_options[__SPLIT_ON_DOMAIN_BOUNDARY]
    assert policy in __POLICIES
    field = parsed_options[__FIELD_NAME]
    field_values = frozenset(map(int, parsed_options[__FIELD_VALUES].split(",")))
    if policy == __FIELD_POLICY:
        field_type = parsed_options[__FIELD_TYPE]
    elif policy == __INTERNAL_SURFACES_POLICY:
        field_type = ""
    vtk_output = vtk_output_parsing.convert(parsed_options)
    vtk_fracture_output = vtk_output_parsing.convert(parsed_options, prefix=__FRACTURE_PREFIX)
    return Options(policy=policy,
                   field=field,
                   field_type=field_type,
                   field_values=field_values,
                   vtk_output=vtk_output,
                   vtk_fracture_output=vtk_fracture_output,
                   split_on_domain_boundary=split_on_domain_boundary)


def display_results(options: Options, result: Result):
    pass
