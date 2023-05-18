import logging

from checks.generate_cube import Options, Result, FieldInfo

from . import vtk_output_parsing, generate_global_ids_parsing, GENERATE_CUBE
from .generate_global_ids_parsing import GlobalIdsInfo


__X, __Y, __Z, __NX, __NY, __NZ = "x", "y", "z", "nx", "ny", "nz"
__FIELDS = "fields"


def convert(parsed_options) -> Options:
    def check_discretizations(x, nx, title):
        if len(x) != len(nx) + 1:
            raise ValueError(f"{title} information (\"{x}\" and \"{nx}\") does not have consistent size.")
    check_discretizations(parsed_options[__X], parsed_options[__NX], __X)
    check_discretizations(parsed_options[__Y], parsed_options[__NY], __Y)
    check_discretizations(parsed_options[__Z], parsed_options[__NZ], __Z)

    def parse_fields(s):
        name, support, dim = s.split(":")
        if support not in ("CELLS", "POINTS"):
            raise ValueError(f"Support {support} for field \"{name}\" must be one of \"CELLS\" or \"POINTS\".")
        try:
            dim = int(dim)
            assert dim > 0
        except ValueError:
            raise ValueError(f"Dimension {dim} cannot be converted to an integer.")
        except AssertionError:
            raise ValueError(f"Dimension {dim} must be a positive integer")
        return FieldInfo(name=name, support=support, dimension=dim)

    gids: GlobalIdsInfo = generate_global_ids_parsing.convert_global_ids(parsed_options)

    return Options(vtk_output=vtk_output_parsing.convert(parsed_options),
                   generate_cells_global_ids=gids.cells,
                   generate_points_global_ids=gids.points,
                   xs=parsed_options[__X],
                   ys=parsed_options[__Y],
                   zs=parsed_options[__Z],
                   nxs=parsed_options[__NX],
                   nys=parsed_options[__NY],
                   nzs=parsed_options[__NZ],
                   fields=tuple(map(parse_fields, parsed_options[__FIELDS])))


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(GENERATE_CUBE,
                              help="Generate a cube and its fields.")
    p.add_argument('--' + __X,
                   type=lambda s: tuple(map(float, s.split(":"))),
                   metavar="0:1.5:3",
                   help="[list of floats]: X coordinates of the points.")
    p.add_argument('--' + __Y,
                   type=lambda s: tuple(map(float, s.split(":"))),
                   metavar="0:5:10",
                   help="[list of floats]: Y coordinates of the points.")
    p.add_argument('--' + __Z,
                   type=lambda s: tuple(map(float, s.split(":"))),
                   metavar="0:1",
                   help="[list of floats]: Z coordinates of the points.")
    p.add_argument('--' + __NX,
                   type=lambda s: tuple(map(int, s.split(":"))),
                   metavar="2:2",
                   help="[list of integers]: Number of elements in the X direction.")
    p.add_argument('--' + __NY,
                   type=lambda s: tuple(map(int, s.split(":"))),
                   metavar="1;1",
                   help="[list of integers]: Number of elements in the Y direction.")
    p.add_argument('--' + __NZ,
                   type=lambda s: tuple(map(int, s.split(":"))),
                   metavar="4",
                   help="[list of integers]: Number of elements in the Z direction.")
    p.add_argument('--' + __FIELDS,
                   type=str,
                   metavar="name:support:dim",
                   nargs="+",
                   required=False,
                   default=(),
                   help="Create fields on CELLS or POINTS, with given dimension (typically 1 or 3).")
    generate_global_ids_parsing.fill_generate_global_ids_subparser(p)
    vtk_output_parsing.fill_vtk_output_subparser(p)


def display_results(options: Options, result: Result):
    logging.info(result.info)
