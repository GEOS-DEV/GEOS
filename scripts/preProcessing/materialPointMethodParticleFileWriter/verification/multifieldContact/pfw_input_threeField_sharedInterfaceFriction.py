"""Oblique three-field collision with two friction constraints sharing one field.

Expected: normal nonpenetration holds, each tangential impulse remains inside
its Coulomb cone (mu=0.30), and total x/y momentum is conserved.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import CELL_SIZE, make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfaceFriction",
    gap=-0.05 * CELL_SIZE,
    velocities=((0.10, 0.20, 0.0), (0.05, -0.20, 0.0), (-0.10, 0.0, 0.0)),
    friction_coefficient=0.30,
)

