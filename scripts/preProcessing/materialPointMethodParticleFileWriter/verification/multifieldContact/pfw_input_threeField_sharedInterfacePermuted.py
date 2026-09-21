"""Baseline overlap with contact-group numbering permuted from (A,B,C) to (2,0,1).

Expected: after mapping fields back to physical bodies, the result matches
sharedInterfaceOverlap.  A difference identifies field/constraint-order bias.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import CELL_SIZE, make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfacePermuted",
    gap=-0.05 * CELL_SIZE,
    velocities=((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0)),
    groups=(2, 0, 1),
)

