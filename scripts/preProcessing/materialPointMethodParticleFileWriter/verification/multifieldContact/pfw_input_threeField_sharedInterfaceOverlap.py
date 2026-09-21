"""Baseline explicit-surface case: two left fields overlap one right field.

Expected: both left-right constraints are active at the junction node, all
normal post-contact relative velocities are nonnegative, and momentum is
conserved.  This is the cleanest three-field PGS integration check.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import CELL_SIZE, make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfaceOverlap",
    gap=-0.05 * CELL_SIZE,
    velocities=((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0)),
)

