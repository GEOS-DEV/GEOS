"""Positive explicit gap with bodies approaching but unable to close during the run.

Expected: contact force remains zero and each field retains its initial velocity.
This catches premature contact caused solely by negative relative velocity.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import CELL_SIZE, make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfacePositiveGap",
    gap=0.25 * CELL_SIZE,
    velocities=((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0)),
    end_time=0.02,
)

