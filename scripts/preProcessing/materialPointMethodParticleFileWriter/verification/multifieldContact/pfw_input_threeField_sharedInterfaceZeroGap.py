"""Exact-zero explicit gap with closing velocities.

Regression target: implicit contact must include the pair when its mapped gap
is within the configured activation tolerance. The unilateral projection must
then stop closing motion without waiting for one step of interpenetration.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfaceZeroGap",
    gap=0.0,
    velocities=((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0)),
)
