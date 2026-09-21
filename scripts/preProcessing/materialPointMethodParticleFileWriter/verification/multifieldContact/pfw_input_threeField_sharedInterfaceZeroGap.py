"""Exact-zero explicit gap with closing velocities.

Regression target: contact should activate when gap <= tolerance.  The current
PGS patch tests gap < 0 exactly, so this deck exposes a one-step activation delay
or missed contact at a mathematically closed interface.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfaceZeroGap",
    gap=0.0,
    velocities=((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0)),
)

