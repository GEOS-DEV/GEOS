"""Interpenetrated explicit surfaces with all pairs separating.

Expected: the unilateral normal multipliers stay zero and no adhesive impulse
is generated.  This isolates complementarity from gap activation.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import CELL_SIZE, make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfaceSeparating",
    gap=-0.05 * CELL_SIZE,
    velocities=((-0.10, 0.0, 0.0), (-0.05, 0.0, 0.0), (0.10, 0.0, 0.0)),
)

