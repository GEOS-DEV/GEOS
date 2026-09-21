"""Coupled three-field overlap with a 1:100:1 density ratio.

Expected: PGS converges within 200 sweeps, conserves normal momentum, and does
not leave a residual above 1e-10.  This targets poor mass-matrix conditioning.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import CELL_SIZE, make_shared_interface_case

pfw = make_shared_interface_case(
    "threeField_sharedInterfaceMassRatio",
    gap=-0.05 * CELL_SIZE,
    velocities=((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0)),
    densities=(1.0, 100.0, 1.0),
)

