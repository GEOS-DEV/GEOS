"""Baseline A|B|C chain: two neighboring constraints share the middle field.

Expected: the common-node normal velocities converge toward the conserved
center-of-mass velocity (zero here), and the PGS solve requires multiple sweeps.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import make_collinear_chain_case

pfw = make_collinear_chain_case("threeField_collinearChain")

