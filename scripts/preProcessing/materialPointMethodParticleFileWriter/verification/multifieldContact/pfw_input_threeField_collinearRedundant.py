"""A|B|C chain with Simple correction, which activates the redundant A-C pair.

Expected: PGS remains finite and converges despite three parallel constraints,
including the geometrically non-neighboring A-C constraint.  This is a focused
reproducer for ordering sensitivity and stagnating residuals.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import make_collinear_chain_case

pfw = make_collinear_chain_case(
    "threeField_collinearRedundant",
    gap_correction="Simple",
)

