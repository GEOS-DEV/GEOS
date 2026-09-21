"""Diagnostic A|B|C chain using the mapped explicit surface positions.

Regression target: one nodal field stores one surface position.  The two faces
of the middle strip can average together at the common node and make both true
neighbor gaps appear positive.  Compare this deck with collinearChain; a missed
contact here identifies a surface-representation limitation, not a PGS failure.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import make_collinear_chain_case

pfw = make_collinear_chain_case(
    "threeField_collinearExplicitSurfaces",
    use_surface_positions=True,
)

