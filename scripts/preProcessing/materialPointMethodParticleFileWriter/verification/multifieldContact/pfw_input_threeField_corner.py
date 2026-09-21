"""Three fields meet at a penetrated 2D corner with nonparallel normals.

Expected: all active projected inequalities converge together without rounding
the corner into a single averaged constraint.  This targets nearly dependent
and competing normals at a triple junction.
"""

#[pfw_dependency] input:pfw_three_field_contact_common.py
from pfw_three_field_contact_common import make_corner_case

pfw = make_corner_case("threeField_corner")

