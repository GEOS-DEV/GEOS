# Lagrangian Contact Mechanics Robustness Cases

This folder contains robustness test cases for the Lagrangian contact mechanics formulations on unit-cube meshes with internal fractures.

All meshes are defined on the domain `[0,1] \times [0,1] \times [0,1]`.

## Case Summary

Different physics configurations (C1, C2, C3) are tested on both Hexahedral (`rc_hex_*`) and Tetrahedral (`rc_tet_*`) meshes.

| Case ID | XML File Pattern | Mesh Type | Physics | Description | BCs/Loads |
| --- | --- | --- | --- | --- | --- |
| **C1** | `rc_hex_c1_lagr_contact.xml`<br>`rc_tet_c1_lagr_contact.xml` | Hex / Tet | **Pure Mechanical**<br>(SolidMechanicsAugmentedLagrangianContact) | Tensile simulation where the boundary is pulled in each canonical direction to open fractures. Tests splits of nodal, edge, facet, and cell-based discretizations. | **Displacement Control:**<br>Faces pulled outwards (+/- 1e-4 m).<br>Zero initial stress. |
| **C2** | `rc_hex_c2_lagr_contact.xml`<br>`rc_tet_c2_lagr_contact.xml` | Hex / Tet | **Pure Flow**<br>(SinglePhaseFVM) | Mixed-dimensional flow simulation. Linear pressure solution expected (exact on Hex). | **Pressure Control:**<br>P_in (xneg) = 2 MPa<br>P_out (xpos) = 1 MPa<br>Linear gradient. |
| **C3** | `rc_hex_c3_lagr_contact.xml`<br>`rc_tet_c3_lagr_contact.xml` | Hex / Tet | **Coupled Poromechanics**<br>(SinglePhasePoromechanicsConformingFracturesALM) | Fully coupled poromechanical simulation testing two-way coupling stability. | **BCs:**<br>Confining Prescribed Displacement (0.0)<br>Injection Pressure in Fracture (40 MPa)<br>Initial Stress (-30 MPa). |

## File Naming Convention

- `rc_{mesh_type}_{case_id}_lagr_contact.xml`
    - `mesh_type`: `hex` or `tet`
    - `case_id`: `c1` (Mech), `c2` (Flow), `c3` (Coupled)

## Notes
- All meshes are unit cubes with a Discrete Fracture Network (DFN).
- Flow cases (C2) use dummy solid/poro couplings in constitutive definitions to satisfy solver requirements.
