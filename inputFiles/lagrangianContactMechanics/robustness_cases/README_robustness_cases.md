# Lagrangian Contact Mechanics Robustness Cases

This folder contains robustness test cases for the Lagrangian contact mechanics formulations on unit-cube meshes with internal fractures.

All meshes are defined on the domain `[0,1] \times [0,1] \times [0,1]`.

## Case Summary

| Case | Base XML | Contact XML | Mesh | Mesh Type | Domain | Fracture | Material | BC/Load | Solver | Robustness Focus |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RC_HEX01_NORM_NF_LAGR | `cases/rc_hex01_norm_nf_lagr_base.xml` | `cases/rc_hex01_norm_nf_lagr_contact.xml` | `meshes/fractured_mesh_hex.vtu` | hex | `[0,1]^3` | single internal fracture (as in mesh) | `ElasticIsotropic(E=1e9, nu=0.25); Coulomb(mu=0, c=0)` | `x0,y0,z0 fixed; compressive normal traction on z1 (ramp 0→1 MPa)` | `SolidMechanicsLagrangeContact, TPFAstabilization, newtonTol=1e-8` | baseline robustness on hex mesh, zero friction |

## Notes

- All robustness cases use zero gravity unless otherwise stated.
- Additional cases can be added under `cases/` and referenced here with their main characteristics.
