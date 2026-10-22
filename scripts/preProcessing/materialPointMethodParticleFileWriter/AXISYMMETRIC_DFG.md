# Axisymmetric ordinary multi-group contact and DFG

This patch is limited to axisymmetry. It removes the solver/PFW DFG rejection
and adapts damage-gradient reconstruction to the meridional section. The
shared separability criterion, its defaults, neighbor-selection policy, field
assignment, contact laws and CZ implementation are unchanged. General
multi-field CZ and separability refinements are being developed separately.

## Geometry and weights

Use x=r, y=axial and z=0. Stress/F/L zz remains the hoop component. For the
normalized damage kernel, use current meridional area A=V/(2*pi*r), where V
is current full-ring volume and r is positive particle radius:

    dbar = sum(A*d*K)/sum(A*K)
    grad(dbar) = sum(A*(d-dbar)*grad(K))/sum(A*K)

K is the existing compact cubic kernel, proportional to
`(1-q)^2*(1+2*q)` for `q=distance/neighborRadius < 1`. Its common normalization
cancels. A reflected sample at (-r,y) enforces even scalar parity across the
axis and zero radial gradient at r=0. Images contribute only inside support;
they do not create particles or mass. For nonnegative query/sample radii,
a reflected point cannot be closer than its real point, so no extra neighbor
search is required. The existing surface-as-damage convention is retained.

Only this geometric reconstruction uses area. Mapped damage remains
mass-weighted; normal-quality tensors and single-field fractions remain
volume-weighted. Mass, momentum, volume, contact forces and impulses keep
full-ring measures. The existing cylindrical contact normal uses
`sum(V*(G_r+T, G_y, 0))`; this patch adds no circumference factor to contact.

Cartesian damage-gradient calculation is unchanged.

## Contact fields and insertion

- Ordinary contact already supports more than two contact groups and loops
  over the occupied field pairs in the standard and FMPM paths.
- With DFG, the existing field index is `side*numContactGroups + group`.
  Three groups therefore allocate six slots; unoccupied slots carry no contact.
- Existing DFG rules decide whether two same-group fields tie or can separate
  and recontact. Their implementation and defaults are untouched.
- Contact-surface insertion does not automatically create cohesive bonds.
- CZ retains its existing two-group, DFG-off restriction. General multi-field
  CZ storage and history are outside this patch.

The existing partition has one reference direction per node and two slots
per contact group. This does not add arbitrary unresolved crack junctions.
Damage evolution must respect axisymmetry: transverse, cylindrical or inclined
meridional cracks can revolve into valid surfaces; localized longitudinal
hoop fractures cannot. Weibull remains deferred.

Solver-reconstructed surface positions/nodal areas, SPH and overlap corrections
remain unsupported. For damage-created surfaces without supplied positions,
use `useSurfacePositionForContact=0` and the existing approximate gap criterion.

## Example and checks

`examples/axisymmetricDFG` supplies four cases (`intact`, `damage`, `contact`,
`combined`) with CPDI core/CPTI skin, FLIP or FMPM, a runner, particle checker
and tracer plots. The damaged band is prescribed to isolate the existing DFG
path; no fracture constitutive model is added. See the example README.

From the PFW directory:

```sh
python3 tests/run_axisymmetric_dfg_cpp_tests.py
python3 tests/run_axisymmetric_dfg_writer_tests.py --output /tmp/axis-dfg-writer
python3 tests/run_axisymmetric_dfg_cpp_tests.py --particle-files \
  /tmp/axis-dfg-writer/intact_FLIP_r1/mpmParticleFile_axisymmetricDFG \
  /tmp/axis-dfg-writer/damage_FLIP_r1/mpmParticleFile_axisymmetricDFG \
  /tmp/axis-dfg-writer/contact_FLIP_r1/mpmParticleFile_axisymmetricDFG \
  /tmp/axis-dfg-writer/combined_FLIP_r1/mpmParticleFile_axisymmetricDFG
```

Use a fresh output directory. The C++ checks compile actual production
reconstruction, mapping, partition, separability and contact function bodies
with lightweight CPU container adapters, strict warnings and address/undefined
behavior sanitizers. They check constant damage, kernel derivatives, axis
parity, radial translation away from the axis, surface flags, three groups/six
slots, and equal-and-opposite contact impulses. Optional example checks use
a serial mapping adapter to exercise the unchanged separation criterion at
the central damage band; they do not run a GEOS timestep.

Writer checks cover four cases, two update methods and 1/3 emulated ranks.
They verify groups, prescribed damage, absence of CZ/explicit surfaces, mixed
geometry, ring volume, momentum, XML options and byte-identical particles
across rank counts. Existing ordinary-contact checks remain applicable.

Local verification passed: 16 DFG writer runs, four ordinary-contact writer
regressions, the new DFG CPU tests and the existing 540 contact-law comparisons.
The generated damaged/combined bands passed the unchanged separation test at
all seven sampled midplane nodes; the intact/contact cases had no same-group
opening there. Source comparison confirms the separability, quality, partition
and pairwise contact functions are byte-identical to the baseline.

A full GEOS build and dynamics run, actual MPI migration/synchronization,
restart, GPU execution and output readback are target-machine acceptance
steps. These lightweight checks do not establish dynamics accuracy or
convergence. The patch updates manual LaTeX, not a generated manual PDF.
