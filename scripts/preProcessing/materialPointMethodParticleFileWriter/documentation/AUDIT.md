# Geometry-aware splitting update (2026-10-21)

The current source supports native tet/wedge runtime subdivision. Documentation
below describing native splitting as unavailable is historical and superseded.
Updated the splitting theory, controls, unified-geometry contract and examples;
corrected full-edge box feasibility/length conventions in the splitting section.
Documented force partition, specific-vs-total energy, melt/history transfer,
FOnly surface frames, pass limits and remaining CZ/CPTI/2D scope. Added standalone
production-header and extracted-solver tests with strict warnings. Nine files
missing from the uploaded ZIP are restored, including the required geometry header.

---

# Current unified-geometry audit (2026-10-21)

Current contract: schema 2 only; shape + full reference/current generators;
temporary CPDI projection; FOnly/FAndGeometry reset policies. Old PFW Python
inputs normalize to new files/XML. Old GEOS files and restarts are unsupported.
The sections below are historical audit records and are superseded on geometry
storage, file compatibility, reset behavior and output by
`../UNIFIED_PARTICLE_GEOMETRY.md` and manual section `unified-particle-geometry`.

Updated: reader/writer, geometry, CPDI mapping, box splitting, reset/healing,
rigid rotation, TypeB CZ anchors, Silo/VTK, PFW aliases, examples, source-derived
option/schema tables, geometry/refinement/output/robustness/CZ theory.

Checks: Python geometry/refinement suite; standalone production geometry and
CPDI kernel bodies; full writer with 1/3/7 emulated ranks and C++ record validation.
Full GEOS build, device execution, actual MPI migration, new restart round-trip,
and coupled CZ mechanics still need validation in the configured GEOS environment.
CPTI, runtime simplex splitting, native 2D and axisymmetry are future extensions.

---

# MPM source and documentation audit

## Current unified geometry patch (GEOS-mpm-4.zip)

`GEOS-mpm-4-unified-geometry.patch` applies to the current v4 upload, including
the earlier merge/build fix and PFW progress/comment updates. It centralizes
affine geometry operations, updates output and box splitting to use them,
separates mapping capacity from geometry vertex count, and adds an opt-in
`particleGeometryFormat="unified"` file schema. See
[`../UNIFIED_PARTICLE_GEOMETRY.md`](../UNIFIED_PARTICLE_GEOMETRY.md).

Existing runtime arrays/checkpoint names remain. Runtime memory reduction, CPTI,
tet splitting, true 2D and axisymmetry are not implemented by this first patch.
The narrative, generated references/source fingerprints, and PDF were rebuilt.

Validation: 32 Python unit tests; 3,000 C++ affine/reset trials; 5,000 domain
checks with 34,000 box children; 96 refined writer executions using emulated
1/3/7-rank collectives; 16 byte-identical comparisons against the pristine v4
writer; eight no-refinement writer executions; and a 2,277-particle complete
example. The production C++ codec reproduced legacy R/native generators for
7,622 particles across 20 corresponding serial cases. Explicit volumes agreed
with legacy determinant volumes within floating-point tolerance.

The standalone tests do not compile the complete GEOS solver or validate its
field registration, Silo/VTK runtime, MPI transport/migration, restart, or GPU
execution. The archive does not contain the required dependency-submodule
contents. Those remain target-build integration checks.

Apply this incremental patch from the unmodified v4 repository root:

    git apply --check /path/to/GEOS-mpm-4-unified-geometry.patch
    git apply /path/to/GEOS-mpm-4-unified-geometry.patch

The earlier v3 audit below is historical; its patch commands describe the
previous merge and must not be reapplied to the v4 baseline.

## Earlier v3 audit record

Audited 2026-10-20 against the selective MPM update based on GEOS-mpm-3.zip.
The source version file still identifies v1.1.0; that identifier alone does not
identify this patch. Generated metadata includes SHA-256 fingerprints of its
principal solver, enum, reader, region, writer, refinement, and restart inputs.

## Applying the complete patch

The single GEOS-mpm-complete-update.patch is relative to the root of the
unmodified GEOS-mpm directory extracted from GEOS-mpm-3.zip. It includes source,
PFW, narrative documentation, generated references, and the rebuilt PDF.
It replaces the earlier selective-merge patch sequence; do not apply both.

From that baseline repository root:

    git apply --check /path/to/GEOS-mpm-complete-update.patch
    git apply /path/to/GEOS-mpm-complete-update.patch

Git binary-patch support is required for the PDF. A plain POSIX patch command
does not apply that binary payload. Keep local changes separate and resolve
any failed baseline check before applying.

## Source scope retained

- Two-pass initial PFW refinement, including Cartesian and local simplex cuts,
  empty-center recovery, material-factor handling, rank slicing, and diagnostics.
- Removal of the legacy weak-trace/mortar implementation. The newer feature
  branch mortar implementation is excluded.
- Approved option A: native wedge/tetrahedron reference geometry, reader checks,
  reference rebasing, and Silo/VTK output.
- Approved option B: origin-centered spherical/circular confining-pressure region.
- Approved option C: portable runner dependency copies, source-root dispatch, and
  paired restart selection.
- Approved option D: general VisIt field selection, region expressions, vector
  overlays, and report rendering support.

The documentation audit adds one PFW parameter-table entry:
applyVelocityGradientWhenSplittingParticles. This already exists in GEOS; PFW
now accepts and emits explicit values. Its PFW default is None, preserving the
solver default of 1 without adding an XML attribute to existing inputs.
No additional solver-mechanics change was made during the documentation audit.

## Corrections

| Area | Current documented behavior |
|---|---|
| Initial refinement | Complete option/default reference, material indexing, factor composition, activation at simplex level 1, fallback/tolerance, rank factories, and diagnostics. |
| Native geometry | Optional ten-column group; native versus equivalent-box output; mechanics remain based on legacy R-vectors; native runtime splitting is guarded. |
| Runtime splitting | Conditional general policy and independent oversized policy; positive ratio requirement; strict rank-count comparison and start-of-step active master count; cohesive exclusions; one bisection per active direction per step. |
| Velocity and conservation | Affine child velocities conserve linear momentum but can add kinetic energy; option 0 copies parent velocities. Mass and current/reference volume partitioning are described separately. |
| PFW validation | Unknown keys raise KeyError. Missing required surface fields raise ValueError. None-valued solver keys are omitted. Obsolete useAPIC examples were removed. |
| Kinematics | Frozen-endpoint matrix exponential with separately integrated volume; FSubcycles is deprecated and ignored. |
| Cohesive setup | Preregistered CohesiveZoneRegions plus ReferenceCohesiveZones; event default Nodal; effective volume-normalization setting is on the region. The registered same-named event value is not forwarded by the current handler. |
| Confining pressure | Positive radius supersedes the box; x/y distance in plane strain, x/y/z in 3D, strict exterior test. |
| Restart | Matching directory/root-file pairs, width search from 16 to 6 digits, root-file modification-time selection, no continuation when no pair exists. |
| Postprocessing | Strict scalar selection, repeatable vector overlays, source renderer preference, fresh-frame preference with existing-frame fallback after a successful render. |
| References | 219 PFW keys, 215 direct solver wrapper registrations, 20 event classes, 17 enums, and 55 particle-file header fields regenerated from source. Conditional geometry fields and serialized enum spellings are explicit. |

Cached MPM schema tables are no longer reused as a competing default source in
the manual. The MPM and ParticleRegion RST references now list current direct
registrations and explicitly identify their inheritance limitation. Generating a
complete inherited XML schema still requires a GEOS build.

## Validation

The documentation pass checked every literal PFW key used in the editable manual
against the current parameter table, verified the generated source fingerprints,
and exercised the complete writer in seven single-process cases with an explicitly
emulated collective API:

1. Default input, confirming omission of the optional child-velocity attribute.
2. Copied child velocity with independent oversized splitting and rank fallback.
3. Explicit affine child velocity.
4. Two-material Cartesian refinement.
5. Native plane-strain wedge metadata with runtime splitting disabled.
6. Native 3D tetrahedral metadata with equivalent-box output.
7. Preregistered cohesive regions and a ReferenceCohesiveZones event.

All seven writer runs and XML assertions passed. They generate files and use
/bin/true as the executable placeholder; they do not execute GEOS.

Earlier selective-port validation, unchanged by the documentation work, included
26 Python geometry/ownership tests, 3,000 randomized production-header geometry
and reference-reset trials, 24 native writer cases across 1/3/7 emulated ranks,
and standalone checks of pressure and diagnostic behavior. These are targeted
checks, not a complete GEOS build.

The LaTeX manual was rebuilt with refreshed metadata, cross-references and index.
Long identifiers and tables were reformatted to avoid clipping. The complete
patch is checked by applying it to a pristine extraction of the original archive
and comparing all resulting file bytes with the intended updated tree.

## Limits

The supplied archive lacks dependency-submodule contents needed for a complete
GEOS build. CPU/GPU solver compilation, real MPI transport, migration/restart,
native Silo/VTK runtime output, and actual VisIt execution remain integration
checks. Native wedge/tetrahedral runtime splitting remains unsupported and errors
when selected. The audit reconciles controls, file formats, and the relevant
implementation paths; it does not independently rederive or numerically validate
every constitutive model or regenerate historical V&V report PDFs.

## Maintaining the reference

From this documentation directory:

    bash build_manual.sh --refresh-data

This updates source metadata, direct-wrapper RST references, generated LaTeX,
and the PDF. Without --refresh-data the build uses the checked-in metadata.
Narrative files remain editable and are never overwritten by generation.
LaTeX intermediates are ignored and removed from the distributed baseline patch.
