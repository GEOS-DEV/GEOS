# MPM source and documentation audit

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
