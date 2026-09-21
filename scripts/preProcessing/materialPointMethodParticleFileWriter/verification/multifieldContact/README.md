# Three-field material-contact verification

These compact plane-strain inputs exercise the CPU-only
`ProjectedGaussSeidel` material-contact path.  Every case uses exactly three
contact groups, PIC transfer, CPDI particles, no DFG, and strict PGS controls:

```python
pfw["contactSolver"] = "ProjectedGaussSeidel"
pfw["contactPGSMaximumIterations"] = 200
pfw["contactPGSVelocityTolerance"] = 1.0e-10
pfw["contactPGSRelaxation"] = 1.0
pfw["contactPGSRequireConvergence"] = 1
```

## Test matrix

| Input | Geometry / loading | Edge case and expected check |
|---|---|---|
| `pfw_input_threeField_collinearChain.py` | A\|B\|C, `+v,0,-v`; mapped centers gate contact | Baseline two-constraint coupled solve. A-B and B-C should activate, A-C should not, the common-node velocities should approach zero, and more than one PGS sweep should be needed. |
| `pfw_input_threeField_collinearRedundant.py` | Same chain with `Simple` correction | Deliberately activates the non-neighbor A-C pair. Checks redundant parallel constraints, order sensitivity, and residual stagnation. |
| `pfw_input_threeField_collinearExplicitSurfaces.py` | Same chain using explicit surface positions | Diagnostic/expected-failure target. One field has two opposing faces but only one mapped nodal surface position; averaging can make both neighbor gaps look positive. This separates a surface-representation failure from a PGS failure. |
| `pfw_input_threeField_sharedInterfaceOverlap.py` | Two left fields push into one right field across a slightly overlapped planar interface | Clean explicit-surface baseline. Both left-right constraints must be active, normal relative velocities must be nonnegative after projection, and total momentum must be conserved. |
| `pfw_input_threeField_sharedInterfaceZeroGap.py` | Same shared interface at exact zero gap | Activation-boundary regression. Contact should be active for closing motion at `gap <= tolerance`; the current exact `gap < 0` condition can delay activation by one step. |
| `pfw_input_threeField_sharedInterfacePositiveGap.py` | Positive gap, approaching velocities, run ends before closure | No premature contact. `gridContactForce` should stay zero and velocities should remain unchanged despite negative relative velocity. |
| `pfw_input_threeField_sharedInterfaceSeparating.py` | Slight geometric overlap but separating velocities | Unilateral complementarity. Normal multipliers must remain zero; contact must not become adhesive merely because the gap is negative. |
| `pfw_input_threeField_sharedInterfacePermuted.py` | Overlap baseline with physical groups numbered `(2,0,1)` | Field-order invariance. Map fields back to bodies and compare with `sharedInterfaceOverlap`; results should agree to the solver tolerance. |
| `pfw_input_threeField_sharedInterfaceMassRatio.py` | Overlap with density ratio `1:100:1` and density-scaled stiffness | Conditioning. The node solve must converge within 200 sweeps while preserving momentum and the projected inequalities. |
| `pfw_input_threeField_sharedInterfaceFriction.py` | Oblique impact, `mu=0.30` | Coupled Coulomb projections. Each tangential impulse must satisfy `|J_t| <= mu J_n`, with x/y momentum conserved. |
| `pfw_input_threeField_corner.py` | Three bodies at a penetrated 2D corner | Nonparallel and competing normals at a triple junction. Checks that all projected inequalities converge without collapsing them into one averaged direction. |

## Important geometry distinction

Three fields need nonzero mapped mass at the **same grid node**, not necessarily
particle centers in the same background cell.  The A\|B\|C case makes the middle
strip one cell wide and centers it on the origin so that node receives CPDI
support from all three fields.  It also fixes the middle field's mapped normal to `+x`:
otherwise its two opposing face normals can cancel at that node.

The shared-interface geometry is the preferred pass/fail integration test for
grid-mapped surface positions.  It gives each field only one relevant contact
surface at the junction and therefore avoids the one-surface-position-per-field
ambiguity present in a literal A\|B\|C sandwich.

## Suggested evaluation order

1. Run `sharedInterfaceOverlap`, `sharedInterfacePositiveGap`, and
   `sharedInterfaceSeparating` first to establish activation and complementarity.
2. Compare `sharedInterfacePermuted` against the overlap baseline.
3. Run `sharedInterfaceMassRatio`, `sharedInterfaceFriction`, and `corner` as
   convergence stress tests.
4. Run the three collinear decks last.  Treat `collinearExplicitSurfaces` as a
   diagnostic rather than a required pass until the middle field can retain
   pair-specific surface positions.

All decks write `gridActive`, `gridMass`, `gridUncontactedVelocity`,
`gridVelocity`, `gridContactForce`, `gridCenterOfMass`, `gridSurfaceNormal`,
and `gridSurfacePosition` to Silo.  The uncontacted velocity is required to
measure the same-step contact impulse and momentum balance without confusing
contact with the rest of the explicit update.
Run on a CPU build: the current PGS implementation intentionally rejects device
builds.

Run one case through the standard verification harness with, for example:

```bash
./runProblem --case threeField_sharedInterfaceOverlap
```

## Automated post-processing

`postProcess_threeFieldContact.py` is both a case-local reducer and a suite
report generator.  The normal verification harness discovers it automatically
after each case, calls VisIt to extract the origin-junction history, and writes
the per-case products into that case's `output/` directory.

After all cases finish, compile the aggregate report with:

```bash
./runPostProcess --visit-cmd /path/to/visit
```

The script discovers run directories from each `*_jobs.json` manifest.  If the
per-case extraction has already run, VisIt is not needed a second time:

```bash
./runPostProcess --no-visit
```

Generated products are:

- `three_field_contact_report.html`: self-contained pass/fail report with an
  expected-versus-actual figure for every case;
- `three_field_contact_report.md`: review-friendly text report;
- `three_field_contact_summary.json`: machine-readable case/check hierarchy;
- `three_field_contact_checks.csv`: one row per numerical check;
- `threeFieldContact_results.tex`: compact status macros for a larger report;
- `figures/*.png`: pre/post/oracle velocity bars, constraint margins, contact
  impulse, and momentum-error histories.

The required checks include finite output, same-step momentum conservation,
the pairwise projected inequalities, contact/no-contact complementarity, the
zero-gap activation boundary, Coulomb cones, two-axis corner response, and the
published local velocity oracle when the sampled state still matches the
prescribed initial nodal masses and velocities.  The permuted case is also
compared directly against the overlap baseline after mapping raw fields back to
physical bodies A/B/C.  `collinearExplicitSurfaces` is reported as
`DIAGNOSTIC`, not as a suite failure, when it reproduces the documented mapped
surface-position limitation.

Exit status is `0` for a completed passing suite (diagnostic classification is
allowed), `1` for a failed numerical check, and `2` for missing/incomplete run
data.  Use `--allow-incomplete` only when generating a preflight report.

The postprocessor's own synthetic regression checks can be run with:

```bash
python3 -m unittest -v test_postProcess_threeFieldContact.py
```

The included `particleFileWriter_pgs_keys.patch` is needed only if the local
particle-file writer has not yet added the five `contactSolver`/`contactPGS*`
keys to its accepted parameter table.
