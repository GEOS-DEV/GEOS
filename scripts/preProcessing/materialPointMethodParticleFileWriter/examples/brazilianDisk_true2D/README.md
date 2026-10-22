# Brazilian disk: native 2D CPDI bulk / CPTI surface

This example requires the applied `GEOS-native-2D.patch`. It uses the current
Brazilian-disk FMPM problem: quartz/CeramicDamage with a 2D Voronoi-Weibull
strength field, horizontal velocity 1, periodic x, and moving y boundaries
that compress the height from 1 to 0.8 over time 0 to 1 with cosine interpolation.
It uses the existing material preset's units and loading; this is a dynamic
compression/fragmentation example, not a calibrated quasistatic strength test.

## Install and run

The supplied `GEOS-brazilianDisk-true2D-example.patch` installs this complete
example next to `brazilianDisk_FMPM` and the other Brazilian-disk examples:

```text
scripts/preProcessing/materialPointMethodParticleFileWriter/
  examples/brazilianDisk_true2D/
```

Apply the full example patch from the GEOS repository root. If you already
applied the earlier input-and-README example patch, apply only
`GEOS-brazilianDisk-true2D-workflow.patch` instead. The ZIP also contains a
complete `brazilianDisk_true2D/` directory that can be copied into PFW's
`examples/` directory. No additional C++ changes are included.

With your normal `userDefs_<user>.py` configured in the PFW root, run:

```bash
cd scripts/preProcessing/materialPointMethodParticleFileWriter/examples/brazilianDisk_true2D
./runProblem
```

This follows the existing FMPM example's workflow: copy the input, current PFW
modules, renderer, user definitions, and case scripts to
`testRunDirectory/brazilianDisk_true2D`; generate particles/XML; submit GEOS;
and submit a dependent Slurm post-processing job. It uses the configured GEOS
executable and account, with `pdebug` and five minutes for post-processing.
Use `--post-partition` and `--post-walltime` to change those post-job settings.
Use `--python /path/to/python3` to select Python for both PFW and post-processing.

The checked-in input keeps `mSubmitJobs=False` for direct PFW use.
`runProblem` enables submission **only in the run-directory copy**.
While a build is in progress, generate and inspect the inputs without submitting:

```bash
./runProblem --generate-only
```

Other options include `--prepare-only` (copy files only), `--no-visit` (reaction
plots only), `--no-postprocess` (submit only GEOS), and `--force` (overwrite
previous run/output products without prompting). `--help` lists all options.
The default solver partitioning is 3 x 3 (nine ranks), with an XY grid count of
108 x 72 and two particles per cell per active direction. Nonperiodic y includes
the normal two padding cells in nJ.

## Controls worth editing

| Setting | Default | Meaning |
|---|---|---|
| `formulation` | `PlaneStrain` | Native XY quad mesh and plane-strain constitutive response |
| `planeStrainThickness` | `1/24` | Physical thickness, fixed independently of grid refinement |
| Object `particleType` | `2` | CPDI parallelogram bulk |
| `particleGeometry` | `fromRefinement` | Retain actual cut triangle geometry |
| `surfaceRefinementMethod` | `simplex` | Local surface cuts; the hexTet preset is 3D only |
| `surfaceRefinementLevel` | `2` | Subdivide surface candidates twice per XY direction; try 4 for finer cuts |
| `simplexParticleType` | `3` | CPTI on actual accepted triangle pieces |
| `surfaceRefinementFallback` | `cartesian` | Retain CPDI quads when a cut cannot be resolved |
| `updateMethod`, `updateOrder` | `FMPM`, `2` | Same transfer choice as brazilianDisk_FMPM |
| `cpdiDomainScalingType` | `homel` | CPDI projection / CPTI uniform diameter clipping |
| `geometryReset` | `FOnly` | Preserve material geometry when F is reset |
| `subdivideParticles`, `oversizedParticleTreatment` | `0`, `none` | Runtime splitting disabled for the initial comparison |

Interior particles are CPDI quads. The refined band contains smaller uncut CPDI
quads, CPTI triangles at accepted cuts, and any CPDI fallback pieces. Setting
`particleType=3` on the whole disk would instead tessellate the bulk as well.
The underlying native disk has no z caps, and its surface normals/offsets are
included explicitly in the particle fields. G/G0 is always written.

The thickness equals the existing FMPM example's default numerical thickness,
but its value is independent of `refine` and `cpp`. The Weibull construction
retains the current native example's seed 1, flaw spacing 6*DX, and minimum
statistical volume DX^2*thickness. Changing grid resolution still changes the
flaw spacing, following the existing example; hold that physical spacing fixed
separately if you want a fixed-flaw resolution study.

To enable runtime geometry-aware splitting, use:

```python
pfw['oversizedParticleTreatment'] = 'split'
pfw['defGradResetMaxParticleDomainToGridCellRatio'] = 0.9
pfw['oversizedParticleResetRankParticleCountThreshold'] = -1
```

The threshold is a full particle-domain diameter / active grid-cell diagonal
ratio. Native CPDI and CPTI children retain their own shapes and interpolation.

## Output and generation check

Silo output uses native quad/triangle cells with 200 output intervals. Reaction
and box-average histories use 2,000 intervals. This directory includes:

- `runProblem`: run-directory preparation, PFW generation, GEOS submission,
  and dependent post-processing submission.
- `postProcess_brazilianDisk_true2D`: GEOS completion checks, y-face reaction
  plots, XY VisIt images of initial StrengthScale and final Damage, and copy-back.
- `plotReactions_brazilianDisk_true2D.py`: standalone y-face reaction PNG/CSV.
- `output/`: collected plots, history CSVs, and post-processing logs; generated
  products are ignored by git.

The post script uses the updated shared `pfw_visit_render.py` copied by the
launcher. Matplotlib is required for reaction plots, and VisIt is required for
Silo images. To rerun post-processing from the installed example directory:

```bash
./postProcess_brazilianDisk_true2D \
  --run-dir /path/to/testRunDirectory/brazilianDisk_true2D \
  --output-dir "$PWD/output" --job-id GEOS_JOB_ID
```

Or plot reactions alone:

```bash
python3 ./plotReactions_brazilianDisk_true2D.py \
  --run-dir /path/to/run --output-dir "$PWD/output"
```

For the supplied defaults, PFW generated:

- 19,156 CPDI quads and 2,266 CPTI triangles: 21,422 particles total.
- Physical volume 0.03272492951355, versus pi*(0.5)^2/24 = 0.03272492347489.
- Positive geometry/volumes, compact 2x2 generators and zero z components.
- Identical particle-file bytes on one and three **emulated** writer ranks.
- Native C2D4 mesh XML, PlaneStrain, and both CPDI/CPTI particle subregions.

The bundle's initial-geometry PNG shows the actual serialized material domains;
it is not a GEOS simulation result. `validation.json` contains the generation
checks and refinement report. Volume agreement does not establish globally
conforming particle interfaces. GEOS dynamics and actual MPI were not run here;
use the target native-2D build for those checks. The launcher was checked with
mocked PFW/Slurm commands, and the post-processing reaction path was exercised
using a synthetic history. Actual Slurm submission and VisIt rendering have
not been run here.
