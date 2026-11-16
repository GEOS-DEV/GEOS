# GEOS MPM — thermal-transport verification suite

This directory contains **32 separately runnable PFW verification cases**, each with
`pfw_input_<case>.py`, `postProcess_<case>.py`, `runTest`, and `test.tex`.
The shared files `thermal_pfw.py`, `thermal_cases.py`, and `thermal_post.py`
implement the test geometry, analytical references and numerical checks.

**Thermal solver API:** Inputs target the `geos_mpm_thermal_extension.patch`
interface (`<Thermal>` child of `SolidMechanics_MPM`, added to PFW by
`pfw["thermal"]`, with `ThermalElasticIsotropic` or
`LinearThermalElasticIsotropic` constitutive materials). This API was verified
against the extension patch in the file library, **not against the exact current
post-patch GEOS working tree on your HPC system**. If that tree has diverged,
adapt `thermal_pfw.py`. In particular, the old ZIP contains no enabled thermal
solver; you need the new extension to execute the generated inputs.

## Installation

From the repository root, use the accompanying Git patch, or unpack the ZIP
such that the directory exists at

```text
scripts/preProcessing/materialPointMethodParticleFileWriter/verification/thermalMPM/
```

A standard GEOS PFW environment with NumPy, mpi4py and the other dependencies
of `particleFileWriter.py` is required to **generate/run actual solver jobs**.
Analytical and post-processing self-tests use Python's standard library.

## Start here

```bash
cd scripts/preProcessing/materialPointMethodParticleFileWriter/verification/thermalMPM
python3 -m unittest -q test_thermal_references
python3 runSuite.py --case cosine --prepare-only
python3 runSuite.py --case source --prepare-only
```

Generate or dispatch a single case:

```bash
cd cosine
python3 runTest --prepare-only
python3 runTest
```

`runTest` uses the repo's `mpm_vv_folder_runner.py`. Depending on your normal
PFW environment, it submits GEOS through Slurm and then schedules
post-processing after the simulation, or executes locally. It may require the
appropriate `GEOS_BANK`, PFW Python environment, executable path and partition
settings (the same as any other existing GEOS PFW verification). It does **not**
automatically turn a failure into a pass.

When you already have a simulation output directory:

```bash
python3 postProcess_cosine.py --run-dir /path/to/thermalVV_cosine \
   --output-dir /tmp/thermal-analysis --raise-on-failure
```

This writes `cosine_metrics.csv`, `cosine_summary.json`,
`cosine_reference.png` (if Matplotlib is installed), and `cosine_results.tex`.
The exit status is nonzero when output is missing, physics checks fail, or the
simulation did not reach its final time.

Run a selection:

```bash
python3 runSuite.py --case uniform --case cosine --case source
```

The full collection can be prepared with `python3 runSuite.py --prepare-only`.
The full case set may consume substantial HPC time; start with simple cases.

## Units and baseline properties

All inputs are in **SI**, with length `m`, time `s`, mass `kg`, thermal energy `J`,
power `W`, and temperature `K`.

| Quantity | Setting |
|---|---:|
| Domain length | `L = 0.010 m` |
| Cross-section height | `Y = 0.00125 m` |
| Mass density | `rho = 1000 kg/m³` |
| Constant heat capacity | `c = 10 J/(kg K)` |
| Isotropic conductivity | `k = 5 W/(m K)` |
| Thermal diffusivity | `alpha = 5e-4 m²/s` |
| Reference temperature | `T0 = 300 K` |
| Thermal difference | `DeltaT = 100 K` |
| Uniform specific heat input | `q_m = 20000 W/kg` |
| Global area-specific resistance | `R'' = 0.025 m² K/W` |

**Do not interpret these as typical physical material properties.** They were
chosen to make analytical thermal transients affordable at explicit MPM time
steps. PFW does not convert units: do not paste the values into your usual
mm/us/mg input without explicit conversion.

## Analytical verification matrix

| Folder | Verification / exact reference |
|---|---|
| `uniform` | Constant initial temperature, zero heat sources: `T(x,t)=T0`. |
| `cosine` | Periodic Fourier mode: `T=T0+(DeltaT/2)cos(2*pi*x/L)exp(-alpha*(2*pi/L)^2*t)`. |
| `cosine_coarse` | Same exact Fourier solution on a coarse mesh. |
| `cosine_fine` | Same exact Fourier solution on a fine mesh. |
| `gaussian` | Early-time infinite-domain Gaussian spreading: variance `sigma²+2*alpha*t`. |
| `step` | Early-time complementary-error-function diffusion across an initial temperature jump. |
| `source` | Uniform adiabatic heating: `T=T0+q_m*t/c`. |
| `sink` | Uniform adiabatic cooling: `T=T0-q_m*t/c`. |
| `source_local` | Conductivity disabled: only the selected particle region rises by `q_m*t/c`. |
| `capacity_linear` | Nonlinear `c(T)=c0+c1*(T-Tr)`; invert `c0*theta+0.5*c1*theta²=q_m*t`. |
| `conductivity_linear` | Steady nonlinear `k(T)=k0+k1*(T-Tr)`; constant Kirchhoff-transform slope. |
| `dirichlet_steady` | Linear 1D profile between finite hot/cold isothermal particle strips. |
| `dirichlet_transient` | Exact 1D Fourier sine series between the selected isothermal strips. |
| `contact_perfect` | Zero contact resistance tends to the continuous-medium error-function profile. |
| `contact_gap` | Finite open gap: no direct thermal contact transfer. |
| `contact_global` | Two finite-capacity thermal lumps, global contact resistance. |
| `contact_material` | Material-pair override of global resistance. |
| `contact_group` | Group-pair override of material-pair and global rules. |
| `contact_combined` | Combined material-and-group rule has highest lookup priority. |
| `contact_same_material` | Same material, separate contact groups, group-specific resistance. |
| `contact_same_group` | Distinct materials sharing a mechanical field, material-pair resistance. |
| `contact_insulating` | Explicit `Insulated` pair prevents contact exchange. |
| `contact_dfg` | Same-material damaged interface, thermal exchange across DFG partitions. |
| `contact_three_field` | Star graph A–C and B–C, A–B explicitly insulated; analytical 3-ODE solution. |
| `energy_closed` | Insulated, spatially varying temperature; total caloric energy remains constant. |
| `thermoelastic_free` | Free strain `epsilon=alpha_T*(T-T0)` with near-zero stress. |
| `plastic_adiabatic` | Model-derived `DeltaT=beta*integral(sigma:d_epsilon_p)/(rho*c)`; **known missing feature**, intentionally red until dissipation coupling is implemented. |
| `timestep` | High-conductivity Fourier decay; thermal substep stability under reduced thermal step. |
| `mpi_reference` | Serial Fourier run for MPI comparison. |
| `mpi_2rank` | Same initial state with two MPI ranks. |
| `gpu_reference` | Host reference for execution-path comparison. |
| `gpu_device` | HIP-enabled build comparison (existing extension currently performs thermal work on host). |

## Analytical limits and acceptance criteria

- `contact_global`, material/group overrides, DFG, and `contact_three_field`
  compare against exact **lumped** conduction ODE solutions, not the exact
  spatially resolved diffusion PDE. The Biot number is explicitly reported.
  Compare at multiple resistances/mesh refinements when gradients matter.
- The `gaussian`, `step`, and `contact_perfect` references are infinite-domain
  or semi-infinite-domain limits. They are compared only at early times and away
  from exterior boundaries.
- Isothermal boundaries are modeled as **finite particle strips**, not a
  zero-thickness mesh-face constraint; analytical coordinates reflect the strip
  boundaries. The transient finite-strip test needs mesh convergence around
  the bath/free-region interfaces.
- These are **numerical tolerances**, mostly normalized RMSE of the measured
  temperature histories against analytic reference. They are intentionally
  looser for coarse spatially resolved MPM. A PASS means a coarse V&V test passed
  its stated threshold, not that continuum convergence has been established.
- For constant heat capacity and fixed total mass, `c*(meanT-T0)` gives
  caloric specific energy. **Do not use `boxAverageHistory.internalEnergy` as
  transported thermal energy**: the current extension keeps that field as
  mechanical stress-work bookkeeping. The suite derives its caloric balance
  from mass-averaged temperature. If available, `MPMThermal` balance residuals
  are also examined in solver logs.
- Different material/thermal-contact rules are tested for precedence and
  symmetry; they do not independently establish correct reconstructed contact
  *area* for nonplanar surfaces. Add interface-area mesh studies as needed.
- The DFG test deliberately uses a fully damaged left block against an
  undamaged right block, explicit opposing surface flags, and the neighbor list.
  Inspect the actual field assignment in results: if the solver leaves the
  interface single-field, its intended DFG functionality was **not exercised**.
- The known extension **does not automatically convert plastic, damage, or
  friction work into heat**. `plastic_adiabatic` is therefore a *forward-looking*
  regression test and should currently fail (or be recorded as expected failure).
  Tracer stress/plastic-strain measurements are needed for its reference.
- These scripts have been syntax-checked, their PFW dictionaries structurally
  checked, and all post-processors exercised against synthetic analytical CSVs.
  **No end-to-end GEOS or MPI/HIP solver test was run here**; a GEOS executable
  and the current post-patch source checkout were unavailable.

## Post-processing and cross-run comparisons

For matching serial / two-rank result directories:

```bash
python3 thermal_compare.py --reference /path/to/mpi_reference \
  --candidate /path/to/mpi_2rank --output-dir compare_mpi
```

For host / HIP-enabled runs:

```bash
python3 thermal_compare.py --reference /path/to/gpu_reference \
  --candidate /path/to/gpu_device --output-dir compare_gpu
```

The comparator checks final particle temperature, time, and tracer coordinates.
It uses absolute and relative temperature tolerances and fails on missing
tracers or mismatched runs. It is an execution-path regression, not proof
that thermal kernels execute on HIP.

To require decreasing Fourier error under refinement:

```bash
python3 thermal_convergence.py \
  --cosine_coarse /path/to/cosine_coarse \
  --cosine /path/to/cosine \
  --cosine_fine /path/to/cosine_fine \
  --output fourier_convergence.json
```

To increase the number of x cells in any case, set
`MPM_THERMAL_NX=66` before PFW generation. The finer and coarser named
cases already choose their own default resolutions.

## Generated PFW conventions and portability

- `pfw_input_*.py` files import `thermal_pfw.build()` via an explicit
  parent-directory `sys.path` entry, and use dependency annotations.
- `pfw["thermal"]` is a structured dictionary with `enabled=1`, diagnostics,
  global contact mode/resistance, `contactResistance`, `heatSources`, and
  `isothermal` lists. It is serialized by the *patched* PFW formatter.
- Pair precedence is combined > group > material > global. The suite
  exercises all levels; `contact_three_field` excludes the A–B pair.
- Sources and isothermal particles are selected using constitutive
  `targetRegions`. They are supported by the patch and avoid relying on
  implementation-specific global particle IDs.
- `materialPropertyString` embeds named `ThermalElasticIsotropic`,
  `LinearThermalElasticIsotropic`, or `ThermalVonMisesJ` instances.
- This suite assumes that the relevant material generators support thermal
  fields, and the geometry code supports callable particle temperature.
- Solver log files and output directory naming vary across LLNL job runners;
  use `--run-dir` for the most reliable postprocessing.

## Developer checks

```bash
python3 -m compileall -q .
python3 -m unittest -q test_thermal_references
python3 -m unittest -q test_thermal_builders
```

`test_thermal_references` generates exact tracer CSVs in temporary directories
and tests the real postprocessing/acceptance code. That does **not** substitute
for running the corresponding GEOS test cases.
