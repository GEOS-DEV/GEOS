Notes
==========

This file is designed to track changes to the integrated test baselines.
Any developer who updates the baseline ID in the .integrated_tests.yaml file is expected to create an entry in this file with the pull request number, date, and their justification for rebaselining.
These notes should be in reverse-chronological order, and use the following time format: (YYYY-MM-DD).

PR #3777 (2026-08-16) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3777-17413-0a5f007.tar.gz>
Move CO2 Brine parameters to xml

PR #4127 (2026-08-23) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4127-17399-63445db.tar.gz>
Rebaseline five restart checks after the TPL update changed VTK/Scotch mesh partitioning. Global mesh topology and fields are unchanged.

PR #3884 (2026-08-16) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3884-17320-39debdc.tar.gz>
=====================
Total stress fix in the thermo-poromechanics model

PR #4114 (2026-08-14) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4114-17283-3bb33cb.tar.gz>
=====================
Stop dumping linear systems from ATS decks (`writeLinearSystem` no longer set). That flag is stored in restart files, so `perf_status_test` restartchecks need a new baseline.

PR #4088 (2026-07-27) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4088-17188-4247846.tar.gz>
=====================
Fluid reset after convergence failure

PR #3972 (2026-07-28) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3972-17154-316e6d8.tar.gz>
=====================
Well model refactor .  Integrated test update due to schema changes

PR #3836 (2026-05-20) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3836-17046-2e89f64.tar.gz>
=====================
Added statistics `Group` objects for each statistics `Task` instance

PR #4040 (2026-06-16) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4040-16993-1393f80.tar.gz>
=====================
Move relperm driver to use new constitutive driver framework

PR #3705 (2026-06-12) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3705-16862-2b262bf.tar.gz>
=====================
Implement compositional enthalpy model

PR #4074 (2026-06-10) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4074-16937-bf66240.tar.gz>
=====================
Change triaxial driver to use restart for checks

PR #4067 (2026-06-10) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4067-16928-3038f6b.tar.gz>
=====================
Add Coulomb friction/cohesion input from vtk mesh

PR #4068 (2026-06-09) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4068-16828-c74157c.tar.gz>
=====================
Add MPI runs for smoke tests with surfaceGenerator

PR #4062 (2026-05-26) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4062-16784-6d8782e.tar.gz>
=====================
Add Porous Solid other than PorousElasticity for ALM solver

PR #4057 (2026-05-21) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4057-16739-5dde641.tar.gz>
=====================
Remove dependency on PVT package

PR #3814 (2026-05-20) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3814-16717-d17ea59.tar.gz>
=====================
Fix RLF coloring

PR #4008 (2026-05-18) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4008-16688-17c55fe.tar.gz>
=====================
Fix Fracture/3D cell co-location in parallel mesh redistribution

PR #4055 (2026-05-18) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4055-16683-0251bb6.tar.gz>
=====================
Trim some fluid model tests

PR #4041 (2026-05-16) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4041-16664-9a29348.tar.gz>
=====================
Fix wellbore nonlinear thermal diffusion

PR #3977 (2026-05-15) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3977-16644-9fc03e7.tar.gz>
=====================
Change face normal, centers and area for local Newell's formula

PR #3999 (2026-05-14) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3999-16602-542583d.tar.gz>
=====================
Fault Perm Update for Contact solvers

PR #4029 (2026-05-06) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4029-16507-7e3d8b5.tar.gz>
=====================
Add single phase viscosity dependency on temperature

PR #3959 (2026-05-04) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3959-16478-faf1698.tar.gz>
=====================
Add reference thermal conductivity

PR #4021 (2026-04-14) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4021-16339-bb862da.tar.gz>
=====================
Add Young Modulus & Poisson import from VTK mesh

PR #3883 (2026-04-10) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3883-16299-3037085.tar.gz>
=====================
Move PVT Driver tests from unit tests to integrated tests

PR #4007 (2026-04-03) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr4007-16213-67a3002.tar.gz>
=====================
Add XML input parameter: "hypredriveInputFile"

PR #3957 (2026-03-30) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3957-16171-da51804.tar.gz>
=====================
Add checkEulerCharacteristic option, rebaseline due to new input.

PR #3967 (2026-03-27) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3967-16106-c0f34de.tar.gz>
=====================
Fix 2D/3D cell co-location in parallel mesh redistribution

PR #3986 (2026-02-27) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3986-15734-7487221.tar.gz>
=====================
Corrected traction boundary conditions

PR #3970 (2026-02-11) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3970-15479-074f42a.tar.gz>
=====================
Bypass well residual calculation for closed wells

PR #3964 (2026-02-09) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3964-15460-26718eb.tar.gz>
=====================
Fix fracture state update for ALM solver

PR #3940 (2026-01-27) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3940-15307-53de7ba.tar.gz>
=====================
Fix the transimissibility calculated between a cell and a surface element

PR #3926 (2026-01-26) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3926-15289-1f676c0.tar.gz>
=====================
MultiPhase Poromechanics ALM solver and a test with curved fractures

PR #3634 (2025-12-31) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3634-15105-6ab70ef.tar.gz>
=====================
Add singlephase reactive transport solver integrated with HPCReact

PR #3795 (2025-12-19) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3795-15047-606f4ac.tar.gz>
=====================
Addition of multiphase initialisation

PR #3925 (2025-12-18) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3925-15032-d307303.tar.gz>
=====================
Add traction update for ALM solver and a test with curved fractures

PR #3829 (2025-11-27) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3829-14878-fdcd94e.tar.gz>
=====================
Fix validation of average region stat needed by option

PR #3522 (2025-11-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3522-14736-f07298f.tar.gz>
=====================
Single Phase Poromechanics Conforming Fractures

PR #3765 (2025-10-30) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3765-14543-498692f.tar.gz>
=====================
Adding boundary fields to handle Dirichlet conditions in compositional MFD.

PR #3849 (2025-10-23) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3849-14514-aaaf0f9.tar.gz>
=====================
Add multiphase contact with wells

PR #3880 (2025-10-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3880-14441-1132122.tar.gz>
=====================
Fix a bug introduced in #3485: mass that is used in accumulation term was not updated with porosity change after mechanics leading to always converged sequential outer loop.

PR #3299 (2025-10-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3299-14426-6c93a0d.tar.gz>
=====================
Add co2 injection into gas with SW EoS and k-value flash.

PR #3279 (2025-10-24) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3279-14414-db55426.tar.gz>
=====================
Output cell-wise average of each stress and strain component.

PR #3851 (2025-10-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3851-14171-9d950d6.tar.gz>
=====================
Enable `FullyImplicit` for `SinglePhaseReservoirPoromechanicsConformingFractures`.

PR #3808 (2025-10-10) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3808-14160-03c20b2.tar.gz>
=====================
Add errorSetMode for the FiedSpecification.

PR #3193 (2025-10-10) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3193-14118-8ee1c34.tar.gz>
=====================
Enable geothermal gradient in HydrostaticEquilibrium for single-phase flow.

PR #3359 (2025-10-07) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3359-14024-3c9ebb4.tar.gz>
=====================
Add functions to connect well perforation to surface elements.

PR #3843 (2025-10-02) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3843-13975-76d125d.tar.gz>
=====================
Fix linear solver issues.

PR #3842 (2025-10-01) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3780-13967-17d92c8.tar.gz>
=====================
Refactor of the single-phase hybrid MFD implementation to remove the upwinding scheme from the discretization.

PR #3842 (2025-10-01) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3842-13948-6f07e39.tar.gz>
=====================
Enable parallel run for some contact tests.

PR #3838 (2025-09-30) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3838-13915-07c068d.tar.gz>
=====================
Fix statistics update logic for coupled solvers.

PR #3790 (2025-09-29) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3790-13899-bb7b286.tar.gz>
=====================
Enable Kozeny-Carman Permeability for PorousSolid.

PR #3779 (2025-09-21) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3779-13734-44eed3f.tar.gz>
=====================
Add new inputs for function input var scaling. Add new Multiscale linear solver parameters XML block.

PR #3821 (2025-09-19) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3821-13726-7dd8089.tar.gz>
=====================
Enable parallel versions for some contact mechanics tests.

PR #3796 (2025-09-18) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3796-13696-bfe23a1.tar.gz>
=====================
Add solid fields.

PR #3801 (2025-09-17) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3801-13669-3ca52ce.tar.gz>
=====================
Tolerance for geometric objects coordinates check.

PR #3813 (2025-09-16) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3813-13627-4228acd.tar.gz>
=====================
Create a separator ("fluid model") for each well. Only schema differences in results.

PR #3629 (2025-09-16) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3711-13370-49f4348.tar.gz>
=====================
Fix some bugs in surface generation communication.

PR #3745 (2025-09-14) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3745-13577-97fabfe.tar.gz>
=====================
Oscillation detection and scaling option.

PR #3776 (2025-09-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3776-13560-88bd98b.tar.gz>
=====================
Constitutive cleanup: rebaseline due to technical diffs, no real results changes.

PR #3349 (2025-09-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3349-13555-41d98d4.tar.gz>
=====================
Configuration loop acceleration for `SolidMechanicsLagrangeContact`.

PR #3812 (2025-09-11) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3812-13496-8dc0bc6.tar.gz>
=====================
Add reset of k-values for compositional fluid model. Improves the convergence of models at the start.

PR #3285 (2025-09-09) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3285-13414-69e5962.tar.gz>
=====================
Add hydrofrac verification cases for leak-off.

PR #3587 (2025-09-08) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3587-13389-99ac8e4.tar.gz>
=====================
Perforation status option. Updates for schema changes and well quantities not being compute if well is closed
PR #3788 (2025-09-07) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3788-13372-e0a1d67.tar.gz>
=====================
Updating txt files for class09_pb3_hystRelperm.

PR #3621 (2025-09-05) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3621-13365-2c193d7.tar.gz>
=====================
Reservoir volume well constraint option.

PR #3629 (2025-09-02) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3680-13304-f49e13a.tar.gz>
=====================
Logic for deciding whether to setup linear solver, new flag `reuseFactorization`.

PR #3629 (2025-08-29) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3629-13262-9e92109.tar.gz>
=====================
Add solver statistics wrapper.

PR #3627 (2025-08-29) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3672-13218-b729958.tar.gz>
=====================
Some fields were not being syncronized as part of the parallel topology change. This PR syncs them and produces a new baseline as a result.

PR #3755 (2025-08-25) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3755-13107-9203370.tar.gz>
=====================
LogInfo cleanup.

PR #3783 (2025-08-25) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3783-13079-916d2a0.tar.gz>
=====================
Update bug in single phase flash handling.

PR #3224 (2025-08-22) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3224-13023-0334b63.tar.gz>
=====================
Add Taper boundary conditions inside second-order wave solvers.

PR #3781 (2025-08-20) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3781-12988-f27ea4b.tar.gz>
=====================
Remove unused from subregions.

PR #2207 (2025-08-19) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr2207-12962-f0dbaad.tar.gz>
=====================
Factoring hysteresis model out of `TableRelativePermeabilityHysteresis`.

PR #2427 (2025-08-09) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr2427-12699-941dbff.tar.gz>
=====================
Change default value of amgNumFunctions from 1 to 3 for solid mechanics solvers. Change in mesh partitioning of PoroElastic_hybridHexPrism_co2 cases due to Scotch version update.

PR #3682 (2025-08-08) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3682-12632-86ad358.tar.gz>
=====================
Add physics-based scaling option.

PR #3662 (2025-08-07) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3662-12617-6b694b1.tar.gz>
=====================
Change to single phase handling of flash for compositional fluid model.

PR #3622 (2025-08-07) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3622-12602-d1632ed.tar.gz>
=====================
Moved "parallelThread" from OutputBase to SiloOutput as it is not useful on any other sub-class

PR #3748 (2025-08-06) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3556-12451-d672aab.tar.gz>
=====================
Use old pore volume in volume balance equation - minor diffs for compositional flow tests.

PR #3746 (2025-08-01) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3556-12451-d672aab.tar.gz>
=====================
Fix initial composition for `2ph_cap_1d_ihu`.

PR #3556 (2025-07-31) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3556-12451-d672aab.tar.gz>
=====================
Enable BartonBandis model add new smoke tests.

PR #3568 (2025-07-31) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3568-12439-ac82cb0.tar.gz>
=====================
Baselines updated due to set reference state for temperature and add new smoke tests.

PR #3740 (2025-07-25) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3740-12360-2177cb4.tar.gz>
=====================
Add missing hydraulic aperture update for sequential poromechanics with conforming fractures.

PR #3732 (2025-07-18) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3732-12212-a92d996.tar.gz>
=====================
Add `numTimestepsSinceLastDtCut` to restart.

PR #3730 (2025-07-18) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3730-12201-ebfc7cf.tar.gz>
=====================
Add 3 tests for compositional Soreide-Whitson EOS.

PR #3659 (2025-07-17) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3517-12189-54b7075.tar.gz>
=====================
Fields and constitutives refactor.

PR #3659 (2025-07-08) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3659-12039-3662bde.tar.gz>
=====================
Add thermal to single-phase well.  Baselines updated due to schema changes.

PR #3635 (2025-06-11) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3635-11765-c0e7e87.tar.gz>
=====================
Add new wave solver (elastic anisotropic TTI).

PR #3679 (2025-05-27) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3679-11653-e066fbb.tar.gz>
=====================
Removed `maxStableDt` and `registerWrapper` for `meshTargets`.

PR #3653 (2025-05-13) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3653-11335-b8096ce.tar.gz>
=====================
Change black oil phase labelling for gas only cells.

PR #3274 (2025-05-08) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3274-11275-3cb35d1.tar.gz>
=====================
New flag `allowNonConvergedLinearSolverSolution` for solvers.

PR #3524 (2025-05-02) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3524-11210-f1b043a.tar.gz>
=====================
Immiscible multiphase flow.

PR #3626 (2025-04-28) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3626-11189-dfa74ce.tar.gz>
=====================
Update in VTK caused change in partitioning for reading vtk meshes. Verified baselines using new scripts/parallelRestartDiff.py

PR #3624 (2025-04-15) <https://storage.googleapis.com/geosx/integratedTests/baseline_integratedTests-pr3624-11053-ae011c7.tar.gz>
=====================
Bugfix for time step cut in sequential, minor time step logic change when a cut occurs.

PR #3537 (2025-04-02)
=====================
Added two attributes to TableFunction: writeCSV and logLevel.

PR #3589 (2024-03-26)
=====================
Hypre update - make co2 cases use direct solver.

PR #3396 (2024-03-21)
=====================
Use solid mechanics solver directly to perform poromechanics initialization.

PR #2125 (2024-03-20)
=====================
Phase-field nucleation model.

PR #3494 (2024-03-11)
=====================
Added more log level documentation

PR #3427 (2024-03-10)
=====================
Well time step selector based on rates/bhp tables and clarify well rates logic.

PR #3485 (2024-03-09)
=====================
Use mass and energy consistently for single phase solvers.

PR #3460 (2024-03-08)
=====================
Refactor single phase constitutive containers.

PR #3525 (2025-03-06)
=====================
Add analytical leakoff feature for hydrofrac solver.

PR #3401 (2025-03-05)
=====================
Bugfix for IHU.

PR #3483 (2025-03-02)
=====================
Remove relative permeability from wells.

PR #3576 (2025-03-01)
=====================
Add an option to skip density and viscosity computes when phase is not present for CO2 fluid update.

PR #3571 (2025-02-28)
=====================
Do not allow negative pressure by default, except for hydrofrac.

PR #3551 (2025-02-19)
=====================
Add Passing Crack to the integrated tests.

PR #3541 (2025-02-18)
=====================
Well control parallel synchronization fix.

PR #3443 (2025-02-17)
=====================
Added tests for overall composition (Z) formulation.

PR #3547 (2025-02-17)
=====================
Multiphase contact bugfix, add test case to ats, removed redundant linear solver params for other tests.

PR #3546 (2025-02-15)
=====================
Fix 1d edfm case and add it to ats.

PR #2968 (2025-02-13)
=====================
Replace array1d<string> with std::vector<string>.

PR #3227 (2025-02-06)
=====================
Add targetRegion for perforations (optional).

PR #3502 (2025-02-04)
=====================
Add array to store the source values in time inside wave solvers.

PR #3395 (2025-01-22)
=====================
Add new fields and change the default input for some tests.

PR #3416 (2025-01-21)
=====================
Refactoring of induced seismicity EQ solvers to add coupling.

PR #3310 (2025-01-21)
======================
Scalable rock toughness required new field.

PR #3228 (2025-01-15)
=====================
deltaVolume added in multiphase.

PR #3495 (2025-01-08)
=====================
Add missing logic to support switching from fixed mass rate injection rate constraint to max injection pressure.

PR #3384 (2025-01-07)
=====================
Added plastic strain output.

PR #3486 (2025-01-06)
=====================
useNewGravity became gravityDensityScheme.

PR #3479 (2024-12-15)
=====================
Refine inputFiles/compositionalMultiphaseFlow: shift reference pressures to initial pressures, make nonlinear tuning more reasonable, minimize output.

PR #3450 (2024-12-14)
=====================
Fix timestep selector flaw in SolidMechanicsLagrangeContact.

PR #3450 (2024-12-08)
=====================
Added test for explicit runge kutta sprinslider.

PR #3480 (2024-12-06)
=====================
Add "logLevel" parameter under /Problem/Outputs in baseline files

PR #3361 (2024-12-03)
=====================
Revert default gravity treatment to old version. Make the way introduced in #3337 optional.

PR #3361 (2024-12-03)
=====================
Baseline diffs after reimplementation of wave equation acoustic gradient for velocity and density parameters: new field "partialGradient2" and "pressureForward" field replacing "pressureDoubleDerivative".

PR #3393 (2024-12-02)
=====================
Fix netToGross bug.

PR #3381 (2024-12-01)
=====================
A few baseline diffs for order FaceElementSubRegion::m_toFacesRelation map. Not sure why this was changed by this PR, but the previous order seems incorrect for a couple of cases.

PR #2957 (2024-11-27)
=====================
Added ExternalDataRepository.

PR #3448 (2024-11-21)
=====================
Switched the FaceElementSubRegion::m_toFacesRelation and FaceElementSubRegion::m_2dElemToElems back to array2d instead of ArrayOfArray. This results in a reordering m_toFacesRelation back to the "correct" assumed order of "original face first". This fixes a bug that failed to remove the CellStencil entry when a FaceElement splits two cells.

PR #2637 (2024-11-21)
=====================
Added numberOfTargetProcesses.

PR #3439 (2024-11-20)
=====================
EDFM bugfixes: derivatives sign, frac/cell element volume, fix apertures inconsistency in test cases.

PR ##3440 (2024-11-18)
=====================
Added Lagrange multiplier with bubble functions stabilization (sli only) and possibility to specify a slip.

PR #3339 (2024-11-14)
=====================
Hypre improvements, rebaseline is due to field value change (amgNumFunctions).

PR #3434 (2024-11-09)
=====================
Bugfix: Fixed output of ArrayOfArray objects to restart files.

PR #3374 (2024-11-09)
====================
Bugfix for gravity treatment in flux for thermal.

PR #3372 (2024-11-09)
====================
Fix a bug related to mass and energy updates for poromechanics.

PR #3426 (2024-11-08)
====================
Bugfix: reset accumulation in fracture when time step cut occurs in hydrofrac solver.

PR #3413 (2024-11-07)
====================
Add tests for poro-thermo-plastic model.

PR #3337 (2024-11-06)
====================
Change density treatment for gravity in multiphase flow solver.

PR #3408 (2024-11-06)
====================
EFEM bugfixes: effective traction + oldStress.

PR #3280 (2024-11-05)
====================
Added Sprig-slider test.

PR #2909 (2024-10-30)
=====================
Add routine for automatic time steps in waveSolvers with new attributes

PR #3156 (2024-10-29)
====================
Restart check errors due to 1) schema node added to enable thermal option in well model and 2) arrays removed/added for option.  Max difference errors due treatment of shutin wells.  Previously non-zero rate value reported for shutin well, new code will set rate arrays to zero.

PR #2878 (2024-10-17)
=====================
Sorted region cellBlocks names alphabetically. Therefore affected ordering of: faceManager/elemSubRegionList, nodeManager/elemList, nodeManager/elemSubRegionList, SurfaceElementSubRegion::fractureElementsToCellSubRegions, field::perforation::reservoirElementSubregion.

PR #3364( 2024-10-15)
=====================
Enable reservoir+wells+contact mechanics. Rebaseline needed because of 'allowNegativePressure' flag added for wells.

PR #3364( 2024-10-01)
=====================
Separate mass and volume residuals for output in compositional flow solver. Baseline update because of minor numerical diffs.

PR #3149( 2024-09-30)
=====================
Added new field "writeCSV"

PR #3163 (2024-09-20)
=====================
Added new fields (krylovStrongestTol, adaptiveGamma, adaptiveExponent) to the LinearSolverParameters for adaptive tolerances.

PR #3338 (2024-09-19)
======================
Updated time-stepping logic. Rebaseline due to new input parameter and minor numerical diffs.

PR #3217 (2024-09-16)
======================
ALM slip and open modes with relative tests.

PR #3318 (2024-09-12)
======================
Modified SeismicityRate poroelastic case.

PR #3322 (2024-09-06)
======================
Print out fracture state for contact model. Rebaseline the corresponding cases.

PR #3302 (2024-09-05)
======================
Added restartcheks to hydrofrac cases and reduced time of cases that were too long.

PR #3135 (2024-09-04)
======================
Temperature dependent single phase thermal conductivity. Rebaseline all thermal cases.

PR #3294 (2024-09-01)
======================
Re-enable enforcement of wave propagation integrated test pass.

PR #3300 (2024-08-28)
======================
Re-enable floating point exceptions. Rebaseline due to minor changing default value of maxRelativeCompDensChange from 1.7976931348623157e+308 to 1.7976931348623157e+208.

PR #3283 (2024-08-22)
======================
Reuse computeSinglePhaseFlux. Rebaseline due to minor numerical diffs.

PR #3249 (2024-08-14)
======================
Two initialization options for poromechanical models. Rebaseline the corresponding cases.

PR #3278 (2024-08-12)
======================
Renamed GEOSX to GEOS in enternal mesh import, so rebaseline to change these names is the baselines.

PR #3202 (2024-08-03)
======================
Acoustic VTI tests needed rebaselining after update in source and receiver location algorithm.

PR #3215 (2024-07-23)
======================
Changed the default value for massCreation and name of the wrapper.

PR #3194 (2024-07-22)
======================
Check pore volume for all element types, also check that default aperture > 0. Rebaseline for modified tests. No real results change.

PR #3213 (2024-07-12)
======================
Added baselines for new tests on Dirichlet boundary conditions for multiphase flow.

PR #3194 (2024-07-10)
======================
Use aperture table in poromechanics with conforming fractures. Rebaseline the corresponding cases.

PR #3006 (2024-07-01)
======================
Added baselines for new tests. Relaxing tolerances for singlePhasePoromechanics_FaultModel_smoke.

PR #3196 (2024-06-28)
======================
Added isLaggingFractureStencilWeightsUpdate to hydrofracture solve. Rebaseline because of the new input.

PR #3177 (2024-06-28)
======================
Added logLevel to TimeHistoryOutput. Rebaseline because of the new input flag.

PR #3181 (2024-06-25)
======================
Decouple debug matrix output from logLevel. Rebaseline because of the new input flag.

PR #3142 (2024-06-20)
======================
Adding output of total strain. Rebaseline because of new inclusion of strain in output.

PR #3170 (2024-06-19)
======================
Fix tutorial example for thermal debonding wellbore problem. Test case modified.

PR #3130 (2024-06-19)
======================
New solver for contact mechanics based on the Augmented Lagrangian Method (ALM). New test case added.

PR #3160 (2024-06-18)
======================
Two experimental options for compositional flow solver. Rebaseline because of the new input flags.

PR #3165 (2024-06-18)
======================
Small bug fix. Rebaseline required due to appearance of useTotalMassEquation in well solver params. No real results change.

PR #3088 (2024-06-17)
======================
Adding temperature-dependent Solid Volumetric Heat Capacity. Rebaseline because of the parameter change in SolidInternalEnergy.

PR #3100 (2024-06-14)
======================
Adding pressure stabilization for single phase poromechanics.

PR #3133 (2024-06-14)
======================
Fix node ordering for faceElements.

PR #3021 (2024-06-13)
======================
Preparatory work for fractures + wells. New test case added.

PR #3152 (2024-06-13)
======================
Some random things. Baseline update because of the new parameter (minScalingFactor).

PR #3138 (2024-06-11)
======================
Properly sync nonlinear solver params for coupled solver. Baseline update mostly due to number of iterations change in baseline files.

PR #3140 (2024-06-11)
======================
Fixed derivative in EzrokhiBrineDensity

PR #3080 (2024-06-07)
=====================
Rebaseline after adding viscoelastic wave propagator.

PR #3075 (2024-06-05)
=====================
Introduce configuration tolerance. Rebaseline because of the new parameter in NonlinearSolverParameters.

PR #3120 (2024-06-05)
=====================
Add missing compositionalMultiphaseFlow tests into ATS and adjust output naming. Rebaseline accordingly.

PR #3113 (2024-06-05)
=====================
Add general version updateConfiguration. Rebaseline of some edfm cases is needed.

PR #3050 (2024-05-20)
=====================
Spatially varying grain bulk modulus. Rebaseline of all poromechanics cases needed.

PR #3141 (2024-05-28)
=====================
Test cashing baselines locally.

PR #3125 (2024-05-16)
=====================
Remove field to store pressure gradient cell-wise for solvers that don't need it.

PR #2110 (2024-05-13)
=====================
new field to store pressure gradient cell-wise.

PR #3060 (2024-05-13)
======================
Rebaselined after addition of elastic VTI wave propagator.

PR #4950 (2024-05-10)
======================
Added smoke tests for SeismicityRate solver in inducedSeismicity.

PR #3086 (2024-05-09)
======================
Added a presure-dependent permeability model and the transmissibility calculation in the CellElementStencil

PR #3105 (2024-05-08)
======================
Added missing derivative for temperature, hence small numerical diffs in thermal tests results and numeracal behavior

PR #2917 (2024-05-07)
======================
New fields for wellsControls: wellControls1_ConstantMassRate_table, targetMassRate, massDensity, ...

PR #3044 (2024-05-02)
======================
Removed old integratedTests submodule
Implemented new baseline storage
Implemented new CI integrated tests
