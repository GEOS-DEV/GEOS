
Notes
==========

This file is designed to track changes to the integrated test baselines.
Any developer who updates the baseline ID in the .integrated_tests.yaml file is expected to create an entry in this file with the pull request number, date, and their justification for rebaselining.
These notes should be in reverse-chronological order, and use the following time format: (YYYY-MM-DD).

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
