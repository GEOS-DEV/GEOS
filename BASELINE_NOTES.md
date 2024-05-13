
Notes
==========

This file is designed to track changes to the integrated test baselines.
Any developer who updates the baseline ID in the .integrated_tests.yaml file is expected to create an entry in this file with the pull request number, date, and their justification for rebaselining.
These notes should be in reverse-chronological order, and use the following time format: (YYYY-MM-DD).

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
