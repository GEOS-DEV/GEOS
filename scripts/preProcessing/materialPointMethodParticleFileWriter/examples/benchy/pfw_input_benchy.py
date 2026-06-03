# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: STL import and prescribed multi-field contact demonstration.
# Solver/PFW features: Copper 3DBenchy boat and steel ball use separate particle regions/contact groups.
# Workflow note: keep this file as a copyable problem definition.  Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

"""3DBenchy FMPM contact example.

Purpose
-------
This input demonstrates three PFW/GEOS-MPM features in a small, copyable
example:

1. Importing a closed STL surface into PFW using ``geom.stl``.
2. Generating explicit surface normals and surface-position vectors from the STL
   surface so GEOS can use prescribed multi-field contact.
3. Assigning two bodies to different ``ContactGroup`` values: a copper 3DBenchy
   boat and a steel impactor ball.

The file is intentionally a normal PFW input, not a test-suite driver.  Users
should be able to copy it, change the STL, materials, grid, and loading, and run
it with ``particleFileWriter.py`` or an example ``runProblem`` wrapper.

Units
-----
The material database in ``pfw_materials.py`` uses the standard MPM validation
unit system: mm, microseconds, milligrams, and Kelvin.  Density is mg/mm^3 and
stress is GPa.  The 3DBenchy STL in this directory is already in millimeters, so
``scale=1.0`` is used.
"""

# PFW copies files listed with #[pfw_dependency] from userDefs_$USER.pfwPath into
# the run directory before importing this file.  The STL is then opened from the
# run directory, which makes the calculation independent of the source checkout.
#[pfw_dependency] pfw_materials.py
#[pfw_dependency] examples/benchy/3DBenchy.stl

import numpy as np

import pfw_geometryObjects as geom
import pfw_materials as matdb

# =============================================================================
# Problem size and run duration
# =============================================================================

pfw = {}
pfw["runDebug"] = True

# This is a demonstration input rather than a production-resolution Benchy run.
# Increase cpp and/or the partition counts for a better-resolved STL surface.
stopTime = 1.0

# Grid partitions.  The total MPI ranks requested by PFW are
# xpar*ypar*zpar.  The values below use 24 ranks, which is convenient for Dane's
# pdebug queue and keeps the STL import example quick enough for iteration.
pfw["xpar"] = 4
pfw["ypar"] = 2
pfw["zpar"] = 3

# Cells per partition in each direction.  Total cell counts are nI,nJ,nK below.
cpp = 12
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = pfw["zpar"] * cpp

# Candidate particles per cell direction before material refinement.  PFW then
# writes CPDI particles using the RVector field requested below.

# The old Benchy input used the STL's natural extents
# to the left for the incoming ball and a small free space around the boat.
pfw["xmin"] = -50.0
pfw["xmax"] = 40.0
pfw["ymin"] = -25.0
pfw["ymax"] = 25.0
pfw["zmin"] = 0.0
pfw["zmax"] = 60.0

# This is a fully 3D problem.  The Brazilian disk examples use planeStrain=1;
# Benchy keeps all three dimensions active.
pfw["planeStrain"] = 0

# Sorting lets PFW skip objects whose x-bounds do not overlap the current x-slice.
# This is useful here because the STL boat and ball occupy different x ranges.
pfw["sortObjects"] = True

# =============================================================================
# Batch settings used by particleFileWriter.py
# =============================================================================

# Bank/account, GEOS executable, default run directory, and default Python command
# are read from userDefs_$USER.py.
pfw["mWallTime"] = "00:10:00"
pfw["mSubmitJobs"] = True

# Keep this example simple: no automatic restart machinery.

# =============================================================================
# GEOS MPM solver controls
# =============================================================================

# Plot only the initial and final states for a compact demonstration output.  The
# restart interval is beyond endTime, so no restart files are written.
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime
pfw["restartInterval"] = 2.0 * stopTime
pfw["outputType"] = "silo"

# Explicit dynamic MPM settings.
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1e-16

# FMPM transfer with CPDI particles.  updateOrder=2 is the typical second-order
# setting used by the compact example problems.
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 2
pfw["cpdiDomainScaling"] = 1

# No damage-field partitioning/DFG for this example.  The goal is STL import and
# multi-field contact, not ceramic damage evolution.
pfw["damageFieldPartitioning"] = 0

# Explicit contact settings.  SurfacePosition and SurfaceNormal are generated
# below from geom.stl and geom.sphere.  ContactGroup separates the boat and
# the ball into different material-point fields for prescribed multi-field
# contact.  No neighbor-list file is needed for this example.
pfw["needsNeighborList"] = 0
pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"] = 1
pfw["explicitSurfaceNormalInfluence"] = 1000.0
pfw["frictionCoefficient"] = 0.25

# Basic diagnostics and guardrails.
pfw["solverProfiling"] = 1
pfw["reactionHistory"] = 1
pfw["boxAverageHistory"] = 1
pfw["maxParticleVelocity"] = 1.0e4
pfw["minParticleJacobian"] = 0.01
pfw["maxParticleJacobian"] = 10.0

# Particle fields written into the initial particle file.  The contact-specific
# fields are ContactGroup, SurfaceFlag, SurfaceNormal, and SurfacePosition.
# RVector stores the CPDI particle-domain vectors.
pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "SurfaceFlag",
    "SurfaceNormal",
    "SurfacePosition",
    "RVector",
]

# =============================================================================
# Materials
# =============================================================================

# Material IDs are assigned by list order.  The copper boat is material 0 and the
# steel ball is material 1.  These dictionaries come directly from pfw_materials.
pfw["materials"] = [matdb.copper["name"], matdb.steel["name"]]
pfw["materialPropertyString"] = matdb.copper["materialString"] + "\n" + matdb.steel["materialString"]

COPPER_BOAT = 0
STEEL_BALL = 1

# =============================================================================
# Geometry
# =============================================================================

# The STL class reads a closed triangulated surface, builds a binned ray-casting
# inside/outside test, and returns nearest-triangle surface normals/positions for
# surface-flagged particles.  The Benchy STL is in millimeters and is already
# positioned with the bottom at z=0, so no scale or translation is needed.
boat = geom.stl(
    "copper_3DBenchy",
    fileName = "3DBenchy.stl",
    scale = 1.0,
    x0 = [0.0, 0.0, 0.0],
    vel = [0.0, 0.0, 0.0],
    mat = COPPER_BOAT,
    group = 0,
    particleType = 2,       # CPDI particles
    rayAxis = 0,            # +x ray through yz bins for inside/outside tests
    binCounts = (96, 96),   # Larger values reduce triangles checked per ray
    kNearest = 64,          # Triangles checked for nearest surface projection
)

# A steel ball starts to the left of the STL and moves in +x toward the boat.
# It uses a different ContactGroup so GEOS can exercise multi-field contact.
ball = geom.sphere(
    "steel_ball",
    [-42.0, 0.0, 26.0],
    5.0,
    vel = [200.0, 0.0, 0.0],
    mat = STEEL_BALL,
    group = 1,
    particleType = 2,
)

# Object order matters for overlapping geometry: the first matching object wins.
# The ball starts outside the boat, but listing it first makes the intended
# material/contact group unambiguous if a user moves it into overlap for testing.
pfw["objects"] = [ball, boat]

# =============================================================================
# Boundary conditions
# =============================================================================

# No prescribed external loading is applied to the grid boundaries in this small
# example.  The ball motion is prescribed through its initial particle velocity.
# Current GEOS expects fTableInterpType as one of Linear, Cosine, or Smoothstep;
# Linear is used here even though no active prescribed f-table is used.
pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0
pfw["fTableInterpType"] = "Linear"

# Boundary-condition order is x-, x+, y-, y+, z-, z+.  Symmetry walls keep the
# compact demonstration contained.  For an open impact calculation, replace some
# faces with Outflow according to the MPM solver boundary-condition enum used in
# other PFW examples.
pfw["boundaryConditionTypes"] = [1, 1, 1, 1, 1, 1]
