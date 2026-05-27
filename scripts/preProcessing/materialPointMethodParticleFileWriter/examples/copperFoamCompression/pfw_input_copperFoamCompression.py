# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Dynamic reverse-ballistic uniaxial compression of 3D copper foam.
# Solver/PFW features: Periodic lateral boundaries, Poisson-disk spherical pores, SinglePointBSpline particles, and copper VonMisesJ plasticity.
# Workflow note: keep this file as a copyable problem definition. Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

"""3D copper foam reverse-ballistic compression.

The foam initially fills a 1:1:6 domain.  The long direction is z, and all
material points start with velocity ``[0, 0, -impactVelocity]`` so the specimen
impacts the symmetry plane at ``z=0``.  The x and y faces are periodic, and the
Poisson-disk pore field is periodic in x and y.

Units are the standard MPM example units: mm, microsecond, milligram, GPa, K.
"""

#[pfw_dependency] pfw_materials.py

import math
import numpy as np

import pfw_geometryObjects as geom
import pfw_materials as matdb

pfw = {}
pfw["runDebug"] = True

# =============================================================================
# User-facing problem controls
# =============================================================================

# Initial aspect ratio is 6:1, longer in the impact direction.
lateralLength = 1.0          # mm
impactAspectRatio = 6.0
domainSize = {
    "x": lateralLength,
    "y": lateralLength,
    "z": impactAspectRatio * lateralLength,
}

# Spherical pore controls.  The generated pore count is integer, so the geometry
# object prints the realized porosity when PFW imports this input.
poreDiameter = 0.15          # mm
poreVolumeFraction = 0.18    # target void volume fraction
poreMinLigament = 0.01       # mm additional clearance between pore surfaces
poreSeed = 20260520

# Reverse-ballistic impact speed.  Positive value here means motion toward the
# negative z symmetry plane.
impactVelocity = 0.25        # mm/us = 250 m/s

# Resolution controls.  cpp, xpar, ypar, and zpar are all scaled by refine.  The
# z partition count is six times the lateral count so the initial cells are
# approximately cubic for the 6:1 specimen.
refine = 1
cpp = 8
xpar = 3
ypar = 3
zpar = 18
ppc = 2                      # particles per cell per direction

# End time estimate: snowplow compaction wave speed in a porous bar.  The compacted
# material is assumed to lock at the matrix density unless lockedPorosity is raised.
lockedPorosity = 0.0
initialRelativeDensity = 1.0 - poreVolumeFraction
lockedRelativeDensity = 1.0 - lockedPorosity
compactionStrain = max(1.0e-12, 1.0 - initialRelativeDensity / lockedRelativeDensity)
compactionWaveSpeed = impactVelocity / compactionStrain
compactionWaveTransitTime = domainSize["z"] / compactionWaveSpeed
stopTime = compactionWaveTransitTime

# =============================================================================
# Domain and grid
# =============================================================================

pfw["xmin"] = -0.5 * domainSize["x"]
pfw["xmax"] =  0.5 * domainSize["x"]
pfw["ymin"] = -0.5 * domainSize["y"]
pfw["ymax"] =  0.5 * domainSize["y"]
pfw["zmin"] = 0.0
pfw["zmax"] = domainSize["z"]

pfw["periodic"] = [True, True, False]
pfw["planeStrain"] = 0

cppRefined = max(1, int(round(refine * cpp)))
pfw["xpar"] = max(3 if pfw["periodic"][0] else 1, int(round(refine * xpar)))
pfw["ypar"] = max(3 if pfw["periodic"][1] else 1, int(round(refine * ypar)))
pfw["zpar"] = max(3 if pfw["periodic"][2] else 1, int(round(refine * zpar)))

# nK includes two ghost cells because z is nonperiodic.  x and y are periodic and
# therefore do not include ghost cells.
pfw["nI"] = pfw["xpar"] * cppRefined
pfw["nJ"] = pfw["ypar"] * cppRefined
pfw["nK"] = pfw["zpar"] * cppRefined + 2
pfw["ppc"] = ppc

pfw["sortObjects"] = True
pfw["outputType"] = "silo"

# =============================================================================
# Batch settings used by particleFileWriter.py
# =============================================================================

pfw["mWallTime"] = "00:20:00"
pfw["mCores"] = pfw["xpar"] * pfw["ypar"] * pfw["zpar"]
pfw["mNodes"] = int(math.ceil(float(pfw["mCores"]) / 112.0))
pfw["mSubmitJobs"] = True

# =============================================================================
# GEOS MPM solver controls
# =============================================================================

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime
pfw["restartInterval"] = 2.0 * stopTime

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-16
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 2

# SinglePointBSpline particles do not use CPDI domain scaling.
pfw["cpdiDomainScaling"] = 0
pfw["damageFieldPartitioning"] = 0
pfw["needsNeighborList"] = 0

pfw["solverProfiling"] = 1
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stopTime / 1000.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime / 1000.0

pfw["maxParticleVelocity"] = max(10.0, 10.0 * impactVelocity)
pfw["minParticleJacobian"] = 0.01
pfw["maxParticleJacobian"] = 10.0

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

pfw["materials"] = [matdb.copper["name"]]
pfw["materialPropertyString"] = matdb.copper["materialString"]
COPPER = 0

# =============================================================================
# Geometry
# =============================================================================

foam = geom.poissonDiskFoam(
    "periodic_copper_foam_3d",
    x0 = [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    x1 = [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    poreDiameter = poreDiameter,
    porosity = poreVolumeFraction,
    seed = poreSeed,
    periodic = pfw["periodic"],
    minLigament = poreMinLigament,
    maxTrialsPerPore = 20000,
    dim = 3,
    # geom.box-style order: [x-, y-, z-, x+, y+, z+].
    flaggedSurfaces = [False, False, True, False, False, True],
    vel = [0.0, 0.0, -impactVelocity],
    mat = COPPER,
    group = 0,
    particleType = 1,  # SinglePointBSpline
)

pfw["objects"] = [foam]

# =============================================================================
# Boundary conditions
# =============================================================================

pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0
pfw["fTableInterpType"] = "Linear"

# Boundary-condition order is x-, x+, y-, y+, z-, z+.
# x/y are periodic; z- is the impact symmetry plane; z+ is an outflow/free face.
pfw["boundaryConditionTypes"] = [0, 0, 0, 0, 1, 0]
