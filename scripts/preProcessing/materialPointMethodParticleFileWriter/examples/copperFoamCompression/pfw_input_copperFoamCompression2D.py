# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Plane-strain reverse-ballistic uniaxial compression of copper foam.
# Solver/PFW features: Periodic lateral x boundary, Poisson-disk circular pores, SinglePointBSpline particles, and copper VonMisesJ plasticity.
# Workflow note: keep this file as a copyable problem definition. Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

"""2D plane-strain copper foam reverse-ballistic compression.

GEOS-MPM plane strain in this PFW branch uses the x-y plane and a single
through-thickness z particle layer.  Therefore the 2D version maps the impact
axis to y: material points start with velocity ``[0, -impactVelocity, 0]`` and
impact the symmetry plane at ``y=0``.  The x faces are periodic, and the
Poisson-disk pore field is periodic in x.  The circular pores are implicitly
extruded through the plane-strain thickness.
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

# Initial aspect ratio is 6:1, longer in the in-plane impact direction y.
lateralLength = 1.0          # mm
impactAspectRatio = 6.0
domainSize = {
    "x": lateralLength,
    "y": impactAspectRatio * lateralLength,
}

# Circular plane-strain pore controls.  The realized area fraction is printed by
# the geometry object when PFW imports this input.
poreDiameter = 0.15          # mm
poreAreaFraction = 0.18      # target void area fraction in the x-y plane
poreMinLigament = 0.01       # mm additional clearance between pore surfaces
poreSeed = 20260521

# Reverse-ballistic impact speed.  Positive value here means motion toward the
# negative y symmetry plane.
impactVelocity = 0.25        # mm/us = 250 m/s

# Resolution controls.  cpp, xpar, and ypar are scaled by refine.  PFW requires
# zpar=1 and nK=3 for planeStrain=1, so through-thickness partitioning is fixed.
refine = 1
cpp = 24
xpar = 3
ypar = 18
zpar = 1
ppc = 2                      # particles per cell per active direction

# End time estimate: snowplow compaction wave speed in a porous bar.  The 2D
# area fraction is used as the relative missing density in this plane-strain demo.
lockedPorosity = 0.0
initialRelativeDensity = 1.0 - poreAreaFraction
lockedRelativeDensity = 1.0 - lockedPorosity
compactionStrain = max(1.0e-12, 1.0 - initialRelativeDensity / lockedRelativeDensity)
compactionWaveSpeed = impactVelocity / compactionStrain
compactionWaveTransitTime = domainSize["y"] / compactionWaveSpeed
stopTime = compactionWaveTransitTime

# =============================================================================
# Domain and grid
# =============================================================================

pfw["xmin"] = -0.5 * domainSize["x"]
pfw["xmax"] =  0.5 * domainSize["x"]
pfw["ymin"] = 0.0
pfw["ymax"] = domainSize["y"]

pfw["periodic"] = [True, False, False]
pfw["planeStrain"] = 1

cppRefined = max(1, int(round(refine * cpp)))
pfw["xpar"] = max(3 if pfw["periodic"][0] else 1, int(round(refine * xpar)))
pfw["ypar"] = max(3 if pfw["periodic"][1] else 1, int(round(refine * ypar)))
pfw["zpar"] = 1

pfw["nI"] = pfw["xpar"] * cppRefined
pfw["nJ"] = pfw["ypar"] * cppRefined + 2
pfw["nK"] = 3
pfw["ppc"] = ppc

# Choose the z thickness so the single plane-strain interior cell is roughly the
# same size as the in-plane cells.
domainThickness = domainSize["y"] / (pfw["nJ"] - 2)
pfw["zmin"] = -0.5 * domainThickness
pfw["zmax"] =  0.5 * domainThickness

pfw["sortObjects"] = True
pfw["outputType"] = "silo"

# =============================================================================
# Batch settings used by particleFileWriter.py
# =============================================================================

pfw["mWallTime"] = "00:10:00"
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
    "periodic_copper_foam_2d",
    x0 = [pfw["xmin"], pfw["ymin"]],
    x1 = [pfw["xmax"], pfw["ymax"]],
    poreDiameter = poreDiameter,
    porosity = poreAreaFraction,
    seed = poreSeed,
    periodic = [pfw["periodic"][0], pfw["periodic"][1]],
    minLigament = poreMinLigament,
    maxTrialsPerPore = 20000,
    dim = 2,
    # geom.box-style order for dim=2: [x-, y-, x+, y+].
    flaggedSurfaces = [False, True, False, True],
    vel = [0.0, -impactVelocity, 0.0],
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
# x is periodic; y- is the impact symmetry plane; y+ is outflow/free; z faces are
# symmetry planes for plane strain.
pfw["boundaryConditionTypes"] = [0, 0, 1, 0, 1, 1]
