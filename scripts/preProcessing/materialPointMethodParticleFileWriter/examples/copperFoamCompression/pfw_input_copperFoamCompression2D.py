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
AR = 6.

sampleY = 1.0          # mm
sampleX = AR*1.0          # mm

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

stopTime = sampleX / impactVelocity

# =============================================================================
# Domain and grid
# =============================================================================

pfw["xmin"] = 0.0
pfw["xmax"] = sampleX
pfw["ymin"] = -0.5*sampleY
pfw["ymax"] = 0.5*sampleY

pfw["periodic"] = [False, True, False]
pfw["planeStrain"] = 1

pfw["xpar"] = int(AR*3*refine)
pfw["ypar"] = int(3*refine)
pfw["zpar"] = 1

cpp = 24
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = 2

# Choose the z thickness so the single plane-strain interior cell is roughly the
# same size as the in-plane cells.
domainZ = (pfw["xmax"] - pfw["xmin"]) / ( pfw["nI"] - 2 )
pfw["zmin"] = -0.5 * domainZ
pfw["zmax"] =  0.5 * domainZ

pfw["outputType"] = "silo"

# =============================================================================
# Batch settings used by particleFileWriter.py
# =============================================================================

pfw["mWallTime"] = "00:10:00"
pfw["mSubmitJobs"] = True

# =============================================================================
# GEOS MPM solver controls
# =============================================================================

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 64
pfw["restartInterval"] = 2.0 * stopTime

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-16
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 2

# Use damage-field partitioning/DFG for usto contact.  Will also allow damage separation
# If the Johnson-Cook damage threshold is reached.
pfw["damageFieldPartitioning"] = 1

# Explicit contact settings.  SurfacePosition and SurfaceNormal are generated
# below from the pore surfaces.
pfw["contactGapCorrection"] = "Implicit" # Activate contact only if gap <=0
pfw["useSurfacePositionForContact"] = 1
pfw["explicitSurfaceNormalInfluence"] = 19.68758675709231

# =============================================================================
# Profile history output
# =============================================================================

pfw["profileHistory"] = 1

# Literal x-direction profile. This slices the domain normal to x and writes
# one CSV file per requested variable.
pfw["profileDirection"] = "x"

# 0 means use the background-grid resolution in x.
# For this copperFoamCompression setup, that corresponds to the x grid spacing.
pfw["profileNumSlices"] = 0

# Write at the same cadence as reactionHistory and boxAverageHistory.
# Do not also set profileCycleInterval > 0.
pfw["profileWriteInterval"] = stopTime / 1000.0
pfw["profileCycleInterval"] = 0

# String form is safest with the current PFW pass-through behavior.
pfw["profileVariables"] = "{ density, volumeFraction, velocityX, velocityY, velocityZ, kineticEnergy, damage, stressXX, stressYY, stressZZ, plasticStrainMagnitude }"

# Limits for particle deletion.  
pfw["maxParticleVelocity"] = 10.0
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
    vel = [ -impactVelocity, 0.0, 0.0],
    mat = COPPER,
    group = 0,
    particleType = 1,  # SinglePointBSpline
)

pfw["objects"] = [foam]

# =============================================================================
# Boundary conditions
# =============================================================================

# Boundary-condition order is x-, x+, y-, y+, z-, z+.
# x is periodic; y- is the impact symmetry plane; y+ is outflow/free; z faces are
# symmetry planes for plane strain.
pfw["boundaryConditionTypes"] = [1, 0, 0, 0, 1, 1]
