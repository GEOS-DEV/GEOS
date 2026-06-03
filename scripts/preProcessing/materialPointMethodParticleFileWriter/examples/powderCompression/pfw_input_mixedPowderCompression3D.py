# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: 3D reverse-ballistic uniaxial compression of a packed spherical granular bed.
# Solver/PFW features: Periodic lateral y/z boundaries, packedSphericalBed particles,
# SinglePointBSpline particles, surface contact fields, and multi-material grains.
# Workflow note: keep this file as a copyable problem definition. Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

"""3D reverse-ballistic compression of a packed spherical bed.

This 3D version uses a full x-y-z domain.  The impact axis is x, matching the
velocity used below:

    vel = [-impactVelocity, 0, 0]

The material starts with velocity toward the negative-x symmetry plane at x=0.
The lateral y and z faces are periodic.  The packedSphericalBed object generates
true 3D spheres inside the box and is periodic in y and z.
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

# Initial aspect ratio is 6:1, longer in the impact direction x.
AR = 6.0

sampleX = AR * 1.0      # mm, impact direction
sampleY = 1.0           # mm, lateral periodic direction
sampleZ = 1.0           # mm, lateral periodic direction

# Packed-bed grain controls.  In 3D, volumeFraction is a true volume fraction.
grainSeed = 17

# Reverse-ballistic impact speed.  Positive value here means motion toward the
# negative x symmetry plane.
impactVelocity = 0.25   # mm/us = 250 m/s

# Resolution controls.  cpp, xpar, ypar, and zpar are scaled by refine.
refine = 1

stopTime = sampleX / impactVelocity

# =============================================================================
# Domain and grid
# =============================================================================

pfw["xmin"] = 0.0
pfw["xmax"] = sampleX

pfw["ymin"] = -0.5 * sampleY
pfw["ymax"] =  0.5 * sampleY

pfw["zmin"] = -0.5 * sampleZ
pfw["zmax"] =  0.5 * sampleZ

# Impact direction x is non-periodic.  Lateral y and z are periodic.
pfw["periodic"] = [False, True, True]

# Full 3D calculation, not plane strain.
pfw["planeStrain"] = 0

pfw["xpar"] = int(AR * 3 * refine)
pfw["ypar"] = int(3 * refine)
pfw["zpar"] = int(3 * refine)

cpp = 24
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = pfw["zpar"] * cpp

# In 3D, ppc is particles per cell per active direction.
# ppc=2 gives 8 material points per full 3D cell before geometry filtering.
pfw["ppc"] = 2

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

# Use damage-field partitioning / DFG for contact and possible damage separation.
pfw["damageFieldPartitioning"] = 1

# Explicit contact settings.  SurfacePosition and SurfaceNormal are generated
# below from the spherical grain surfaces.
pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"] = 1
pfw["explicitSurfaceNormalInfluence"] = 1000.0

# =============================================================================
# Profile history output
# =============================================================================

pfw["profileHistory"] = 1

# Slice normal to the impact direction.
pfw["profileDirection"] = "x"

# 0 means use the background-grid resolution in x.
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

pfw["materials"] = [matdb.copper["name"], matdb.quartz["name"]]
pfw["materialPropertyString"] = matdb.copper["materialString"] + matdb.quartz["materialString"]

COPPER = 0
QUARTZ = 1

# =============================================================================
# Geometry
# =============================================================================

bed = geom.packedSphericalBed(
    "bed",

    x0=[pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    x1=[pfw["xmax"], pfw["ymax"], pfw["zmax"]],

    dim=3,

    # Periodic in the lateral directions only.
    periodic=[False, True, True],

    materials=[
        {
            "mat": COPPER,
            "volumeFraction": 0.45,
            "meanDiameter": 0.25,
            "stdDiameter": 0.025,
        },
        {
            "mat": QUARTZ,
            "volumeFraction": 0.15,
            "meanDiameter": 0.15,
            "stdDiameter": 0.015,
        },
    ],

    # Initial material velocity: reverse-ballistic motion toward x = 0.
    vel=[-impactVelocity, 0.0, 0.0],

    overlapPercent=2.0,
    method="auto",
    seed=grainSeed,

    # Uniform contact group and particle type for the whole bed.
    group=0,
    particleType=1,
)

pfw["objects"] = [bed]

# =============================================================================
# Boundary conditions
# =============================================================================

# Boundary-condition order is x-, x+, y-, y+, z-, z+.
# Integer options in the MPM solver are:
#   Outflow  = 0
#   Symmetry = 1
#   Moving   = 2
#   Contact  = 3
#
# x- is the impact symmetry plane.
# x+ is outflow/free.
# y and z are periodic, so their boundary-condition entries are effectively
# bypassed by the periodic setting.
pfw["boundaryConditionTypes"] = [1, 0, 0, 0, 0, 0]