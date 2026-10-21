# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Two-disk collision example.
# Solver/PFW features: Multiple geometry objects, moving material points, and explicit dynamic update.
# Workflow note: keep this file as a copyable problem definition.  Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
# [pfw_dependency] pfw:pfw_materials.py
import importlib
matdb = importlib.import_module('pfw_materials')

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
stopTime = 2.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.93 # mass density mg/mm^3
E = 10.0     # in-plane stiffness (GPa)
nu = 0.27    # in plane poisson's ratio

# Domain ---------------------------------------------------------------------------------
domainWidth = 2.0
domainHeight = 2.0

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainHeight # mm

pfw["formulation"] = "PlaneStrain"

cpp = 20 # cells per partition in each direction

refine = 2
pfw["xpar"] = refine  # grid partitions
pfw["ypar"] = refine
pfw["zpar"] = 1

pfw["nI"] = pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"] = pfw["ypar"]*cpp  # grid cells in the y-direction

domainLength = domainHeight/refine/cpp

# Physical out-of-plane measure; independent of object end caps and z cells.
pfw["planeStrainThickness"] = domainLength

pfw["outputType"] = "silo"
pfw["gridFieldNames"] = [
    "gridMass",
    "gridVelocity",
    "gridInternalForce",
    "gridDamageGradient",
]

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime/200
pfw["restartInterval"] = stopTime*100 # Don't need restarts for now

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1e-16
pfw["cpdiDomainScaling"] = 1
pfw["damageFieldPartitioning"] = 1

pfw["solverProfiling"] = 0
pfw["needsNeighborList"] = 0
pfw["reactionHistory"] = 1
pfw["boxAverageHistory"] = 1
pfw["frictionCoefficient"] = 0.25

pfw["maxParticleVelocity"] = 10.0
pfw["minParticleJacobian"] = 0.01
pfw["maxParticleJacobian"] = 10.0

pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 4

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [ 0, 0, 0, 0, 0, 0 ]
pfw["prescribedBoundaryFTable"] = 0

# Define all the geometric objects -------------------------------------------------------
r1 = 0.2*domainHeight
r2 = r1
disk1 = geom.disk('cyl1', [pfw['xmin'] + r1, pfw['ymin'] + r1], r1, vel=[1.0, 1.0, 0.0], mat=0, group=0)
disk2 = geom.disk('cyl2', [pfw['xmax'] - r2, pfw['ymax'] - r2], r2, vel=[-1.0, -1.0, 0.0], mat=0, group=0)
pfw["objects"] = [disk1,disk2]

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mWallTime"] = "00:10:00"
pfw["mSubmitJobs"] = True

# GEOS MPM i/o parameters ---------------------------------------------------------------

# UNPORTED MPM SOLVER PARAMETERS

pfw["materials"] = [ matdb.elasticDemo["name"] ]
pfw["materialPropertyString"] = matdb.elasticDemo["materialString"]
