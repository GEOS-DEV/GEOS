# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Compressed elastic disk example.
# Solver/PFW features: ElasticIsotropic response, F-table boundary motion, and reaction-history output.
# Workflow note: keep this file as a copyable problem definition.  Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
# [pfw_dependency] pfw_materials.py
import importlib
matdb = importlib.import_module('pfw_materials')

pfw = {}
pfw["runDebug"] = True
stopTime = 20

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.93 # mass density mg/mm^3
E = 45.0     # in-plane stiffness (GPa)
nu = 0.27    # in plane poisson's ratio

# Domain ---------------------------------------------------------------------------------

pfw["formulation"] = "PlaneStrain"

refine = 4
cpp = 10 # cells per partition in each direction
pfw["xpar"] = refine  # grid partitions
pfw["ypar"] = refine
pfw["zpar"] = 1

pfw["nI"] = pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"] = pfw["ypar"]*cpp  # grid cells in the y-direction

domainWidth = 1.0
domainHeight = domainWidth
domainLength = domainHeight*(3-2)/(pfw["nJ"]-2)  # m, to get cubic cells

# Define all the geometric objects -------------------------------------------------------

pfw["xmin"] = -0.5*domainWidth   # m
pfw["xmax"] = 0.5*domainWidth   # m
pfw["ymin"] = 0.0   # m
pfw["ymax"] = domainHeight  # m
# Physical out-of-plane measure; independent of object end caps and z cells.
pfw["planeStrainThickness"] = domainLength

pfw["outputType"] = "silo"

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime/200
pfw["restartInterval"] = stopTime*10 # Don't need restarts for now

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
pfw["updateOrder"] = 2

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [ 2, 2, 2, 2, 1, 1 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [[0,            1.0,  1.0,    1.0],
               [stopTime,     1.0,  0.5,    1.0]]

# Define all the geometric objects -------------------------------------------------------

disk = geom.disk('disk', [0.0, domainHeight / 2], 0.5 * domainHeight, vel=[0.0, 0.0, 0.0], mat=0, group=0)
pfw["objects"] = [disk]

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mWallTime"] = "12:00:00"
pfw["mCores"] = pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"] = False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = [ matdb.elasticDemo["name"] ]
pfw["materialPropertyString"] = matdb.elasticDemo["materialString"]
