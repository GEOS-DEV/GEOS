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
# [pfw_dependency] pfw_materials.py
import importlib
matdb = importlib.import_module('pfw_materials')

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = True
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

pfw["planeStrain"] = 1

cpp = 20 # cells per partition in each direction

refine = 2
pfw["xpar"] = refine  # grid partitions
pfw["ypar"] = refine
pfw["zpar"] = 1

pfw["nI"] = pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"] = pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"] = 3          # grid cells in the z-direction

domainLength = domainHeight/refine/cpp

pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

pfw["outputType"] = "silo"

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

pfw["updateMethod"] = "XPIC"
pfw["updateOrder"] = 2

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [ 0, 0, 0, 0, 0, 0 ]
pfw["prescribedBoundaryFTable"] = 0

# Define all the geometric objects -------------------------------------------------------
r1 = 0.2*domainHeight
r2 = r1
disk1 = geom.cylinder('cyl1',[pfw["xmin"] + r1,pfw["ymin"]+r1,pfw["zmin"]-2*domainLength],[pfw["xmin"] + r1,pfw["ymin"]+r1,pfw["zmax"]+2*domainLength],r1,vel=[1.0,1.0,0.],mat=0,group=0)
disk2 = geom.cylinder('cyl2',[pfw["xmax"] - r2,pfw["ymax"]-r2,pfw["zmin"]-2*domainLength],[pfw["xmax"] - r2,pfw["ymax"]-r2,pfw["zmax"]+2*domainLength],r2,vel=[-1.0,-1.0,0.],mat=0,group=0)
pfw["objects"] = [disk1,disk2]

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mWallTime"] = "12:00:00"
pfw["mCores"] = pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"] = int(np.ceil(float(pfw["mCores"])/36.))
pfw["mSubmitJobs"] = False

# GEOS MPM i/o parameters ---------------------------------------------------------------

# UNPORTED MPM SOLVER PARAMETERS

pfw["materials"] = [ matdb.elasticDemo["name"] ]
pfw["materialPropertyString"] = matdb.elasticDemo["materialString"]
