# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import pfw_materials as matdb

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs using the buckling capability with non-monotonic cascading
# compaction.
#

pfw = {}
pfw["runDebug"] = False
stopTime = 5.0

maxCompressiveStrain = 0.40

# material direction:
matDir = np.array([1.,0.,1.])
matDir = matDir / np.linalg.norm(matDir)

# DOMAIN ---------------------------------------------------------------------------------

sampleX = 1.0  # mm
sampleY = 1.0 # mm
sampleZ = 1.0 # mm

domainX = 1.2*sampleX  # This would be increased for unconfined compression.
domainY = 1.2*sampleY
domainZ = sampleZ

pfw["xmin"] = 0.0             # mm
pfw["xmax"] = domainX    # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainY # mm
pfw["zmin"] = 0.0 # mm
pfw["zmax"] = domainZ # mm

refine = 1  # partitions in each direction
cpp = 8    # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# These fields are needed
pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "StrengthScale",       # needed for weibull
                             "SurfaceFlag",         # needed for CZ and contact
                             "MaterialDirection",   # needed for graphite model
                             "RVector"]             # needed for cpdi and plotting bspline

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.02
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"] = stopTime/1000
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime/1000

pfw["solverProfiling"]=0
pfw["frictionCoefficient"]=0.25

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

# MATERIAL PROPERTIES --------------------------------------------------------------------
pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-0.5*sampleX,-0.5*sampleY,pfw["zmin"]],[0.5*sampleX,0.5*sampleY,pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)

blockWDir = geom.materialDirectionWrapper(name='blockWDir',subObject=block,matDir=matDir)

pfw["objects"]=[blockWDir]

# DEFORMATION -----------------------------------------------------------------------------

# Ftable only controls z-direction
pfw["boundaryConditionTypes"]=[ 0, 0, 0, 0, 2, 2 ]
pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"]=[
    [0,	       1.0,	1.0, 1.0],
    [stopTime, 1.0, 1.0, np.exp(-maxCompressiveStrain) ]
    ]
