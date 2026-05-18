# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import pfw_materials as matdb

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs. 
#
# TODO: add some actual verification, so make sure the response is correct.


pfw = {}
pfw["runDebug"] = False
stopTime = 100.0

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 1.0  # mm
sampleHeight = 1.0 # mm
sampleLength = 1.0 # mm

domainWidth = sampleWidth  # This would be increased for unconfined compression.
domainHeight = sampleHeight
domainLength = sampleLength

pfw["xmin"] = 0.0             # mm
pfw["xmax"] = domainWidth    # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainHeight # mm
pfw["zmin"] = 0.0 # mm
pfw["zmax"] = domainLength # mm

refine=1  # partitions in each direction
cpp=3     # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:05:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25  
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

# MATERIAL PROPERTIES --------------------------------------------------------------------
pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]

# GEOMETRY OBJECTS -------------------------------------------------------
block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------

pfw["fTableInterpType"]='Smoothstep'
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,             1.000, 1.000, 1.000],
               [.1*stopTime, 0.995, 0.995, 0.995],
               [.2*stopTime, 0.998, 0.998, 0.998],
               [.3*stopTime, 0.980, 0.980, 0.980],
               [.4*stopTime, 0.985, 0.985, 0.985],
               [.5*stopTime, 0.960, 0.960, 0.960],
               [.6*stopTime, 0.965, 0.965, 0.965],
               [.7*stopTime, 0.940, 0.940, 0.940],
               [.8*stopTime, 0.945, 0.945, 0.945],
               [.9*stopTime, 0.920, 0.920, 0.920],
               [stopTime,    0.925, 0.925, 0.925]]

# prescribed deformation (moving pistons) at all faces
pfw["boundaryConditionTypes"]=[2, 2, 2, 2, 2, 2]