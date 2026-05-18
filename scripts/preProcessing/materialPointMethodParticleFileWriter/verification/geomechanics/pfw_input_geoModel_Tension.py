# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import pfw_materials as matdb

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs. 
#
# Initializes material with a confining stress (should be p < -p0/3 ) and then ramps
# load in 1 direction with that confinement maintained through stress control on
# lateral boundaries.
# TODO: add some actual verification, so make sure the response is correct.


pfw = {}
pfw["runDebug"] = False
stopTime = 1000.0

physicalStopTime = 360*24*3600.e6 # actual end of 1-year test time, micro-seconds.
timeScale = physicalStopTime / stopTime  # strain rate multiplier to match 1-year creep test in stopTime

#confiningPressure = 0.005  # confining pressure GPa
#maxCompressiveStrain = 0.05 # stress control end point (GPa)
confiningPressure = 0.00  
maxCompressiveStrain = -0.1

fractureEnergyReleaseRate = 1.e-10
fractureStress = 0.005

kp = 0.05
ki = 0.00
kd = 0.005
pfw["cflFactor"]=0.005



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
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"

pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"] = stopTime/1000
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime/1000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="PIC"
pfw["updateOrder"]=2

pfw["useEvents"] = 1

# MATERIAL PROPERTIES --------------------------------------------------------------------
pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------
# Ftable only controls z-direction
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]
pfw["fTableInterpType"]='Cosine'
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"]=[
    [0,	          1.000,	    1.,	1.],
    [stopTime,	    1.0,1.0, np.exp(-maxCompressiveStrain) ]
    ]

