# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as matdb

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs.
#
# TODO: add some actual verification, so make sure the response is correct.


pfw = {}
pfw["runDebug"] = False
stopTime = 1000.0

physicalStopTime = 360*24*3600.e6 # actual end of 1-year test time, micro-seconds.
timeScale = physicalStopTime / stopTime  # strain rate multiplier to match 1-year creep test in stopTime

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

refine = 1  # partitions in each direction
cpp = 3     # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=1               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100

# Silo output is required by the verification-suite VisIt smoke renderer.
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity"]
pfw["restartInterval"]=stopTime/25
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/1000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/1000

pfw["solverProfiling"]=0
pfw["frictionCoefficient"]=0.25

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2
# pfw["constantTimeStepValue"] = stopTime/3000  # obsolete in current GEOS-MPM

# MATERIAL PROPERTIES --------------------------------------------------------------------
pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------

pfw["fTableInterpType"] = "Smoothstep"
pfw["prescribedBoundaryFTable"]=1

# copy strain data from hydrostatic test for ghareb
volStrainBauer=np.array([0,-0.049,-0.044,-0.138,-0.132,-0.194,-0.188,-0.221,-0.21,-0.24,-0.23,-0.26,-0.27,-0.265,-0.297,-0.309,-0.317,-0.22])
pfw["fTable"]=[[ ii/(len(volStrainBauer)-1)*stopTime, np.exp(v)**(1/3),np.exp(v)**(1/3),np.exp(v)**(1/3)] for ii,v in enumerate(volStrainBauer)]


# prescribed deformation (moving pistons) at all faces
pfw["boundaryConditionTypes"]=[2, 2, 2, 2, 2, 2]
