# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 2.0

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 0.01  # mm
sampleHeight = 0.01 # mm
sampleLength = 0.01 # mm

domainWidth = 1.1*sampleWidth
domainHeight = sampleHeight
domainLength = 1.1*sampleLength

pfw["xmin"] = -0.5*domainWidth  # mm
pfw["xmax"] =  0.5*domainWidth  # mm
pfw["ymin"] =  0.0              # mm
pfw["ymax"] =  domainHeight     # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

refine=1 # partitions in each direction
cpp=3    # cells per partition in each direction
pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["useEvents"]=1
pfw["frictionCoefficient"]=0.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "mat1" ]
pfw["materialPropertyString"] ="""
<ElasticIsotropic
	name="mat1"
	defaultDensity="2.7"
	defaultBulkModulus="67"
	defaultShearModulus="26"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-sampleWidth/2,0.0,-sampleLength/2],[sampleWidth/2,sampleHeight,sampleLength/2],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
			   [stopTime/2,   1.00,  0.95,  1.00],
			   [stopTime, 	  1.00,  0.95,  1.00]]

# MPM EVENTS -------------------------------------------------------------------------------

pfw["mpmEventsString"]="""
<MPMEvents>
	<Anneal 
	time=""" + '"' + str(stopTime/2) + '"' + """
	interval=""" + '"' + str(stopTime/2) + '"' + """
	source="all" />
</MPMEvents>
"""
