# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 2.0

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth=0.01
sampleHeight=sampleWidth
sampleLength=sampleWidth

domainWidth = sampleWidth
domainHeight = sampleHeight
domainLength = sampleLength

pfw["xmin"] =-0.5*domainWidth	# m
pfw["xmax"] = 0.5*domainWidth	# m
pfw["ymin"] =-0.5*domainHeight	# m
pfw["ymax"] = 0.5*domainHeight	# m
pfw["zmin"] =-0.5*domainLength  # m
pfw["zmax"] = 0.5*domainLength  # m

periodic = [True, True, True]

refine=5
cpp=5

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine


pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp	# grid cells in the z-direction
pfw["ppc"]=2   		        # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=0.001
pfw["restartInterval"]=stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.5    
pfw["initialDt"]=1e-16
pfw["solverProfiling"]=0
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0

pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
			   
# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2.648"
	defaultBulkModulus="36.3"
	defaultShearModulus="26.0" />
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('box',[pfw["xmin"], pfw["ymin"], pfw["zmin"]],[pfw["xmax"], pfw["ymax"], pfw["zmax"]], vel=[0.0,0.0,0.0], mat=0, group=0)
pfw["objects"]=[block]

# DEFORMATION ----------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 0, 0, 0, 0 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedFTable"]=1
pfw["prescribedBoundaryFTable"]=0
pfw["fTable"]=[[0,        1.00, 1.00, 1.00],
		       [stopTime, 0.95, 1.00, 1.00]]
