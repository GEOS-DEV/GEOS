# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 1.0

# DOMAIN ---------------------------------------------------------------------------------

refine=5
cpp=5

pfw["periodic"] = [True, True, False]

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=3  			# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainWidth = 1.0
domainHeight = 1.0 # m
domainLength = domainWidth*(pfw["nK"]-2)/(pfw["nI"]-2)  # m, to get cubic cells

pfw["xmin"] =-0.5*domainWidth	# m
pfw["xmax"] = 0.5*domainWidth	# m
pfw["ymin"] =-0.5*domainHeight	# m
pfw["ymax"] = 0.5*domainHeight	# m
pfw["zmin"] =-0.5*domainLength 	# m
pfw["zmax"] = 0.5*domainLength 	# m

lx = pfw["xmax"] - pfw["xmin"]
ly = pfw["ymax"] - pfw["ymin"]
lz = pfw["zmax"] - pfw["zmin"]

pfw["planeStrain"] = 1

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:06:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False
# pfw["autoRestart"]=True

# GEOSX MPM SOLVER PARAMETERS ------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.5    
pfw["initialDt"]=1e-16

pfw["solverProfiling"]=0
pfw["neighborRadius"]=-1.01
pfw["needsNeighborList"]=0

pfw["boundaryConditionTypes"]=[0, 0, 0, 0, 1, 1 ]  

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2.7"
	defaultBulkModulus="70.0"
	defaultShearModulus="24.0"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

disk1 = geom.cylinder('disk1',[0.0, 0.0, pfw["zmin"]],[0.0, 0.0, pfw["zmax"]], 0.4, vel=[-2.0,-2.0,0.0], mat=0, group=0)
pfw["objects"]=[disk1]
