# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 0.01

# Domain ---------------------------------------------------------------------------------
refine=2
cpp=5

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=3  			# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainLength = 1.0 # m
domainThickness = domainLength*(pfw["nK"]-2)/(pfw["nI"]-2)  # m, to get cubic cells

# Define all the geometric objects -------------------------------------------------------
pfw["xmin"] =-0.5*domainLength	# m
pfw["xmax"] = 0.5*domainLength	# m
pfw["ymin"] =-0.5*domainLength	# m
pfw["ymax"] = 0.5*domainLength	# m
pfw["zmin"] =-0.5*domainThickness 	# m
pfw["zmax"] = 0.5*domainThickness 	# m

lx = pfw["xmax"] - pfw["xmin"]
ly = pfw["ymax"] - pfw["ymin"]
lz = pfw["zmax"] - pfw["zmin"]

pfw["planeStrain"] = 1

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

# GEOSX MPM PARAMETERS ------------------------------------------------------------------

pfw["endTime"]=stopTime				# seconds
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime/4

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.5    
pfw["initialDt"]=1e-16

pfw["solverProfiling"]=0
pfw["neighborRadius"]=-1.01
pfw["needsNeighborList"]=1
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1
pfw["contactGapCorrection"] = "Implicit"
pfw["frictionCoefficient"]=0.25

pfw["prescribedBoundaryFTable"]=0

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[0, 0, 0, 0, 0, 0 ]  

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

disk1 = geom.cylinder('disk1',[0.0, 0.0, pfw["zmin"]],[0.0, 0.0, pfw["zmax"]], 0.3, vel=[-120.0,-120.0,0.0], mat=0, group=0)
pfw["objects"]=[disk1]

pfw["materials"] = [ "aluminum", "aluminum10x" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e8"
	defaultShearModulus="24.0e8"/>
<ElasticIsotropic
	name="aluminum10x"
	defaultDensity="27000"
	defaultBulkModulus="70.0e9"
	defaultShearModulus="24.0e9"/>
"""
