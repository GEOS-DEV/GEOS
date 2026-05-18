# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 0.01

# DOMAIN ---------------------------------------------------------------------------------

refine=3 # partitions in each direction
cpp=16 # grid cells per direction for each partition

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=pfw["ypar"]*cpp  	# grid cells in the z-direction
pfw["ppc"]=2   		        # particles per cell in each direction

domainLength = 1.0 # m
domainThickness = domainLength*(pfw["nK"]-2)/(pfw["nI"]-2)  # m, to get cubic cells

pfw["xmin"] =-0.5*domainLength	   # m
pfw["xmax"] = 0.5*domainLength	   # m
pfw["ymin"] =-0.5*domainLength	   # m
pfw["ymax"] = 0.5*domainLength	   # m
pfw["zmin"] =-0.5*domainThickness  # m
pfw["zmax"] = 0.5*domainThickness  # m

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS ------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/1000.0001
pfw["restartInterval"]=stopTime/4

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25    
pfw["initialDt"]=1e-16

pfw["contactGapCorrection"] = "Implicit"
pfw["frictionCoefficient"]=0.0 

# MATERIAL PROPERTIES ---------------------------------------------------------

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

# GEOMETRY OBJECTS -------------------------------------------------------

ball1 = geom.sphere('ball1',[pfw["xmax"]/2,pfw["ymax"]/2,pfw["zmax"]/2],0.15,vel=[-60.0,-60.0,-60.0],mat=0,group=0)
ball2 = geom.sphere('ball2',[-pfw["xmax"]/2,-pfw["ymax"]/2,-pfw["zmax"]/2],0.15,vel=[6.0,6.0,6.0],mat=1,group=1)
pfw["objects"]=[ball1,ball2]

