# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 0.1

# Domain ---------------------------------------------------------------------------------
refine=6
cpp=15
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=3		 	# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainWidth= 10.0 # m
domainHeight= domainWidth
domainLength = domainWidth*(pfw["nK"]-2)/(pfw["nI"]-2)  # m, to get cubic cells

# Define all the geometric objects -------------------------------------------------------

pfw["xmin"] =-0.5*domainWidth  # m
pfw["xmax"] = 0.5*domainWidth  # m
pfw["ymin"] =-0.5*domainHeight # m
pfw["ymax"] = 0.5*domainHeight # m
pfw["zmin"] =-0.5*domainLength # m
pfw["zmax"] = 0.5*domainLength # m

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=True

# GEOS MPM i/o parameters ---------------------------------------------------------------

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25   
pfw["initialDt"]=1e-16
pfw["planeStrain"]=1

pfw["solverProfiling"]=0
pfw["cpdiDomainScaling"]=1       
pfw["contactGapCorrection"]=1
pfw["frictionCoefficient"]=0.25

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 1, 1, 1, 1, 1, 1 ]    

pfw["prescribedBoundaryFTable"]=0
pfw["fTableInterpType"]=0  

# Define all the geometric objects -------------------------------------------------------

crop=1.0
plate1 = geom.box('plate1',[pfw["xmin"],crop*pfw["ymin"],pfw["zmin"]],[0,crop*pfw["ymax"],pfw["zmax"]],vel=[100,0,0],mat=0,group=0)
plate2 = geom.box('plate2',[0,crop*pfw["ymin"],pfw["zmin"]],[pfw["xmax"],crop*pfw["ymax"],pfw["zmax"]],vel=[-100,0,0],mat=0,group=0)
hole = geom.cylinder('hole',[0,0,pfw["zmin"]],[0,0,pfw["zmax"]],2,[0,0,0],0,0,0)
plateWithHole1 = geom.difference(plate1,hole)
plateWithHole2 = geom.difference(plate2,hole)

pfw["objects"]=[plateWithHole1,plateWithHole2]

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e8"
	defaultShearModulus="24.0e8"/>
"""
