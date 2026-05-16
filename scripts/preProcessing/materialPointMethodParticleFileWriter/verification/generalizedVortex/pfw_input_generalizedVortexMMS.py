# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 4

# DOMAIN ---------------------------------------------------------------------------------

pfw["planeStrain"] = 1

refine=5 # partitions in each direction
cpp=20   # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp
pfw["nJ"]=pfw["ypar"]*cpp
pfw["nK"]=3
pfw["ppc"]=2

domainHeight = 4
domainWidth = 4 

pfw["xmin"] = -0.5*domainHeight # mm
pfw["xmax"] = 0.5*domainHeight # mm
pfw["ymin"] =-0.5*domainWidth # mm
pfw["ymax"] = 0.5*domainWidth # mm

dx = domainWidth/(pfw["nI"]-2)
dy = domainHeight/(pfw["nJ"]-2)

domainLength = 3*dx

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# PATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 100
pfw["restartInterval"] = stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.05   
pfw["initialDt"]=1e-16

pfw["cpdiDomainScaling"] = 0
pfw["bodyForce"] = [ 0, 0, 0 ]
pfw["generalizedVortexMMS"] = 1
pfw["reactionHistory"] = 0
pfw["boxAverageHistory"] = 0

pfw["boundaryConditionTypes"] = [ 0, 0, 0, 0, 0, 0 ]
pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0  

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1000
firstLameConstant = 577
shearModulus = 384.615384615384585

pfw["materials"] = [ "hyperelasticMMS" ]
pfw["materialPropertyString"]="""
<HyperelasticMMS
	name="hyperelasticMMS"
	defaultDensity=""" + '"' + str(density) + '"' + """
	defaultLambda=""" + '"' + str(firstLameConstant) + '"' + """
	defaultShearModulus=""" + '"' + str(shearModulus) + '"' + """/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('test_specimen',[pfw["xmin"]+2*dx,pfw["ymin"]+2*dx,pfw["zmin"]],[pfw["xmax"]-2*dx,pfw["ymax"]-2*dx,pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]
