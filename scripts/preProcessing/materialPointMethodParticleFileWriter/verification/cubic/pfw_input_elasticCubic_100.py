# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 10.0

# DOMAIN ---------------------------------------------------------------------------------

sampleHeight = 1.0 # mm
sampleWidth = 1.0  # mm
sampleLength = 1.0 # mm 

domainHeight = sampleHeight
domainWidth = sampleWidth
domainLength = sampleLength

refine=1 # partitions in each direction
cpp=5   # cells per partition in each direction

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=1   # particles per cell in each direction

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] =  0.0 # mm
pfw["ymax"] =  domainHeight # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"] = stopTime
pfw["plotInterval"] = 1.0
pfw["restartInterval"] = stopTime*5.0

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.05    
pfw["initialDt"]=1e-16

pfw["solverProfiling"]=0
pfw["planeStrain"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["cpdiDomainScaling"]=0
pfw["damageFieldPartitioning"]=0

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
							 "MaterialDirection"]

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-sampleWidth/2, 0.0, -sampleLength/2],[sampleWidth/2, sampleHeight, sampleLength/2],vel=[0.0,0.0,0.0],mat=0,group=0)
blockWDir = geom.materialDirectionWrapper('blockWDir',block,[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
pfw["objects"]=[blockWDir]


# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.329
c11 = 165.64
c12 = 63.94
c44 = 79.51

pfw["materials"] = ["silicon"]
pfw["materialPropertyString"]="""
<ElasticCubic
	name="silicon"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultC11="""+'"'+str(c11)+'"'+"""
	defaultC12="""+'"'+str(c12)+'"'+"""
    defaultC44="""+'"'+str(c44)+'"'+"""/>
"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[0, 0, 2, 2, 0, 0]   

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[
[0,  1.0, 1.0,	1.0],
[10, 1.0, 1.10, 1.0]]
