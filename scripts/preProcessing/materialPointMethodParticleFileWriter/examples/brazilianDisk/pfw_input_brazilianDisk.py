# -*- coding: utf-8 -*-
#
# 2D Brazilian Disk fragmentation under dynamic loading.
# The disk is in a rectangular domain with periodic lateral boundaries
# and an additional horizontal velocity, subject to a compression via
# F-table driven boundary motion.  
#
# This is a demonstration of many advanced capabilities and is used as a robustness
# test for new algorithms, a general demonstration of problems where MPM excels
# and for profiling performance in planeStrain with a realistic feature set.
#
# Ceramic Damage model
# Damage-field gradient partitioning
# Crack-tip correction
# Periodic boundary conditions
# FMPM with incremental contact
# Weibull variability (strength scale)
# --------------------------------------------------------------------------------------
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import math


pfw = {} 
pfw["runDebug"] = True
stopTime = 2.0


# MATERIAL PROPERTIES ----------------------------------------------------------------
# Import material properties from validation materials file --------------------------
# [pfw_dependency] /pfw_materials.py
import importlib
matFile = importlib.import_module('pfw_materials')
quartz = matFile.quartz

# Assign imported model to pfw dictionary.
pfw["materials"] = [quartz["name"]]
pfw["materialPropertyString"]=quartz["materialString"]

# Domain ---------------------------------------------------------------------------------

refine = 6  				# Partitions in each direction.
cpp = 64    				# Cells per partition
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=round(pfw["xpar"]*cpp*1.5)   	# grid cells in the x-direction

pfw["nJ"]=pfw["ypar"]*cpp 	# grid cells in the y-direction
pfw["nK"]=3  			# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainHeight = 1.0
domainWidth = 1.5*domainHeight
domainLength = domainHeight*(pfw["nK"])/(pfw["nJ"])  # m, to get cubic cells

# grid cell size
DX = domainWidth / pfw["nI"]

# Define all the geometric objects -------------------------------------------------------
pfw["xmin"] =-0.5*domainWidth	# m
pfw["xmax"] = 0.5*domainWidth	# m
pfw["ymin"] = 0.0	# m
pfw["ymax"] = domainHeight	# m
pfw["zmin"] =-0.5*domainLength 	# m
pfw["zmax"] = 0.5*domainLength 	# m

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:30:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

# Options are 'silo' (for VisIt) or 'vtk' (for Paraview)
pfw["outputType"]="silo"

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1
pfw["planeStrain"] = 1

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

pfw["solverProfiling"]=1
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["frictionCoefficient"]=0.25

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "StrengthScale",
                             "SurfaceFlag",
                             "RVector"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]
pfw["periodic"]=[ True, False, False]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,        1.00, 1.00, 1.00],
		       [stopTime, 1.00, 0.80, 1.00]]

# Define all the geometric objects -------------------------------------------------------

disk1 = geom.cylinder('disk1',[0,domainHeight/2,pfw["zmin"]],[0,domainHeight/2,pfw["zmax"]],domainHeight/2,vel=[1.0,0,0],mat=0,group=0)


# Assign a voronoi-tesselated weibull distribution of strength to the grain.

weibullFlawSize = 6*DX

weibullSample = geom.voronoiWeibullBoxWrapper('weibullSubstrate',                                                     
subObject=disk1,       
x0 = np.array( [ pfw["xmin"] - weibullFlawSize, pfw["ymin"] - weibullFlawSize, pfw["zmin"] - weibullFlawSize ] ),
x1 = np.array( [ pfw["xmax"] + weibullFlawSize, pfw["ymax"] + weibullFlawSize, pfw["zmax"] + weibullFlawSize ] ),
flawSize=weibullFlawSize,
weibullVolume=quartz["weibullReferenceVolume"],
weibullModulus=quartz["weibullModulus"] ,
weibullSeed=1,
vMin=(DX)**3.,
vpts=None,
dim=3,
randomMatDir = False
)
    

pfw["objects"]=[weibullSample]
