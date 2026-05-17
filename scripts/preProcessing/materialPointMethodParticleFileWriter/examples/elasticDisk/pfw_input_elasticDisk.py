# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
# [pfw_dependency] pfw_materials.py
import importlib
matdb = importlib.import_module('pfw_materials')

pfw = {} 
pfw["runDebug"] = True
stopTime = 20

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.93 # mass density mg/mm^3
E = 45.0     # in-plane stiffness (GPa)
nu = 0.27    # in plane poisson's ratio

# Domain ---------------------------------------------------------------------------------

pfw["planeStrain"] = 1

pfw["periodic"] = [False, False, False]

refine = 4
cpp = 10 # cells per partition in each direction
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3          # grid cells in the z-direction
pfw["ppc"]=2         # particles per cell in each direction

domainWidth = 1.0
domainHeight = domainWidth
domainLength = domainHeight*(pfw["nK"]-2)/(pfw["nJ"]-2)  # m, to get cubic cells

# Define all the geometric objects -------------------------------------------------------

pfw["xmin"] =-0.5*domainWidth   # m
pfw["xmax"] = 0.5*domainWidth   # m
pfw["ymin"] = 0.0   # m
pfw["ymax"] = domainHeight  # m
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

pfw["outputType"]="silo"

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime             
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime*10 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["frictionCoefficient"]=0.25

# pfw["plotUnscaledParticles"]=1

# pfw["overlapCorrection"]=2
# pfw["overlapThreshold1"]=1.05
# pfw["overlapThreshold2"]=1.10

pfw["maxParticleVelocity"]=10.0               
pfw["minParticleJacobian"]=0.01                       
pfw["maxParticleJacobian"]=10.0  

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------
pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 1, 1 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,            1.0,  1.0,    1.0],
               [stopTime,     1.0,  0.5,    1.0]]

# Define all the geometric objects -------------------------------------------------------

disk=geom.cylinder('disk',[0.0,domainHeight/2,pfw["zmin"]],[0.0,domainHeight/2,pfw["zmax"]],0.5*domainHeight,vel=[0.,0.,0.],mat=0,group=0)
pfw["objects"]=[disk]

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = [ matdb.elasticDemo["name"] ]
pfw["materialPropertyString"] = matdb.elasticDemo["materialString"]
