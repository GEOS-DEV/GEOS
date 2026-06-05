# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = True
stopTime = 10.0

# Domain ---------------------------------------------------------------------------------
domainWidth = 1.0
domainHeight = 1.0

sampleWidth = 1.0
sampleHeight = 0.5*domainHeight

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]

refine = 1
cpp = 12 # cells per partition in each direction

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3                # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

domainLength = domainHeight/(pfw["nJ"]-2)

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

dx=domainWidth/(pfw["nI"]-2)
dy=domainHeight/(pfw["nJ"]-2)

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

fluidPressure = 0.1

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1

pfw["lastRestartBufferInSeconds"]=600
pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

pfw["implicitFluidPressure"]=fluidPressure

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
							 "SurfaceNormal",
							 "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

wall=geom.box('wall',[-sampleWidth/2, 0.0, pfw["zmin"]],[sampleWidth/2, sampleHeight, pfw["zmax"]],vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2)
pfw["objects"]=[wall]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 1, 1 ]

pfw["prescribedFTable"]=0
pfw["prescribedBoundaryFTable"]=0

pfw["stressControl"]=[1,1,0]
pfw["stressControlKp"] = 1.0
pfw["stressControlKi"] = 0.0
pfw["stressControlKd"] = 0.0
pfw["fTableInterpType"] = "Cosine"
pfw["stressTable"]=[[0.0, -fluidPressure, -fluidPressure, -fluidPressure],
                    [stopTime, -fluidPressure, -fluidPressure, -fluidPressure]]

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False #True
pfw["lastRestartBufferInSeconds"] = 600

# GEOS MPM i/o parameters ---------------------------------------------------------------

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.00 # mass density mg/mm^3
E = 10.0   # stiffness (GPa)
nu = 0.17   # poisson's ratio (-)

pfw["materials"] = ["elasticIsotropic", "silicaNoDmg", "silica"]
pfw["materialPropertyString"]="""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultYoungModulus=""" + '"' + str(E) + '"' + """
    defaultPoissonRatio=""" + '"' + str(nu) + '"' + """/>"""
