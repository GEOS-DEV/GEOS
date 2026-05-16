# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = True
stopTime = 20.0

# Domain ---------------------------------------------------------------------------------

sampleWidth = 1.0
sampleHeight = 1.0

domainWidth = 1.25*sampleWidth
domainHeight = sampleHeight

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [False, False, False]

refine = 2 # grid partitions
cpp = 10 # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3                # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

domainLength = domainHeight/(pfw["nJ"]-2)

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

dy_p = domainHeight/( pfw["ppc"]*(pfw["nJ"]-2) ) # Particle height
dy_p = domainHeight/( pfw["ppc"]*(pfw["nJ"]-2) ) # Particle height

overlapStrain = dy_p / ( 2 * domainHeight ) # strain to overlap half the particle height of the cohesive zone

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0 # Explicitly defined particle fields (e.g. 0 and 1)

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageWriteInterval"]=stopTime/2000
pfw["frictionCoefficient"]=0.25

pfw["plotGridFields"]=1

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["enableContact"]=0 # Disabled to allow cz to interpentrate during initialization
pfw["contactGapCorrection"]="Implicit"
pfw["explicitSurfaceNormalInfluence"]= 1000
pfw["useSurfacePositionForContact"]= 1
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"]=1

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

pfw["useEvents"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
                             "SurfaceNormal",
                             "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

box1=geom.box('box1',[-sampleWidth/2, -sampleHeight/2, pfw["zmin"]],[sampleWidth/2, 0.0, pfw["zmax"]],vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, True])
box2=geom.box('box2',[-sampleWidth/2, 0.0, pfw["zmin"]],[sampleWidth/2, sampleHeight/2, pfw["zmax"]],vel=[0.0, 0.0, 0.0], mat=0, group=1, dim=2, flaggedSurfaces=[False, True, False, False])

boxWFlag1 = geom.surfaceFlagWrapper('boxWFlag1',box1, 3)
boxWFlag2 = geom.surfaceFlagWrapper('boxWFlag2',box2, 3)

pfw["objects"]=[boxWFlag1,boxWFlag2]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,        1.00,  1.00,  1.00],
               [stopTime/4, 1.00,  1.00 - overlapStrain,  1.00],
               [stopTime/2, 1.00,  1.00 - overlapStrain , 1.00],
               [stopTime,   1.00,  1.10, 1.00]]

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBank"]="imcomp"
pfw["mBatch"]=True
pfw["mWallTime"]="03:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.))
pfw["mSubmitJobs"]=False
# pfw["autoRestart"]=False

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.46 # mass density mg/mm^3
K = 38.67 # bulk modulus (GPa)
G = 29.0 # shear modulus (GPa)

pfw["materials"] = ["elasticIsotropic"] # Include only the particle constitutive models
pfw["materialPropertyString"]="""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultBulkModulus=""" + '"' + str(K) + '"' + """
    defaultShearModulus=""" + '"' + str(G) + '"' + """/>
<CoupledCohesiveZone
    name="cz1"
    maxNormalStress="0.1"
    maxShearStress="0.1"
    characteristicNormalDisplacement="0.05"
    characteristicTangentialDisplacement="0.05"
    maxTangentialDisplacement="0.1"
    maxNormalDisplacement="0.1"/>"""

# EVENTS ------------------------------------------------------------------------------------

pfw["mpmEventsString"]="""
<CohesiveZone
    name="cz1" 
    startTime=""" + '"' + str(stopTime/4) + '"' + """
    endTime=""" + '"' + str(stopTime) + '"' + """
    constitutiveModel="cz1"
    computeNormalsAndPositions=""" + '"' + str(0) + '"' + """
    normalsAndPositionsMethod=""" + '"' + str("LogisticRegression") + '"' + """/>
"""
