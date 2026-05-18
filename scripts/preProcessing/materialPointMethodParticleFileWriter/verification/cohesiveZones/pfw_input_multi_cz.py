# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = False
stopTime = 20.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.46 # mass density mg/mm^3
K = 38.67 # bulk modulus (GPa)
G = 29.0 # shear modulus (GPa)

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

refine=5 # grid partitions
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

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime*100 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageWriteInterval"]=stopTime/2000
pfw["frictionCoefficient"]=0.25

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

pfw["useEvents"]=1
pfw["plotGridFields"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "CZTag",
                             "RVector",
                             "SurfaceNormal",
                             "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------
numBoxes = 3
box_dy = sampleHeight/numBoxes 

def tag(pt):
    if pt[1] < 0.0:
        return 0
    else:
        return 1   

pfw["objects"] = []
for b in range(numBoxes):
    x0 = np.array([-sampleWidth/2, -sampleHeight/2, pfw["zmin"]]) + np.array([0.0, b*box_dy, 0.0])
    x1 = x0 + np.array([sampleWidth, box_dy, domainLength])
    box=geom.box('box',x0,x1,vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, True, False, True])
    boxWSurfFlag = geom.surfaceFlagWrapper('boxWFlag', box, 3)
    boxWTag = geom.czTagWrapper('boxWTag', boxWSurfFlag, tag )
    pfw["objects"].append(boxWTag)

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.00,  1.30,  1.00]]

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"]="03:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/112.))
pfw["mSubmitJobs"]=False
# pfw["autoRestart"]=True

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"]="""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultBulkModulus=""" + '"' + str(K) + '"' + """
    defaultShearModulus=""" + '"' + str(G) + '"' + """/>
<CoupledCohesiveZone
    name="cz1"
    defaultMaxNormalStress="0.1"
    defaultMaxShearStress="0.1"
    characteristicNormalDisplacement="0.05"
    characteristicTangentialDisplacement="0.05"
    maxTangentialDisplacement="0.1"
    maxNormalDisplacement="0.1"/>
<CoupledCohesiveZone
    name="cz2"
    defaultMaxNormalStress="0.05"
    defaultMaxShearStress="0.05"
    characteristicNormalDisplacement="0.05"
    characteristicTangentialDisplacement="0.05"
    maxTangentialDisplacement="0.1"
    maxNormalDisplacement="0.1"/>
"""

pfw["mpmEventsString"]="""
<CohesiveZone
    name="czEvent" 
    startTime=
""" + '"' + str(0.0) + '"' + """
    regionNames="{cz1, cz2}"
    constitutiveModels="{cz1, cz2}"
    czTags="{0, 1}"
    czVolumeNormalization="1"/>
"""

