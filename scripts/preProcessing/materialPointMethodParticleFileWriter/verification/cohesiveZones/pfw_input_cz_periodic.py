# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = True
stopTime = 20.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.46 # mass density mg/mm^3
K = 38.67 # bulk modulus (GPa)
G = 29.0 # shear modulus (GPa)

# Domain ---------------------------------------------------------------------------------

sampleWidth = 1.0
sampleHeight = 1.0

domainWidth = sampleWidth
domainHeight = sampleHeight

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [True, False, False]

refine=3 # grid partitions
cpp = 11 # cells per partition in each direction

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
pfw["restartInterval"]=stopTime*2

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

pfw["preventCZInterpenetration"]=1

pfw["useEvents"]=1
pfw["plotGridFields"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
                             "SurfaceNormal",
                             "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

box=geom.box('box',[-domainWidth/2, -domainHeight/2, pfw["zmin"]],[domainWidth/2, domainHeight/2, pfw["zmax"]],vel=[0.1, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, False])

interface1 = geom.surfaceWrapper('interface1', box, [0.0, -domainHeight/2, 0.0], [0.0, domainHeight/2, 0.0], 3)
interface2 = geom.surfaceWrapper('interface2', interface1, [-domainWidth/2, 0.0, 0.0], [-0.1, 0.0, 0.0], 3)
interface3 = geom.surfaceWrapper('interface2', interface2, [0.1, 0.0, 0.0], [domainWidth/2, 0.0, 0.0], 3)

pfw["objects"]=[interface3]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.00,  1.20,  1.00]]

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBank"]="MAHEM"
pfw["mBatch"]=True
pfw["mWallTime"]="00:30:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"]=f"""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity="{density}"
    defaultBulkModulus="{K}"
    defaultShearModulus="{G}"/>
<CoupledCohesiveZone
    name="coupled"
    defaultMaxNormalStress="0.1"
    defaultMaxShearStress="0.1"
    characteristicNormalDisplacement="0.05"
    characteristicTangentialDisplacement="0.05"
    maxTangentialDisplacement="0.1"
    maxNormalDisplacement="0.1"/>"""

pfw["cohesiveZoneRegions"]="""
<CohesiveZoneRegion
    name="coupled"
    constitutiveModel="coupled"
    tag="0"/>
"""

pfw["mpmEventsString"]="""
<CohesiveZone
    name="cz"
    regionNames="{coupled}"
    startTime="0.0"
    czVolumeNormalization="1"
    computeNormalsAndPositions="0"
    normalsAndPositionsMethod="LogisticRegression"/>
"""

