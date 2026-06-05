# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = False
stopTime = 40.0

# Domain ---------------------------------------------------------------------------------

sampleWidth = 0.5
sampleHeight = 1.0

domainWidth = 1.5
domainHeight = 1.5

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [False, False, False]

refine = 1 # grid partitions
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
pfw["restartInterval"]=stopTime*100 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["plotUnscaledParticles"]=1
pfw["plotGridFields"]=1

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["updateMethod"]="FLIP"
# pfw["updateOrder"]=2

pfw["useEvents"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
                             "MaterialDirection",
                             "SurfaceNormal",
                             "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Radially varying velocity (spinning disk)
maxVelocity = 0.05
def getVelocity(self, pt):
    pt = np.array(pt)
    norm = np.linalg.norm(pt)
    pt = pt / norm
    vel = maxVelocity * ( norm / ( sampleHeight / 2 ) ) * np.cross(np.array(pt),np.array([0,0,1]))
    return vel

box1=geom.box('box1',[-sampleWidth/2, -sampleHeight/2, pfw["zmin"]],[sampleWidth/2, 0.0, pfw["zmax"]],vel=getVelocity, mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, True])
box2=geom.box('box2',[-sampleWidth/2, 0.0, pfw["zmin"]],[sampleWidth/2, sampleHeight/2, pfw["zmax"]],vel=getVelocity, mat=0, group=1, dim=2, flaggedSurfaces=[False, True, False, False])

angle = np.deg2rad(10.0)

boxWMatDir1 = geom.materialDirectionWrapper('boxWMatDir1', box1, np.array([1.0, 0.0, 0.0]) )
boxWMatDir2 = geom.materialDirectionWrapper('boxWMatDir2', box2, np.array([np.cos(angle), np.sin(angle), 0.0]))

boxWFlag1 = geom.surfaceFlagWrapper('boxWFlag1',boxWMatDir1, 3)
boxWFlag2 = geom.surfaceFlagWrapper('boxWFlag2',boxWMatDir2, 3)

pfw["objects"]=[boxWFlag1, boxWFlag2]

# BOUNDARY CONDITIONS_---------------------------------------------------------------------

pfw["boundaryConditionTypes"]=[ 0, 0, 0, 0, 1, 1 ]

pfw["prescribedBoundaryFTable"]=0

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.46 # mass density mg/mm^3
K = 38.67 # bulk modulus (GPa)
G = 29.0 # shear modulus (GPa)

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"] = f"""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity="{density}"
    defaultBulkModulus="{K}"
    defaultShearModulus="{G}"/>
<CoupledCohesiveZone
    name="cz"
    defaultMaxNormalStress="2.0"
    defaultMaxShearStress="2.0"
    characteristicNormalDisplacement="0.05"
    characteristicTangentialDisplacement="0.05"/>"""

pfw["cohesiveZoneRegions"] = """
<CohesiveZoneRegion
    name="cz"
    constitutiveModel="cz"
    tag="0"/>"""

# EVENTS ------------------------------------------------------------------------------------

pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="cz"
    regionNames="{cz}"
    czVolumeNormalization="1"
    startTime="0.0"/>
"""
