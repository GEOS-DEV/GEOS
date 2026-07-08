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

density = 2.00 # mass density mg/mm^3
E = 10.0   # stiffness (GPa)
nu = 0.17   # poisson's ratio (-)

# Domain ---------------------------------------------------------------------------------
sampleWidth = 1.0
sampleHeight = 0.5

domainWidth = sampleWidth
domainHeight = 2*sampleHeight

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [False, False, False]

cpp = 12 # cells per partition in each direction

refine = 1
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

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
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
pfw["frictionCoefficient"]=0.0

# pfw["plotUnscaledParticles"]=1
# pfw["overlapCorrection"]=2
# pfw["overlapThreshold1"]=1.05
# pfw["overlapThreshold2"]=1.10

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

pfw["enablePrescribedBoundaryTransverseVelocities"]=[1,1,1,1,1,1]
pfw["prescribedBoundaryTransverseVelocities"]=[[0.0, 0.0],
                                               [0.0, 0.0],
                                               [0.0, 0.0],
                                               [0.0, 0.0],
                                               [0.0, 0.0],
                                               [0.0, 0.0]]


pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"]=1
pfw["explicitSurfaceNormalInfluence"]=141.4213562373095

pfw["useEvents"]=1
pfw["plotGridFields"] = 1
pfw["preventCZInterpenetration"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
							 "SurfaceNormal",
							 "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

gap = sampleWidth/4

box1=geom.box('box1',[-sampleWidth/2+gap, -sampleHeight/2, pfw["zmin"]],[sampleWidth/2, 0.0, pfw["zmax"]],vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2,flaggedSurfaces=[False, False, False, True])
box2=geom.box('box2',[-sampleWidth/2, 0.0, pfw["zmin"]],[sampleWidth/2-gap, sampleHeight/2, pfw["zmax"]],vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, True, False, False])

boxWFlag1 = geom.surfaceFlagBoxWrapper('boxWFlag1', [-sampleWidth/2+gap, -sampleHeight/4,  pfw["zmin"]], [sampleWidth/2-gap, sampleHeight/4, pfw["zmax"]], 3, box1)
boxWFlag2 = geom.surfaceFlagBoxWrapper('boxWFlag2', [-sampleWidth/2+gap, -sampleHeight/4,  pfw["zmin"]], [sampleWidth/2-gap, sampleHeight/4, pfw["zmax"]], 3, box2)

pfw["objects"]=[boxWFlag1, boxWFlag2]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 0, 0, 1, 1 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.25,  1.00,  1.00]]

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"] = f"""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity="{density}"
    defaultYoungModulus="{E}"
    defaultPoissonRatio="{nu}"/>
<CoupledCohesiveZone
    name="cz"
    defaultMaxNormalStress="0.1"
    defaultMaxShearStress="0.1"
    characteristicNormalDisplacement="0.1"
    characteristicTangentialDisplacement="0.1"
    maxNormalDisplacement="0.1"
    maxTangentialDisplacement="0.1"/>"""

pfw["cohesiveZoneRegions"] = """
<CohesiveZoneRegion
    name="cz"
    constitutiveModel="cz"
    tag="0"/>"""

pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="cz"
    regionNames="{cz}"
    czVolumeNormalization="1"
    startTime="0.0"/>
"""
