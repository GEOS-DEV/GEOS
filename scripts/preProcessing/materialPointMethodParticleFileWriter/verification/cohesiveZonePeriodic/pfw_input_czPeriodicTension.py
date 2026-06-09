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

domainWidth = 1.0
domainHeight = 1.0

sampleWidth = domainWidth/2.0
sampleHeight = 2*domainHeight

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [True, False, False]

refine = 5 # grid partitions
cpp = 8  # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3          # grid cells in the z-direction
pfw["ppc"]=2         # particles per cell in each direction

domainLength = domainHeight/(pfw["nJ"]-2)

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200

# Silo output is required by the verification-suite VisIt smoke renderer.
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity"]
pfw["restartInterval"]=stopTime

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

# pfw["plotUnscaledParticles"]=1
# pfw["overlapCorrection"]=2
# pfw["overlapThreshold1"]=1.05
# pfw["overlapThreshold2"]=1.10

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
                             "RVector",
                             "SurfaceNormal",
                             "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]

x_vel = 0.1
rel_vel = 0.0
pfw["enablePrescribedBoundaryTransverseVelocities"]=[0,0,1,1,0,0]
pfw["prescribedBoundaryTransverseVelocities"]=[[0.0, 0.0],
                                               [0.0, 0.0],
                                               [0.0, x_vel+rel_vel/2],
                                               [0.0, x_vel-rel_vel/2],
                                               [0.0, 0.0],
                                               [0.0, 0.0]]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.00,  1.15,  1.00]]


# Define all the geometric objects -------------------------------------------------------

box1=geom.box('box1',[-sampleWidth/2, -2*sampleHeight/2, pfw["zmin"]],[sampleWidth/2, 0.0, pfw["zmax"]], vel=[x_vel, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, True])
box2=geom.box('box2',[-sampleWidth/2, 0.0, pfw["zmin"]], [sampleWidth/2, 2*sampleHeight/2, pfw["zmax"]], vel=[x_vel, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, True, False, False])

angle=-45
M = geom.rotate(np.deg2rad(angle))
rotatedBox1 = geom.transform("rotatedBox1", box1, M)
rotatedBox2 = geom.transform("rotatedBox2", box2, M)

boxWFlag1 = geom.surfaceFlagWrapper('boxWFlag1', rotatedBox1, 3)
boxWFlag2 = geom.surfaceFlagWrapper('boxWFlag2', rotatedBox2, 3)

geometry = [boxWFlag1,boxWFlag2]

pfw["objects"]=[boxWFlag1,boxWFlag2]

# Add periodic images shifted by +- domainWidth
for g in geometry:
    pfw["objects"].append(geom.transform("image", g, geom.translate([-domainWidth, 0.0 , 0.0])))

for g in geometry:
    pfw["objects"].append(geom.transform("image2", g, geom.translate([domainWidth, 0.0 , 0.0])))

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:15:00"
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"] = f"""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity="{density}"
    defaultBulkModulus="{K}"
    defaultShearModulus="{G}"/>
<CoupledCohesiveZone
    name="cz1"
    defaultMaxNormalStress="0.1"
    defaultMaxShearStress="0.1"
    characteristicNormalDisplacement="0.05"
    characteristicTangentialDisplacement="0.05"
    maxTangentialDisplacement="0.1"
    maxNormalDisplacement="0.1"/>"""

pfw["cohesiveZoneRegions"] = """
<CohesiveZoneRegion
    name="cz1"
    constitutiveModel="cz1"
    tag="0"/>"""

pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="cz1"
    startTime="0.0"
    regionNames="{cz1}"
    czVolumeNormalization="1"/>
"""
