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
domainLength = 1.0

sampleWidth = 0.5*domainWidth
sampleHeight = 2.0*domainHeight
sampleLength = 0.5*domainLength

pfw["xmin"] = -0.5*domainWidth  # mm
pfw["xmax"] =  0.5*domainWidth  # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

pfw["planeStrain"] = 0

pfw["periodic"] = [False, False, False]

refine = 5 # grid partitions
cpp = 8 # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2         # particles per cell in each direction

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
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

# Define all the geometric objects -------------------------------------------------------

box1=geom.box('box1',[-sampleWidth/2, -sampleHeight, -sampleLength/2], [sampleWidth/2, 0.0,          sampleLength/2], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=3, flaggedSurfaces=[False, False, False, False, True, False])
box2=geom.box('box2',[-sampleWidth/2,  0.0,          -sampleLength/2], [sampleWidth/2, sampleHeight, sampleLength/2], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=3, flaggedSurfaces=[False, True, False, False, False, False])

def transform_y_to_axis(axis, origin=(0.0, 0.0, 0.0)):
    """
    Create a 4x4 homogeneous transform whose local +Y direction aligns
    with the given 3D axis vector.

    Convention:
        column vectors, i.e. p_world = T @ p_local_homogeneous

    Therefore:
        T[:3, :3] @ [0, 1, 0] == normalized(axis)

    Parameters
    ----------
    axis : array-like, shape (3,)
        Target direction for the local +Y axis.
    origin : array-like, shape (3,), optional
        Translation part of the transform.

    Returns
    -------
    T : ndarray, shape (4, 4)
        Homogeneous transformation matrix.
    """

    axis = np.asarray(axis, dtype=float)
    origin = np.asarray(origin, dtype=float)

    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        raise ValueError("axis must be a non-zero 3D vector")

    target = axis / norm
    source = np.array([0.0, 1.0, 0.0])

    c = np.dot(source, target)

    # Already aligned
    if np.isclose(c, 1.0):
        R = np.eye(3)

    # Opposite direction: rotate 180 degrees around X
    elif np.isclose(c, -1.0):
        R = np.array([
            [1.0,  0.0,  0.0],
            [0.0, -1.0,  0.0],
            [0.0,  0.0, -1.0],
        ])

    else:
        v = np.cross(source, target)
        s = np.linalg.norm(v)

        K = np.array([
            [0.0,   -v[2],  v[1]],
            [v[2],   0.0, -v[0]],
            [-v[1],  v[0],  0.0],
        ])

        # Rodrigues rotation formula
        R = np.eye(3) + K + K @ K * ((1.0 - c) / (s ** 2))

    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = origin

    return T

M = transform_y_to_axis([1.0, 1.0, 1.0])
rotatedBox1 = geom.transform("rotatedBox1", box1, M)
rotatedBox2 = geom.transform("rotatedBox2", box2, M)

boxWFlag1 = geom.surfaceFlagWrapper('boxWFlag1', rotatedBox1, 3)
boxWFlag2 = geom.surfaceFlagWrapper('boxWFlag2', rotatedBox2, 3)

pfw["objects"]=[boxWFlag1,boxWFlag2]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.0 + 0.15/np.sqrt(3),  1.0 + 0.15/np.sqrt(3),  1.0 + 0.15/np.sqrt(3)]]

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:15:00"
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=True

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
