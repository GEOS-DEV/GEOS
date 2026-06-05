import pfw_geometryObjects as geom
import numpy as np
from sklearn.neighbors import KDTree

pfw = {}
pfw["runDebug"] = False
stopTime = 10

# DOMAIN ---------------------------------------------------------------------------------

refine = 5 # partitions per direction
cpp = 10   # cells per partition in each direction

domainHeight = 2.0
domainWidth = 2.0
domainLength = 2.0

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth  # mm
pfw["ymin"] = 0.0              # mm
pfw["ymax"] = domainHeight     # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

pfw["planeStrain"] = 0
pfw["periodic"] = [False, False, False]

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

ppvc=10               # average number of particles across voronoi cell.

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8
thirdInvariantDependence = 1
damagedMaterialFrictionalSlope = 0.10

pfw["materials"] = ["sand"] # Only include constitutive models for the particles in this list in the order they appear below
# Include both particle and cohesive zone models here
pfw["materialPropertyString"] = f"""
<CeramicDamage
    name="sand"
    defaultDensity="{density}"
    defaultBulkModulus="{bulk}"
    defaultShearModulus="{shear}"
    tensileStrength="{tensileStrength}"
    compressiveStrength="{compressiveStrength}"
    maximumStrength="{maximumStrength}"
    crackSpeed="{crackSpeed}"
    thirdInvariantDependence="{thirdInvariantDependence}"/>
<PolymerCohesiveZone
    name="cz1"
    thickness="0.01"
    bulkModulus="1.0"
    bulkModulusA="0.0"
    bulkModulusB="1.0"
    bulkModulusT0="0.0"
    shearModulus="0.25"
    shearModulusA="0.0"
    shearModulusB="1.0"
    shearModulusT0="0.0"
    yieldStrength0="0.005"
    yieldStrengthA="0.0"
    yieldStrengthB="1.0"
    yieldStrengthT0="0.0"
    r0="0.025"
    r0A="0.0"
    r0B="1.0"
    r0T0="0.0"
    r1="0.01"
    r2="0.15"
    Gr="0.002"
    GrA="0.0"
    GrB="1.0"
    GrT0="0.0"
    maximumStretch="3.0"
    maximumStretchA="0.0"
    maximumStretchB="1.0"
    maximumStretchT0="0.0"/>
"""

# BATCH PARAMETERS  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "01:00:00"
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20
pfw["solverProfiling"]=0

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["needsNeighborList"]=1
pfw["damageFieldPartitioning"]=1
pfw["treatFullyDamagedAsSingleField"]=1
pfw["resetDefGradForFullyDamagedParticles"]=1
pfw["plotUnscaledParticles"]=1
pfw["FSubcycles"]=10

pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["frictionCoefficient"]=0.01
pfw["contactGapCorrection"] = "Implicit"
pfw["explicitSurfaceNormalInfluence"]= 1000.0
pfw["useSurfacePositionForContact"]= 1
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"]=1

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

pfw["useEvents"]=1

pfw["preventCZInterpenetration"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
							 "StrengthScale",
                             "SurfaceFlag",
                             "RVector",
							 "SurfaceNormal",
							 "SurfacePosition"]

# GEOMETRY OBJECTS -------------------------------------------------------
grainDiameter = 0.225
radialBias = 1.0
seed = 1436359
porosity = 0.0

dx = domainHeight/(pfw["nI"]-2)/pfw["ppc"]
dy = domainWidth/(pfw["nJ"]-2)/pfw["ppc"]
dz = domainWidth/(pfw["nK"]-2)/pfw["ppc"]

# Weibull variability is applied per each granule
vMin=dx*dy*dz
flawSize = dx*ppvc     # average diameter of voronoi cell.

weibullModulus = 5
weibullVolume = 100

radius = 0.75
height = domainHeight

def make_objects():
    tesselation=geom.czPrill("tesselation", [-radius, 0.0, -radius], [radius, height, radius], grainDiameter=grainDiameter, porosity=porosity, radialBias=radialBias, seed=seed, vel=[0.0,0.0,0.0], mat=0, group=0, dim=3)
    cylinder=geom.cylinder("cylinder", [0.0, 0.0, 0.0], [0.0, height, 0.0], radius, vel=[0.0,0.0,0.0], mat=0, group=0)
    czCylinder=geom.intersection("czCylinder", tesselation, cylinder)
    czCylinderWStrength = geom.voronoiWeibullBoxWrapper("prill1WStrength", czCylinder, [-radius, 0.0, -radius], [radius, height, radius], flawSize, weibullModulus, weibullVolume, seed, vMin, vpts=tesselation.vpts, dim=3)

    return [czCylinderWStrength]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 0, 0 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,            1.0,  1.0,    1.0],
               [stopTime,     1.0,  0.9,    1.0]]

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
