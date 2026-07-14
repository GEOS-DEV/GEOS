# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 500.0

# Domain ---------------------------------------------------------------------------------

refine=1

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

# We define this to be a single cell in the loading direction
# so we aren't affected by localization or shear banding.

pfw["nI"]=5     # grid cells in the x-direction
pfw["nJ"]=3     # grid cells in the y-direction
pfw["nK"]=5     # grid cells in the z-direction
pfw["ppcx"]=1   # particles per cell in each direction
pfw["ppcy"]=1   # particles per cell in each direction
pfw["ppcz"]=1   # particles per cell in each direction

# Define all the geometric objects -------------------------------------------------------
sampleWidth = 0.1  # mm
sampleHeight = 0.1 # mm
sampleLength = 0.1 # mm

domainWidth = 2.0*sampleWidth  # This would be increased for unconfined compression.
domainHeight = sampleHeight
domainLength = 2.0*sampleLength

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] =  -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:05:00"
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/1000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/1000
pfw["useInternalForceAsFaceReaction"]=0

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1   # Would need to be zo for large stretch tension test if not single element or high ppc
pfw["damageFieldPartitioning"]=0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["useEvents"]=1
pfw["frictionCoefficient"]=0.0

pfw["updateMethod"]="PIC"


# Deformation ---------------------------------------------------------------------------------
# Prescribed stretch in the loading direction, outflow boundary laterally to allow free 
# contraction or expansion.  There is a significant gap so the material shouldn't
# interact with lateral boundary.
pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 0, 0 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,            1.0,   1.0,   1.0],
			   [0.8*stopTime, 1.0,   0.606, 1.0],
			   [stopTime,     1.0,   0.637, 1.0]]

# Define all the geometric objects -------------------------------------------------------

block = geom.box('block',[-sampleWidth/2, pfw["ymin"] ,-sampleLength/2],[sampleWidth/2,pfw["ymax"] ,sampleLength/2], vel=[0.0,0.0,0.0], mat=0, group=0)
pfw["objects"]=[block]

# Generic FKM/Viton-like parameters for this legacy StrainHardeningPolymer
# validation input.  These are deliberately generic validation values and should
# be replaced by a project material card for material-specific studies.
viton_density = 1.85
viton_young_modulus = 0.01577
viton_poisson_ratio = 0.49
viton_bulk_modulus = viton_young_modulus / (3.0 * (1.0 - 2.0 * viton_poisson_ratio))
viton_shear_modulus = viton_young_modulus / (2.0 * (1.0 + viton_poisson_ratio))
viton_yield_strength = 0.0030
viton_softening_magnitude = 0.0030
viton_softening_r1 = 0.30
viton_softening_r2 = 1.25
viton_hardening_slope = 0.0020
viton_maximum_stretch = 2.60

pfw["materials"] = [ "polymer" ]
pfw["materialPropertyString"]="""
<StrainHardeningPolymer
    name="polymer"
    defaultDensity="{density}"
    defaultBulkModulus="{bulk}"
    bulkModulusA="0.0"
    bulkModulusB="0.0"
    bulkModulusT0="300.0"
    defaultShearModulus="{shear}"
    shearModulusA="0.0"
    shearModulusB="0.0"
    shearModulusT0="300.0"
    defaultYieldStrength="{yield_strength}"
    yieldStrengthA="0.0"
    yieldStrengthB="0.0"
    yieldStrengthT0="300.0"
    strainHardeningSlope="{hardening}"
    strainHardeningSlopeA="0.0"
    strainHardeningSlopeB="0.0"
    strainHardeningSlopeT0="300.0"
    shearSofteningMagnitude="{softening}"
    shearSofteningMagnitudeA="0.0"
    shearSofteningMagnitudeB="0.0"
    shearSofteningMagnitudeT0="300.0"
    shearSofteningShapeParameter1="{r1}"
    shearSofteningShapeParameter2="{r2}"
    maximumStretch="{maximum_stretch}"
    maximumStretchA="0.0"
    maximumStretchB="0.0"
    maximumStretchT0="300.0"
  />
""".format(
    density=viton_density,
    bulk=viton_bulk_modulus,
    shear=viton_shear_modulus,
    yield_strength=viton_yield_strength,
    hardening=viton_hardening_slope,
    softening=viton_softening_magnitude,
    r1=viton_softening_r1,
    r2=viton_softening_r2,
    maximum_stretch=viton_maximum_stretch,
)
