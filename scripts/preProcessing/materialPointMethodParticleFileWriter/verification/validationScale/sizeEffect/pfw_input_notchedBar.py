# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = False
# GEOMETRY ----------------------------------------------------------------------------
D = 1 # sample thickness (mm) this is "D" int he paper"

sampleX = 3.5*D    # ratio of weibull referene volume to sample volume
sampleY = 1.0*D
sampleZ = 1.0*D
crackLength = 0.2*D  # this is for the 1/5th

refine = 1

# MESH ----------------------------------------------------------------------------
cpvc=5               # average number of grid cells across voronoi cell.
ppc=2                # particles per cell in each direction

# MATERIAL ----------------------------------------------------------------------------
# Will be read from input file for plotting.
density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8e-32
thirdInvariantDependence=1
damagedMaterialFrictionSlope=0.2

weibullSeed = 1
weibullModulus = 10.0
weibullVolume = 1.0

enableEnergyFailureCriterion = 1
fractureToughness = 1.591*(.001)*np.sqrt(1000.)     # fracture toughness, KIc:  GPa mm^1/2
constrainedModulus = bulk + 4./3.*shear
fractureEnergyReleaseRate = fractureToughness * fractureToughness / constrainedModulus

# LOADING ----------------------------------------------------------------------------
waveSpeed = np.sqrt( (bulk+1.333*shear)/density )
stopTime = 1000.*D/waveSpeed
# Indicate whether an F table is being used to apply boundary conditions or not. 0 = no, 1 = yes


# Domain ---------------------------------------------------------------------------------
domainX = 1.4*sampleX
domainY = 1.4*sampleY
domainZ = 1.4*sampleZ

cppx = 8
cppy = 8
cppz = 8

pfw["xpar"]=7*refine
pfw["ypar"]=2*refine
pfw["zpar"]=2*refine

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] =-0.5*domainX # mm
pfw["xmax"] = 0.5*domainX # mm
pfw["ymin"] =-0.5*domainY # mm
pfw["ymax"] = 0.5*domainY # mm
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

DX = domainX/pfw["nI"]

# Global roller/notch gap used by getDamage and make_objects.
gap = 0.5*(domainY-sampleY)


# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# END BATCH PARAMETERS ---------------------------------------------------------------
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 100
pfw["restartInterval"] = stopTime*5.0

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------
pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1
pfw["separabilityMinDamage"]=0.5
pfw["maxSingleFieldStateFractionForSeparability"]=0.999
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=1
pfw["plotUnscaledParticles"]=1
pfw["frictionCoefficient"]=0.05
pfw["updateMethod"]="FMPM"
pfw["updateOrder"]="2"
pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "StrengthScale",
                             "Damage",
                             "RVector"]

# END GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------
# GEOMETRY OBJECTS ------------------------------------------------------------------------------
# function to define crack:
def getDamage(self,pt):
  damage = 0.
  if ( ( -0.5*DX < pt[0] and pt[0] < 0.5*DX ) and pt[1] < pfw["ymin"]+gap+crackLength ):
    damage = 1.
  return damage

def make_objects():
  bar = geom.box('box',
                 x0 = [ -0.5*sampleX, -0.5*sampleY, -0.5*sampleZ ],
                 x1 = [  0.5*sampleX,  0.5*sampleY,  0.5*sampleZ ],
                 vel=[0.,0.,0],
                 mat=0,
                 flaggedSurfaces=[True,True,True,True,True,True],
                 group=0)

  notchedBar = geom.damageWrapper('notchedBar',subObject=bar, damage=getDamage)

  DX = domainX/pfw["nI"]
  weibullFlawSize = 6.0*DX
  weibullBar = geom.voronoiWeibullBoxWrapper('weibullSubstrate',
    subObject=notchedBar,
    x0 = np.array( [-0.5*sampleX - weibullFlawSize, -0.5*sampleY - weibullFlawSize, -0.5*sampleZ - weibullFlawSize ] ),
    x1 = np.array( [ 0.5*sampleX + weibullFlawSize, 0.5*sampleY + weibullFlawSize, 0.5*sampleZ + weibullFlawSize ] ),
    flawSize=weibullFlawSize,
    weibullVolume=weibullVolume,
    weibullModulus=weibullModulus,
    weibullSeed=weibullSeed,
    vMin=(DX)**3.,
    vpts=None,
    dim=3,
    randomMatDir = False
    )


  leftRoller = geom.cylinder("left",
                            x1=[ -1.25*D, pfw["ymin"], pfw["zmin"] ],
                            x2=[ -1.25*D, pfw["ymin"], pfw["zmax"] ],
                            r=gap,
                            vel=[0.,0.,0.],
                            mat=1,
                            group=0)
  rightRoller = geom.cylinder("left",
                            x1=[ 1.25*D, pfw["ymin"], pfw["zmin"] ],
                            x2=[ 1.25*D, pfw["ymin"], pfw["zmax"] ],
                            r=gap,
                            vel=[0.,0.,0.],
                            mat=1,
                            group=0)
  topRoller = geom.cylinder("left",
                            x1=[ 0., pfw["ymax"], pfw["zmin"] ],
                            x2=[ 0., pfw["ymax"], pfw["zmax"] ],
                            r=gap,
                            vel=[0.,0.,0.],
                            mat=1,
                            group=0)

  objects = [weibullBar,leftRoller,rightRoller,topRoller]
  return objects


# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "ceramic", "steel" ]
pfw["materialPropertyString"]="""
<CeramicDamage
	name="ceramic"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
	compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
	maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
	thirdInvariantDependence="""+'"'+str(thirdInvariantDependence)+'"'+"""
  damagedMaterialFrictionSlope="""+'"'+str(damagedMaterialFrictionSlope)+'"'+"""
  fractureToughness="""+'"'+str(fractureToughness)+'"'+"""
  enableEnergyFailureCriterion="""+'"'+str(enableEnergyFailureCriterion)+'"'+"""
  fractureEnergyReleaseRate="""+'"'+str(fractureEnergyReleaseRate)+'"'+"""/>

 <VonMisesJ
  name="steel"
  defaultDensity="8.0"
  defaultBulkModulus="160."
  defaultShearModulus="80."
  defaultYieldStrength="0.20"
  />
"""


# Deformation ---------------------------------------------------------------------------------

pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 0, 0 ]

pfw["fTableInterpType"] = "Smoothstep"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,        1.00, 1.00, 1.00],
		       [stopTime, 1.00, 0.95, 1.00]]
