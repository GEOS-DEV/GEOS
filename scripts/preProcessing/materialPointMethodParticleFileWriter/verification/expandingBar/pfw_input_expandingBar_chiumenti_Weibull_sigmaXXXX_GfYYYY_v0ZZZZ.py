# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = False
stopTime = 5.

sampleX = 1.0
domainAspectRatio = 56

initialVelocityGradient = ZZZZ #0.01592
innerDiameter = 2.*(7.96) - 0.25
outerDiameter = 2*(7.96) + 0.25

density = 2.75						# mg/mm^3
bulk = 275./( 3.*( 1. - 2.*0.3) )	# GPa
shear = 275./( 2.*( 1. + 0.3) )		# GPa
failureStrength = XXXX #0.30				# GPa
energyReleaseRate = YYYY #1.0e-5 			# mg/us^2

weibullModulus = 6.



# Domain ---------------------------------------------------------------------------------
domainX = 1.2*sampleX
domainY = domainX/domainAspectRatio
domainZ = domainX/domainAspectRatio

cppx=14   # cells per partition in X direction
cppy=14   # cells per partition in Y direction
cppz=14   # cells per partition in Z direction

refine=1

pfw["xpar"]=refine*domainAspectRatio  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

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

DY = domainY/pfw["nJ"]
DZ = domainZ/pfw["nK"]

sampleY = domainY - 2.*DY
sampleZ = domainZ - 2.*DZ

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True
pfw["mBank"]="imcomp" #"MAHEM"
pfw["mWallTime"]="02:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
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
pfw["treatFullyDamagedAsSingleField"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=0
pfw["boxAverageHistory"]=0

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0
pfw["FSubcycles"]=10

pfw["plotUnscaledParticles"]=0

pfw["frictionCoefficient"]=0.05

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]="2"

pfw["particleFileFields"] = ["Velocity",
                    "MaterialType",
                    "ContactGroup",
                    "SurfaceFlag",
                    "StrengthScale",
                    "RVector"]

# END GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------

# GEOMETRY OBJECTS -------------------------------------------------------------------------------------------

# Linear varying velocity for cylinder:
def getVelocity(self,pt):
    x = np.array(pt)
    v = [initialVelocityGradient*x[0],0.,0.]
    return v

def make_objects():
  bar = geom.box('box', 
                 x0 = [ -0.5*sampleX, -0.5*sampleY, -0.5*sampleX ],
                 x1 = [  0.5*sampleX,  0.5*sampleY,  0.5*sampleX ],
                 vel=getVelocity,
                 mat=0,
                 group=0)
  weibullVolume = 1.0
  
  weibullBar = geom.voronoiWeibullBoxWrapper('weibullbar',
         subObject=bar,
         x0=[pfw["xmin"],pfw["ymin"],pfw["zmin"]],
         x1=[pfw["xmax"],pfw["ymax"],pfw["zmax"]],
         flawSize=4*DX,
         weibullVolume=weibullVolume,
         weibullModulus=weibullModulus,
         weibullSeed=1,
         vMin=(DX)**3.
         )
  objects=[weibullBar]
  return objects


# MATERIAL PROPERTIES --------------------------------------------------------------------



pfw["materials"] = ["sand"]
pfw["materialPropertyString"]="""
<Chiumenti
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	failureStrength="""+'"'+str(failureStrength)+'"'+"""
	energyReleaseRate="""+'"'+str(energyReleaseRate)+'"'+"""/>
"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["boundaryConditionTypes"]=[0, 0, 0, 0, 0, 0]   
