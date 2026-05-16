# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 5.

innerDiameter = 7.96 - 0.5
outerDiameter = 7.96 + 0.5
initialRadialVelocity = 0.01592


# Domain ---------------------------------------------------------------------------------
domainX = 1.1*outerDiameter
domainY = 1.1*outerDiameter

cppx=20   # cells per partition in each direction
cppy=20
cppz=5

refine=6
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] =-0.5*domainX # mm
pfw["xmax"] = 0.5*domainX # mm
pfw["ymin"] =-0.5*domainY # mm
pfw["ymax"] = 0.5*domainY # mm

DZ = ( pfw["xmax"] - pfw["xmin"] ) / ( pfw["nI"] - 2 )
domainZ = 3.*DZ
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True
pfw["mBank"]="mahem" #"MAHEM"
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
							 "StrengthScale",
							 "SurfaceFlag",
                             "RVector"]

# END GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------

# GEOMETRY OBJECTS -------------------------------------------------------------------------------------------

# Radially varying velocity for cylinder:
def getVelocity(self,pt):
    x = np.array(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    if ( r > 1.e-12):
      er = xr/r
      v = initialRadialVelocity*er
    else:
      v = np.array([0.,0.,0.])
    return v

def make_objects():
    ring = geom.cylinder('ring',
        x1=[0.0,0.0,-0.5*DZ],
        x2=[0.0,0.0,0.5*DZ],
        ri=0.5*innerDiameter,
        r=0.5*outerDiameter,
        v=getVelocity,
        mat=0,
        group=0)

    weibullring = geom.voronoiWeibullBoxWrapper('weibullring',
        subObject=ring,
        x0=[pfw["xmin"],pfw["ymin"],pfw["zmin"]],
        x1=[pfw["xmax"],pfw["ymax"],pfw["zmax"]],
        flawSize=5*DZ,
        weibullVolume=weibullVolume,
        weibullModulus=weibullModulus,
        weibullSeed=1,
        vMin=(DZ)**3.
        )

    objects=[weibullring]
    return objects

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.75
bulk = 275./( 3.*( 1. - 2.*0.3) )
shear = 275./( 2.*( 1. + 0.3) )
tensileStrength = 0.30
compressiveStrength = 8.*tensileStrength
maximumStrength = 5.0*compressiveStrength
crackSpeed = 1.8
thirdInvariantDependence=1

weibullVolume = DZ*0.25*np.pi*( outerDiameter**2 - innerDiameter**2 )
weibullModulus = 6.

pfw["materials"] = ["sand"]
pfw["materialPropertyString"]="""
<CeramicDamage
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
	compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
	maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
	crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
	thirdInvariantDependence=""" + '"' + str(thirdInvariantDependence) + '"' + """/>
"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["boundaryConditionTypes"]=[0, 0, 0, 0, 0, 0]   
