# -*- coding: utf-8 -*-

import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True

# Unit system :  mg, microsecond, mm,  stress-> GPA, velocity mm/us=km/s

impacterDiameter = 3 # SS sphere, mm
sampleDiameter = 152.4 # target diameter
sampleLength = 152.4 # target length (cylinder) , mm
impactVelocity = 1.515 # mm/us

stopTime = sampleLength/impactVelocity 

steelDensity = 7.85
steelYield = .25*1.2 #20 percent increased for strain rate hardening
steelYoungsModulus = 210.
steelPoissonsRatio = 0.3

density = 2.4
bulk = 40./( 3.*( 1. - 2.*0.3) )
shear = 40./( 2.*( 1. + 0.3) )
tensileStrength = 0.005
compressiveStrength = 8.*tensileStrength
maximumStrength = 5.0*compressiveStrength
crackSpeed = 1.8
thirdInvariantDependence=1

weibullVolume = 10 # mm^3
weibullModulus = 18.
# Domain ---------------------------------------------------------------------------------

domainX = sampleLength+impacterDiameter #additional 0.5 to catch the ejecta 
domainY = 1*sampleDiameter #changed
domainZ = 1*sampleDiameter #changed

cppx=240   # cells per partition in each direction
cppy=240
cppz=12

refine=1

pfw["xpar"]=2*refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1  #changed from refine to 1

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
#pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["nK"]=3   # grid cells in the z-direction #changed to 3
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] = 0.0 # mm
pfw["xmax"] = domainX # mm
pfw["ymin"] = -0.5*domainY # mm
pfw["ymax"] = 0.5*domainY # mm
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

DX = (pfw["xmax"] - pfw["xmin"])/pfw["xpar"]/cppx

# BATCH PARAMETERS --------------------------------------------------------

pfw["mBatch"]=True
pfw["mBank"]="imcomp" #"MAHEM"
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False #This prevents from automatic submission of Job
pfw["autoRestart"]=False

# END BATCH PARAMETERS ---------------------------------------------------------------

pfw["endTime"] = stopTime*4.0
pfw["plotInterval"] = stopTime / 5
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
pfw["frictionCoefficient"]=0.25

pfw["planeStrain"] = 1 ##added for plane strain condition

 
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
 

def make_objects():

    impactor = geom.sphere('impactor',

        x0=[ pfw["xmin"] + sampleLength + 0.5*impacterDiameter,0.0,0.0],
        r=0.5*impacterDiameter,
        v=[-impactVelocity,0.,0.],
        mat=0,
        group=0,
        particleType=2)

    target = geom.cylinder('target',

        x1=[ pfw["xmin"],0.0,0.0],
        x2=[ pfw["xmin"]+sampleLength,0.0,0.0],
        r=0.5*sampleDiameter,
        v=[0.,0.,0.],
        mat=1,
        group=0,
        particleType=2)

 

    weibullTarget = geom.voronoiWeibullBoxWrapper('weibullWrapper',
        subObject=target,
        x0=[ pfw["xmin"], pfw["ymin"], pfw["zmin"] ],
        x1=[ pfw["xmax"], pfw["ymax"], pfw["zmax"] ],
        flawSize=5*DX,
        weibullVolume=weibullVolume,
        weibullModulus=weibullModulus,
        weibullSeed=1,
        vMin=(DX)**3.)

    objects=[impactor,weibullTarget]

    return objects

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = ["steel","target"]
pfw["materialPropertyString"]="""
<VonMisesJ

    name="steel"

    defaultDensity=""" + '"' + str(steelDensity) + '"' + """
    defaultYoungModulus=""" + '"' + str(steelYoungsModulus) + '"' + """
    defaultPoissonRatio=""" + '"' + str(steelPoissonsRatio) + '"' + """
    defaultYieldStrength=""" + '"' + str(steelYield) + '"' + """/>

 <CeramicDamage

    name="target"
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

pfw["boundaryConditionTypes"]=[4, 0, 4, 4, 0, 0]  
pfw["absorbingDampingFactor"]=[1000,0,1000,1000,0,0] 
pfw["plottableFields"]=["particleStrengthScale","target_damage", "gridVelocity","particleVelocity","particleWavespeed", "particleVolume"]