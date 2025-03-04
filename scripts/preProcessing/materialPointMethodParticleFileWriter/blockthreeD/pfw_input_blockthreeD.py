# -*- coding: utf-8 -*-

import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = False

# Unit system :  mg, microsecond, mm,  stress-> GPA, velocity mm/us=km/s

impacterDiameter = 3 # aluminum sphere, mm
sampleDiameter = 38.2/2 # target diameter
sampleLength = 38.2/8 # target length (cylinder) , mm
impactVelocity = 1.5 # mm/us

stopTime = sampleLength/impactVelocity 

aluminumDensity = 2.7
aluminumYield = .030
aluminumYoungsModulus = 70.
aluminumPoissonsRatio = 0.3

density = 2.75
bulk = 275./( 3.*( 1. - 2.*0.3) )
shear = 275./( 2.*( 1. + 0.3) )
tensileStrength = 0.30
compressiveStrength = 8.*tensileStrength
maximumStrength = 5.0*compressiveStrength
crackSpeed = 1.8
thirdInvariantDependence=1

weibullVolume = 1.0 # mm^3
#weibullVolume = np.pi*sampleLength*(sampleDiameter**2)/4
weibullModulus = 6.
# Domain ---------------------------------------------------------------------------------

domainX = sampleLength+impacterDiameter
domainY = sampleDiameter #changed to sampleDiameter from 1.1sampleDiameter
domainZ = sampleDiameter #same as above

cppx=30   # cells per partition in each direction #originaally it was 30,30,30
cppy=30
cppz=30

refine=1

pfw["xpar"]=2*refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=4   # particles per cell in each direction

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
pfw["mWallTime"]="48:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False #This prevents from automatic submission of Job
pfw["autoRestart"]=True #Trues

# END BATCH PARAMETERS ---------------------------------------------------------------

pfw["endTime"] = stopTime/2
pfw["plotInterval"] = stopTime/100
pfw["restartInterval"] = stopTime*5

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

    new_target = geom.box('new_target',

        x0=[ pfw["xmin"], pfw["ymin"], pfw["zmin"] ],
        x1=[ pfw["xmin"]+sampleLength, pfw["ymax"], pfw["zmax"] ],
        v=[0.,0.,0.],
        mat=1,
        group=0,
        particleType=2,
        dim=3,
        flaggedSurfaces=[False, False, False, False, False, False])

    weibullTarget = geom.voronoiWeibullBoxWrapper('weibullWrapper',
        subObject=new_target,
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

pfw["materials"] = ["aluminum","new_target"]
pfw["materialPropertyString"]="""
<VonMisesJ

    name="aluminum"

    defaultDensity=""" + '"' + str(aluminumDensity) + '"' + """
    defaultYoungModulus=""" + '"' + str(aluminumYoungsModulus) + '"' + """
    defaultPoissonRatio=""" + '"' + str(aluminumPoissonsRatio) + '"' + """
    defaultYieldStrength=""" + '"' + str(aluminumYield) + '"' + """/>

 <CeramicDamage

    name="new_target"
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

pfw["boundaryConditionTypes"]=[4, 0, 4, 4, 4, 4]  
pfw["absorbingDampingFactor"]=[500, 500, 1000, 1000, 1000, 1000]  
pfw["plottableFields"]=["particleStrengthScale","target_damage", "gridVelocity","particleVelocity","particleWavespeed","gridMass", "particleVolume", "particleMass"]
