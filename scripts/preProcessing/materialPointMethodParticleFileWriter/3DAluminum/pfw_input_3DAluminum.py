# -*- coding: utf-8 -*-

import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = False

# Unit system :  mg, microsecond, mm,  stress-> GPA, velocity mm/us=km/s

impacterDiameter = 2 # aluminum sphere, mm
sampleDiameter = 152.4/2 # target diameter
sampleLength = 152.4 # target length (cylinder) , mm ##changed to 38.1 mm
impactVelocity = 2.077# mm/us. ##note the velocity changes here

stopTime = sampleLength/impactVelocity 

aluminumDensity = 2.784 #These are properties for aluminum
aluminumYield = .270
aluminumYoungsModulus = 70.
aluminumPoissonsRatio = 0.33

density = 2.1
bulk = 25./( 3.*( 1. - 2.*0.2) )  #poisson ratio changed from 0.3 to 0.2, E from 40 to 25
shear = 25./( 2.*( 1. + 0.2) ) 
#tensileStrength = 0.00625  #recipe R2
compressiveStrength = 0.025 # 25 MPa confirmed by Zhifei, Ryan
tensileStrength = compressiveStrength/8
maximumStrength = 5.0*compressiveStrength
crackSpeed = 0.8 # See Curbach paper
thirdInvariantDependence=1
refStrainRate=1e-10
rateSensitivity=0.2 # gives best results
m2=0.2

# weibullVolume = 100.0 # mm^3
# #weibullVolume = np.pi*sampleLength*(sampleDiameter**2)/4
# weibullModulus = 6.

weibullVolume = 27.5361
weibullModulus = 7.3692


# Domain ---------------------------------------------------------------------------------

domainX = sampleLength+impacterDiameter
domainY = sampleDiameter #changed to sampleDiameter from 1.1sampleDiameter #changed to 0.5 to get quarter symmetry
domainZ = sampleDiameter #same as above

cppx=16   # cells per partition in each direction #originaally it was 30,30,30
cppy=16
cppz=16

refine=6

pfw["xpar"]=2*refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] = 0.0 # mm
pfw["xmax"] = domainX # mm
pfw["ymin"] = 0 # mm
pfw["ymax"] = 0.5*domainY # mm
pfw["zmin"] = 0 # mm #changed from -0.5*domainZ
pfw["zmax"] = 0.5*domainZ # mm

DX = (pfw["xmax"] - pfw["xmin"])/pfw["xpar"]/cppx

# BATCH PARAMETERS --------------------------------------------------------

pfw["mBatch"]=True
pfw["mBank"]="imcomp" #"MAHEM"
pfw["mWallTime"]="4800000:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False #This prevents from automatic submission of Job
pfw["autoRestart"]=True #Trues

# END BATCH PARAMETERS ---------------------------------------------------------------

pfw["endTime"] = 50
pfw["plotInterval"] = 5
pfw["restartInterval"] = stopTime*10

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

pfw["materials"] = ["aluminum","target"]
pfw["materialPropertyString"]="""
<VonMisesJ

    name="aluminum"

    defaultDensity=""" + '"' + str(aluminumDensity) + '"' + """
    defaultYoungModulus=""" + '"' + str(aluminumYoungsModulus) + '"' + """
    defaultPoissonRatio=""" + '"' + str(aluminumPoissonsRatio) + '"' + """
    defaultYieldStrength=""" + '"' + str(aluminumYield) + '"' + """/>

 <CeramicDamage

    name="target"
    defaultDensity="""+'"'+str(density)+'"'+"""
    defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
    defaultShearModulus="""+'"'+str(shear)+'"'+"""
    tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
    compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
    maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
    crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
    thirdInvariantDependence=""" + '"' + str(thirdInvariantDependence) + '"' + """
    refStrainRate="""+'"'+str(refStrainRate)+'"'+"""
    rateSensitivity="""+'"'+str(rateSensitivity)+'"'+"""
    m2="""+'"'+str(m2)+'"'+"""/>

"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["boundaryConditionTypes"]=[1, 0, 1, 0, 1, 0]  
#pfw["absorbingDampingFactor"]=[500, 500, 1000, 1000, 1000, 1000]  
pfw["plottableFields"]=["particleStrengthScale","target_damage", "particleVelocity"]

#make the projectile look better
pfw["particleRefinement"]=[8,1]
