# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree # nearest neighbor search with KDTree
import math

pfw = {}
pfw["runDebug"] = True

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 1.0
sampleHeight = 1.0

domainWidth = sampleWidth*1.25
domainHeight = sampleHeight

pfw["xmin"] = -0.5*domainWidth  # mm
pfw["xmax"] = 0.5*domainWidth   # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainHeight  # mm

pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]


refine = 4
cpp = 10

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3                # grid cells in the z-direction
pfw["ppcx"]=2               # particles per cell in x-direction
pfw["ppcy"]=4               # particles per cell in y-direction

dx = domainWidth / (pfw["nI"])
dy = domainHeight / (pfw["nJ"]-2)

domainLength = 5*dx

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBank"]="MAHEM"
pfw["mBatch"]=True
pfw["mWallTime"]="06:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/56.))
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# MATERIAL PROPERTIES ----------------------------------------------------------------

density = 1.93 # mass density mg/mm^3
E = 45.0     # in-plane stiffness (GPa)
nu = 0.27    # in plane poisson's ratio

#~~~~~~~~~~~~~~~~~~~~~~~ Temperature dependent strain hardening polymer parameters ~~~~~~~~~~~~~~~~~
# K - bulk modulus
bulkModulusDefault=0.49
bulkModulusA= 1.7
bulkModulusB= -0.1
bulkModulusT0= 309.

# G - shear modulus
shearModulusDefault=0.05
shearModulusA= 60.
shearModulusB= 0.17
shearModulusT0= 253.

# sigma_y0 - yield strength
yieldStrengthDefault=0.0008
yieldStrengthA= 60.
yieldStrengthB= 0.07
yieldStrengthT0= 270.

maximumStretch=0.5
maximumStretchA= 1.5
maximumStretchB= -0.1
maximumStretchT0= 283.15

# Gr - strainhardening slope
strainHardeningSlopeDefault=0.015
strainHardeningSlopeA= -1.
strainHardeningSlopeB= -0.08
strainHardeningSlopeT0= 280.

# r0 - shearsofteningmagnitude
shearSofteningMagnitudeDefault=0.1
shearSofteningMagnitudeA= -1.0
shearSofteningMagnitudeB= -0.15
shearSofteningMagnitudeT0= 260.

#r1
shearSofteningShapeParameter1 = 0.1    

#r2
shearSofteningShapeParameter2 = 1.0

pfw["materials"] = ["isotropic", "polymer"]
pfw["materialPropertyString"]="""
<ElasticIsotropic
    name="isotropic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultYoungModulus=""" + '"' + str(E) + '"' + """
    defaultPoissonRatio=""" + '"' + str(nu) + '"' + """/>    
<StrainHardeningPolymer
    name="polymer"
    defaultDensity="2.648"
    defaultBulkModulus=""" + '"' + str(bulkModulusDefault)  + '"' + """
    bulkModulusA=""" + '"' + str(bulkModulusA)  + '"' + """
    bulkModulusB=""" + '"' + str(bulkModulusB)  + '"' + """
    bulkModulusT0=""" + '"' + str(bulkModulusT0)  + '"' + """
    defaultShearModulus=""" + '"' + str(shearModulusDefault)  + '"' + """
    shearModulusA=""" + '"' + str(shearModulusA)  + '"' + """
    shearModulusB=""" + '"' + str(shearModulusB)  + '"' + """
    shearModulusT0=""" + '"' + str(shearModulusT0)  + '"' + """   
    defaultYieldStrength=""" + '"' + str(yieldStrengthDefault)  + '"' + """
    yieldStrengthA=""" + '"' + str(yieldStrengthA)  + '"' + """
    yieldStrengthB=""" + '"' + str(yieldStrengthB)  + '"' + """
    yieldStrengthT0=""" + '"' + str(yieldStrengthT0)  + '"' + """
    maximumStretch=""" + '"' + str(maximumStretch)  + '"' + """ 
    maximumStretchA=""" + '"' + str(maximumStretchA)  + '"' + """
    maximumStretchB=""" + '"' + str(maximumStretchB)  + '"' + """
    maximumStretchT0=""" + '"' + str(maximumStretchT0)  + '"' + """                        
    strainHardeningSlope=""" + '"' + str(strainHardeningSlopeDefault)  + '"' + """
    strainHardeningSlopeA=""" + '"' + str(strainHardeningSlopeA)  + '"' + """
    strainHardeningSlopeB=""" + '"' + str(strainHardeningSlopeB)  + '"' + """
    strainHardeningSlopeT0=""" + '"' + str(strainHardeningSlopeT0)  + '"' + """
    shearSofteningMagnitude=""" + '"' + str(shearSofteningMagnitudeDefault)  + '"' + """
    shearSofteningMagnitudeA=""" + '"' + str(shearSofteningMagnitudeA)  + '"' + """                    
    shearSofteningMagnitudeB=""" + '"' + str(shearSofteningMagnitudeB)  + '"' + """                  
    shearSofteningMagnitudeT0=""" + '"' + str(shearSofteningMagnitudeT0)  + '"' + """                   
    shearSofteningShapeParameter1=""" + '"' + str(shearSofteningShapeParameter1)  + '"' + """                      
    shearSofteningShapeParameter2=""" + '"' + str(shearSofteningShapeParameter2)  + '"' + """/>  
"""

# EVENTS ----------------------------------------------------------------

# loadTime = 4.0
# stopTime = 4.4

# pfw["useEvents"]=1
# pfw["mpmEventsString"]="""
# <MPMEvents>
#     <Anneal 
#         time=""" + '"' + str(loadTime) + '"' + """
#         interval=""" + '"' + str(0.1) + '"' + """
#         targetRegion="all"/>
#     <PolymerHeal
#         time=""" + '"' + str(loadTime+0.1) + '"' + """
#         interval=""" + '"' + str(0.1) + '"' + """
#         targetRegion="all"/>
# </MPMEvents>
# """

loadTime = 5.0
unloadTime = loadTime + 5.0

annealTimeInterval = 0.1

crystalHealTime = unloadTime + annealTimeInterval
crystalHealTimeInterval = 0.2

testTime = crystalHealTime + crystalHealTimeInterval

stopTime = testTime + 5.0

pfw["useEvents"]=1
pfw["mpmEventsString"]="""
<MPMEvents>
    <Anneal 
        time=""" + '"' + str(unloadTime) + '"' + """
        interval=""" + '"' + str(annealTimeInterval) + '"' + """
        targetRegion="all"/>
    <PolymerHeal
        time=""" + '"' + str(crystalHealTime) + '"' + """
        interval=""" + '"' + str(crystalHealTimeInterval) + '"' + """
        targetRegion="all"/>
</MPMEvents>
"""

# DEFORMATION ---------------------------------------------------------------------------------


pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1]

pfw["prescribedFTable"]=0
pfw["prescribedBoundaryFTable"]=1
pfw["fTableInterpType"]="Cosine"
pfw["fTable"]=[[0.0,         1.00,  1.00,      1.00],
               [loadTime,    1.00,  1.00+0.5,  1.00],
               [unloadTime,  1.00,  1.00-0.1,  1.00],
               [testTime,    1.00,  1.00-0.1,  1.00],
               [stopTime,    1.00,  1.00+0.5,  1.00]]

# pfw["fTable"]=[[0.0,         1.00,  1.00,      1.00],
#                [loadTime,    1.00,  1.00+0.2,  1.00],
#                [stopTime,    1.00,  1.00+0.2,  1.00]]

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20
pfw["reactionWriteInterval"] = stopTime/2000
pfw["boxAverageWriteInterval"]= stopTime/2000

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1
pfw["separabilityMinDamage"]=0.5               
pfw["treatFullyDamagedAsSingleField"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["frictionCoefficient"]=0.25

pfw["plotUnscaledParticles"]=1
pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0
pfw["FSubcycles"]=10

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

pfw["contactGapCorrection"]="Implicit"
pfw["useSurfacePositionForContact"]= 0
pfw["explicitSurfaceNormalInfluence"]= 0

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "Damage",
                             "StrengthScale",
                             "Temperature",
                             "RVector",
                             "MaterialDirection",
                             "SurfaceNormal",
                             "SurfacePosition"]

pfw["plottableFields"]=["particleID",
                        "particleMass",
                        "particleVolume",
                        "particleDamage",
                        "particleDensity",
                        "particleMaterialType",
                        "particleGroup",
                        "particleSurfaceFlag",
                        "particleStrengthScale",
                        "particleTemperature",
                        "particleCenter",
                        "particleReferencePosition",
                        "particleVelocity",
                        "particleStress",
                        "particleSurfaceNormal",
                        "particleSurfacePosition",
                        "particleMaterialDirection",
                        "particleDeformationGradient",
                        "particlePlasticStrain"]

# GEOMETRY OBJECTS -------------------------------------------------------
binder_thickness = 0.2

def make_objects():
    top_platen = geom.box("top_platen", [-sampleWidth/2, domainHeight - (domainHeight-binder_thickness)/2, 0.0], [ sampleWidth/2, domainHeight, 0.0], mat=0, group=0,dim=2, flaggedSurfaces=[False, False, False, False])
    binder = geom.box("binder", [-sampleWidth/2, domainHeight - (domainHeight-binder_thickness)/2, 0.0], [ sampleWidth/2, (domainHeight-binder_thickness)/2, 0.0], mat=1, group=0,dim=2, flaggedSurfaces=[False, False, False, False])
    bottom_platen = geom.box("bot_platen", [-sampleWidth/2, 0.0, 0.0], [ sampleWidth/2, (domainHeight-binder_thickness)/2, 0.0], mat=0, group=0,dim=2, flaggedSurfaces=[False, False, False, False])
    
    return [top_platen, binder, bottom_platen]

