# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree # nearest neighbor search with KDTree

# Initialize data structure:
pfw = {}
pfw["runDebug"] = False
pfw['particleFileFields'] = ["Velocity",  # +3
                                        "MaterialType", # +1
                                        "ContactGroup", # +1
                                        "SurfaceFlag", # +1
                                        "RVector",
                                        "MaterialDirection"]

# multiplier on partitions in each direction
refine=1                    

# PARAMETERS
# ------------------------------

ringInnerDiameter = (0.50)*25.4 # mm
ringOuterDiameter = (0.50)*25.4 # mm
# Sample dimensions used by the Voronoi material-direction cylinder below.
sampleDiameter = ringOuterDiameter
sampleLength = (0.30)*25.4 # mm

grainSize = 0.2*sampleDiameter
maxCompressiveStrain = 0.015

# Total simulation time, micro seconds (approximately 100 transit times)
stopTime = 100.*sampleLength/np.sqrt(15./1.92) 

# DOMAIN ---------------------------------------------------------------------------------
cpp=10  # 10   # cells per partition in each direction

domainZ = sampleDiameter
domainX = (0.60)*25.4 # mm  
domainY = (0.30)*25.4 # mm  

pfw["xmin"] = -0.5*domainX # mm
pfw["xmax"] = 0.5*domainX # mm
pfw["ymin"] = -0.5*domainY # mm
pfw["ymax"] = 0.5*domainY # mm
pfw["zmin"] = 0.0 # mm
pfw["zmax"] = domainZ

pfw["xpar"]=6*refine  # grid partitions
pfw["ypar"]=3*refine
pfw["zpar"]=5*refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True
pfw["mBank"]="mahem" #"MAHEM"
pfw["mWallTime"]="02:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20

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

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0
pfw["FSubcycles"]=10

pfw["plotUnscaledParticles"]=0

pfw["frictionCoefficient"]=0.05

pfw["contactGapCorrection"]="Implicit"
pfw["explicitSurfaceNormalInfluence"]= 1000.  # 1000
pfw["useSurfacePositionForContact"]= 1  # 1

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]="2"

pfw["plottableFields"]=["particleID",
                        "particleMass",
                        "particleVolume",
                        "particleDamage",
                        "particleInitialVolume",
                        "particleDensity",
                        "particleMaterialType",
                        "particleGroup",
                        "particleSurfaceFlag",
                        "particleStrengthScale",
                        "particleCenter",
                        "particleReferencePosition",
                        "particleVelocity",
                        "particleStress",
                        "particleSurfaceNormal",
                        "particleSurfacePosition",
                        "particleMaterialDirection",
                        "particleDeformationGradient",
                        "particlePlasticStrain",
                        "mass",
                        "damage",
                        "position",
                        "velocity",
                        "acceleration",
                        "surfaceNormal",
                        "surfacePosition",
                        "contactForce"]

# mpmSolverParameterString="""                               
# optimizeBinSort="0"                        
# binSizeMultiplier="4"                     
# neighborRadius="-1.0"                      
# overlapCorrection="0"
# overlapThreshold1="1.05"
# overlapThreshold2="1.10"                   
# """                 

# MATERIAL PROPERTIES ----------------------------------------------------------------

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8
thirdInvariantDependence=1
                                                  
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

# GEOMETRY OBJECTS -------------------------------------------------------------------------------------------
def make_objects():
    sample=geom.cylinder('cyl',
        x1=[0.0,-0.5*sampleLength,0.5*pfw["zmax"]],
        x2=[0.0,0.5*sampleLength,0.5*pfw["zmax"]],
        r=sampleDiameter/2,
        vel=[0.,0.,0.],
        mat=0,
        group=0)

    polycrystallineSample=geom.voronoiMatDirBoxWrapper('polycyl',
        subObject=sample,
        x0=[pfw["xmin"],pfw["ymin"],pfw["zmin"]],
        x1=[pfw["xmax"],pfw["ymax"],pfw["zmax"]],
        flawSize=grainSize,
        seed=1,
        dim=3
        )

    objects=[polycrystallineSample]
    return objects

# DEFORMATION ---------------------------------------------------------------------------------
pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]= [ 0, 0, 0, 0, 2, 2 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0, 1.00, 1.00, 1.00],
               [stopTime, 1.00, 1.00, 1.00 - maxCompressiveStrain]]

