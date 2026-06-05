# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 10.0

# DOMAIN ---------------------------------------------------------------------------------

domainWidth = 1.0
domainHeight = 1.0
domainLength = 1.0

pfw["xmin"] = -0.5*domainWidth  # mm
pfw["xmax"] =  0.5*domainWidth  # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

refine = 1 # partitions in each direction
cpp = 8   # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

dx = domainWidth / (pfw["nI"]-1)
dy = domainHeight / (pfw["nJ"]-1)
dz = domainLength / (pfw["nK"]-1)

# BATCH PARAMETERS --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1

pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["minParticleJacobian"]=1e-8
pfw["maxParticleJacobian"]=1e12
pfw["subdivideParticles"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
                             "SurfaceNormal",
                             "SurfacePosition",
							 "SurfaceTraction"]

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
                        "particleCohesiveZoneFlag",
                        "particleCohesiveForce",
                        "particleRVectors",
                        "particleReferenceRVectors",
                        "gridMass",
                        "gridDamage",
                        "gridPosition",
                        "gridVelocity",
                        "gridAcceleration",
                        "gridSurfaceNormal",
                        "gridSurfacePosition",
                        "gridContactForce",
                        "gridCohesiveNode"]

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "gas" ]
pfw["materialPropertyString"] ="""
<Gas
	name="gas"
	defaultDensity="0.1"
	referencePressure="0.101325"
    referenceTemperature="325"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-dx/2,-dy/2,-dz/2],[dx/2,dy/2,dz/2],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 1, 1, 1, 1, 1, 1 ]
