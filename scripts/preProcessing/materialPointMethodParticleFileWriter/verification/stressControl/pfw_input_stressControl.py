# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True

preloadTime = 10.0
stressRate = 0.1 / 20
pressure = 0.1
maxStress = -0.2

testTime = preloadTime + abs( maxStress ) / stressRate
stopTime = max( abs(-pressure + maxStress), abs(pressure) ) / stressRate + testTime

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 0.1  # mm
sampleHeight = 0.1 # mm
sampleLength = 0.1 # mm

domainWidth = sampleWidth
domainHeight = sampleHeight
domainLength = sampleLength

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] =  0.0             # mm
pfw["ymax"] = domainHeight     # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

refine=1 # paritions in each direction
cpp=5    # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime           
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.05 
pfw["initialDt"]=1e-16

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/2000
pfw["frictionCoefficient"]=0.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.93 # mass density (mg/mm^3)
E = 45.0       # Young's modulus (GPa)
nu = 0.27      # Poisson's ratio (-)

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"]="""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultYoungModulus=""" + '"' + str(E) + '"' + """
    defaultPoissonRatio=""" + '"' + str(nu) + '"' + """/>"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]

pfw["fTableInterpType"]='Cosine'
pfw["prescribedFTable"]=0
pfw["prescribedBoundaryFTable"]=0

pfw["stressControl"]=[ 1, 1, 1]
pfw["stressTableInterpType"] = 'Cosine'
pfw["stressControlKp"] = 1.0
pfw["stressControlKi"] = 0.0
pfw["stressControlKd"] = 0.0
pfw["stressTable"]=[[0.0,      	      0.0,       0.0,                   0.0],
					[preloadTime,     -pressure, -pressure,             -pressure],
					[testTime,        -pressure, -pressure + maxStress, -pressure],
					[stopTime,        0.0,       0.0,                   0.0]]
