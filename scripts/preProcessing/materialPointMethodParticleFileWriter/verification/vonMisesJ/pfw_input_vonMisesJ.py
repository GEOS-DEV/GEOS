# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 10.0

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 1.0
sampleHeight = 1.0
sampleLength = 1.0

domainWidth = 1.25*sampleWidth
domainHeight = sampleHeight
domainLength = 1.25*sampleLength

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

pfw["periodic"] = [False, False, False]

refine=1 # partitions in each direction
cpp = 10 # cells per partition in each direction

pfw["xpar"]=refine  # grid partitions
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
pfw["plotInterval"]=stopTime/10
pfw["restartInterval"]=stopTime*10 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.00 # mass density mg/mm^3
E = 5.0   # Young's modulus (GPa)
nu = 0.22   # Poisson's ratio (-)
strength = 0.8 # yield strength (GPa)

pfw["materials"] = ["elasticPlastic"]
pfw["materialPropertyString"]="""
<VonMisesJ
    name="elasticPlastic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultYoungModulus=""" + '"' + str(E) + '"' + """
    defaultPoissonRatio=""" + '"' + str(nu) + '"' + """
    defaultYieldStrength=""" + '"' + str(strength) + '"' + """/>"""

# GEOMETRY OBJECTS -------------------------------------------------------

box=geom.box('box',[-sampleWidth/2, -sampleHeight/2, -sampleLength/2],[sampleWidth/2, sampleHeight/2, sampleLength/2], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=3)
pfw["objects"]=[box]

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.00,  1.20,  1.00]]
