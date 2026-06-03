# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {}
pfw["runDebug"] = True
stopTime = 0.00001

# Domain ---------------------------------------------------------------------------------

sampleWidth = 1.0
sampleHeight = sampleWidth

domainWidth = sampleWidth
domainHeight = sampleHeight

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]

cpp = 25
refine = 6
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3                # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

domainLength = domainHeight/(pfw["nJ"]-2)

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

dy = domainHeight/(pfw["nJ"]-2)

# MATERIAL PROPERTIES --------------------------------------------------------------------

rho = 1.0
E = 10.0
nu = 0.22

pfw["materials"] = ["elastic"]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="elastic"
	defaultDensity=""" + '"' + str(rho) + '"' + """
	defaultYoungModulus=""" + '"' + str(E) + '"' + """
	defaultPoissonRatio=""" + '"' + str(nu) + '"' + """/>"""

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/2
pfw["restartInterval"]=stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["frictionCoefficient"]=0.0

pfw["plotGridFields"]=1
pfw["needsNodalNeighborList"]=1
pfw["computeParticleSurfaceNormalsAndPositions"]=1
pfw["normalAndPositionMethod"]="LogisticRegression"

pfw["updateMethod"]="PIC"

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
							 "Damage",
							 "StrengthScale",
                             "SurfaceFlag",
                             "RVector",
							#  "SurfaceNormal",
							#  "SurfacePosition"
							 ]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------
# radius = 0.2
# plate = geom.foam('plate',[-domainWidth/2, -domainHeight/2,pfw["zmin"]], [domainWidth/2, domainHeight/2, pfw["zmax"]], [[radius, 0.0, 0.0, 0.0]], vel=[0.0,0.0,0.0], mat=0, group=0)
# sphere = geom.sphere('sphere', [0.0,0.0,0.0], radius, vel=[0.0,0.0,0.0], mat=0, group=1)
# pfw["objects"]=[plate, sphere, fill]

#approximate a sphere with a very smooth polygon
dtheta = 2*np.pi/72
thetas = np.arange(0.0, 2*np.pi - dtheta, dtheta)
points = []
radius = 0.15
x = -0.3
y = 0.3
for theta in thetas:
	points.append([radius*np.cos(theta)+x, radius*np.sin(theta)+y])

print(points)

i = [[[0.0,0.0],
	  [0.3,0.0],
	  [0.3,0.3],
	  [0.0,0.3]],
	 [[-0.3,-0.3],
	  [-0.0,-0.3],
	  [-0.15,-0.0]]]
i.append(points)
inclusions = geom.polygonInclusions("inclusions",i,fillGroup=0,inclusionGroup=1)
pfw["objects"]=[inclusions]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 1, 1 ]
pfw["prescribedBoundaryFTable"]=0

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mBank"]="MAHEM"
pfw["mWallTime"]="00:15:00"
pfw["mSubmitJobs"]= True
