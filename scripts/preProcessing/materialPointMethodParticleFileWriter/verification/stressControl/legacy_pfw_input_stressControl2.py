# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True

stopTime = 10.

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

refine = 1 # paritions in each direction
cpp = 5    # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime
pfw["restartInterval"]=stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.05
pfw["initialDt"]=1e-16

pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/2000

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.93 # mass density (mg/mm^3)
E = 1.0       # Young's modulus (GPa)
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


pfw["fTableInterpType"] = "Smoothstep"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,	 1,	    1,	1],
               [.4*stopTime, 1.,	 0.8,	1],
               [.5*stopTime, 1.,	 0.8,	1],
               [stopTime,	 1,	    1.2,	1]
               ]

# This puts a cap on mDomain L, it can't be compressive if stress is below table value
# can't be tensile if above, but deformation will otherwise be fTable-controlled.
pfw["stressControl"]=[ 0, 2, 0]
pfw["stressTableInterpType"] = 'Linear'

maxCompressiveStress = 0.2
maxTensileStress = 0.1

pfw["stressTable"]=[[0.0,      	      0.0,      -maxCompressiveStress,  0.0],
					[.5*stopTime,    0.0,  -maxCompressiveStress,       0.0],
     				[.50000001*stopTime,    0.0,  maxTensileStress,     0.0],
					[stopTime,        0.0,       maxTensileStress,      0.0]]
