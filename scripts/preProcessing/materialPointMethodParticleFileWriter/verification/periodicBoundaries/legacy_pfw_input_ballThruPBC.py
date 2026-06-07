# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 5.0

# DOMAIN ---------------------------------------------------------------------------------

refine = 1
cpp = 5
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["periodic"] = [True, True, True]

# For periodic boundaries don't add extra ghost cell
pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp		# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainWidth = 1.0 # m
domainHeight = domainWidth
domainLength = domainWidth

pfw["xmin"] =-0.5*domainWidth
pfw["xmax"] = 0.5*domainWidth
pfw["ymin"] =-0.5*domainHeight
pfw["ymax"] = 0.5*domainHeight
pfw["zmin"] =-0.5*domainLength
pfw["zmax"] = 0.5*domainLength

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]= stopTime/100
pfw["restartInterval"]=stopTime*10

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.5
pfw["initialDt"]=1e-16
pfw["solverProfiling"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 0, 0, 0, 0 ]

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e8"
	defaultShearModulus="24.0e8"/>
"""

# GEOMETRY OBJECTS  -------------------------------------------------------

ball = geom.sphere('disk1',[0.0, 0.0, 0.0], 0.8*domainWidth/2, vel=[1.0,1.0,1.0], mat=0, group=0)
pfw["objects"]=[ball]
