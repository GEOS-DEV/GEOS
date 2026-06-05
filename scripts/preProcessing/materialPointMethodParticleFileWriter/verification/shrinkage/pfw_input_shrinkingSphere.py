# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 2.0

# Domain ---------------------------------------------------------------------------------

refine = 1
cpp = 8  # cells per partition in each direction
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

# Define all the geometric objects -------------------------------------------------------
sampleWidth = 0.1  # mm
sampleHeight = 0.1 # mm
sampleLength = 0.1 # mm

substrateThickness = 0.1 # mm

domainWidth = 2.0*sampleWidth  # This would be increased for unconfined compression.
domainHeight = 2.0*sampleHeight
domainLength = 2.0*sampleLength

pfw["xmin"] =  -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] =  -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] =  -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

dx = domainWidth / ( pfw["nI"] - 2)
dy = domainHeight / ( pfw["nJ"] - 2)
dz = domainLength / ( pfw["nK"] - 2)

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0
pfw["frictionCoefficient"]=0.0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

pfw["useEvents"]=1
# TODO: Add passthrough to solver.
pfw["initialTemperature"] = 300.0
# pfw["shrinkageTable"]=[[293.0,  0.000],
#                        [373.0,  0.000],
# 					   [473.0,  0.000],
# 					   [573.0,  0.213],
# 					   [673.0,  0.307],
# 					   [773.0,  0.359],
# 					   [873.0,  0.407],
# 					   [973.0,  0.420],
# 					   [1073.0, 0.420],
# 					   [1173.0, 0.420],
# 					   [1273.0, 0.420],
# 					   [1373.0, 0.420],
# 					   [1473.0, 0.420],
# 					   [1573.0, 0.420]]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 1, 1, 1, 1, 1, 1 ]

# Define all the geometric objects -------------------------------------------------------

pillar = geom.sphere('pillar', [0.0, 0.0, 0.0], sampleWidth/2, vel=[0.0, 0.0, 0.0], mat=0, group=0)
pillarWFlag = geom.shrinkageFlagWrapper('pillarWFlag', pillar, 1)

pfw["objects"]=[ pillarWFlag ]

pfw["materials"] = [ "FusedSilica" ]
pfw["materialPropertyString"] ="""
<ElasticIsotropic
	name="FusedSilica"
	defaultDensity="2.12"
	defaultBulkModulus="35"
	defaultShearModulus="31"
	defaultThermalExpansionCoefficient="0.5e-6"/>
"""

pfw["mpmEventsString"]="""
<TemperatureProfile
	time=
""" + '"' + str(0.0) + '"' + """
	interval=""" + '"' + str(stopTime) + '"' + """
	temperatureTable=""" + '"' + """{{0.0, """ + str(pfw["initialTemperature"]) + '}, {'  + str(stopTime) + """, 1273.0}}""" + '"' + """/>
</MPMEvents>+
"""
