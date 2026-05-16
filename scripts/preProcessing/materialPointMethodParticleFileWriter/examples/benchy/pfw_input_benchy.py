# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# START!

# Domain ---------------------------------------------------------------------------------
xpar=10  # grid partitions
ypar=5
zpar=8

cpp=15
nI=xpar*cpp  	# grid cells in the x-direction
nJ=ypar*cpp  	# grid cells in the y-direction
nK=zpar*cpp 	# grid cells in the z-direction
ppc=2   		# particles per cell in each direction

domainLength = 60.001 # m
domainWidth = 31.004 # m
domainHeight = 48 # m

# Currently there is a bug that will cause a crash if the problem is initialized with
# an empty patch.  This option will create one particle of material 0 in each grid cell.
# This currently does nothing!
mpiParticles = False

# Define all the geometric objects -------------------------------------------------------
xmin =-29.176	# m
xmax = 30.825	# m
ymin =-15.502	# m
ymax = 15.502	# m
zmin = 0.0 	# m
zmax = 48 	# m

# Batch parameters for GEOS runs.  --------------------------------------------------------
# An error will result if there are too many cores for
# a low resolution simulation.  If there is insufficient run-time to obtain a signal
# for a given run, that run will have its results ommited from the Hugoniot analysis.
mBatch=False

# read in the default bank:
import sys
import importlib
import getpass
username = getpass.getuser()
userDefsFile = str('userDefs_'+str(username))
userDefs = importlib.import_module(userDefsFile)
mBank = userDefs.defaultBank

mWallTime="00:05:00"
mCores=xpar*ypar*zpar
mNodes=int(np.ceil(float(mCores)/36.)) 
mSubmitJobs=False

# GEOSX MPM input parameters ---------------------------------------------------------------
endTime="0.01"				# seconds
writePlot="1"				# this does nothing for now
writeRestart="1"			# this does nothing for now
plotInterval="0.0001"
restartInterval="0.0025"

# specify an array with all objects to be included, order matters. for overlapping objects, the first one listed will be assigned at each point.
# "fill" must be last on the list.

ball = geom.sphere('ball',[-20,0.0,27],5,[100.0,0.0,0.0],0,1,0)

objects=[ball]

# Material Properties
# Notes:
# 1. It DOES NOT matter if the 'materials' array matches the order of 'materialPropertyString'.
# 2. It's okay if there are more materials specified than used.
# 3. The 'materials' array is a map from material name to material ID in the particle file.
#    i.e. if the first material in 'materials' is 'sand', then 'sand' is material ID 0.
#    This matters when using integers to flag objects with materials.
materials = [ "aluminum" ]
materialPropertyString="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e8"
	defaultShearModulus="24.0e8"/>
"""

mpmSolverParameterString="""
timeIntegrationOption="ExplicitDynamic"
cflFactor="0.25"    
initialDt="1e-16"

prescribedBcTable="0"
prescribedBoundaryFTable="0"
fTableInterpType="0"    

solverProfiling="0"         

contactGapCorrection="1"
frictionCoefficient="0.25"

boundaryConditionTypes="{ 1, 1, 1, 1, 1, 1 }"    
"""
