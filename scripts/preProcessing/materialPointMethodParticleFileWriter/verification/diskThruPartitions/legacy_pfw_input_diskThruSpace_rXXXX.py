# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# START!

# Domain ---------------------------------------------------------------------------------
refine=XXXX
xpar=1  # grid partitions
ypar=1
zpar=1

cpp=2*refine
nI=2+xpar*cpp  	# grid cells in the x-direction
nJ=2+ypar*cpp  	# grid cells in the y-direction
nK=2+zpar*cpp		# grid cells in the z-direction
ppc=2   		# particles per cell in each direction

domainLength = 1.0 # m

# Define all the geometric objects -------------------------------------------------------
xmin =-0.5*domainLength	# m
xmax = 0.5*domainLength	# m
ymin =-0.5*domainLength	# m
ymax = 0.5*domainLength	# m
zmin =-0.5*domainLength 	# m
zmax = 0.5*domainLength 	# m

lx = xmax - xmin
ly = ymax - ymin
lz = zmax - zmin


# Batch parameters for GEOS runs.  --------------------------------------------------------
# The suite runner fills the scheduler bank when staging jobs; keep source inputs
# independent of user-specific userDefs modules so preflight checks are portable.
mBank = ""

minutes = int(2+(refine*refine*refine)*0.0005)
hours=int((minutes-minutes%60)/60)
minutes -= 60*hours
mWallTime=str(hours)+":"+("00"+str(minutes))[-2:]+":00"

mCores=xpar*ypar*zpar
mSubmitJobs=True

# GEOSX MPM input parameters ---------------------------------------------------------------
endTime=str(2.0/refine)		# seconds
writePlot="0"				# this does nothing for now
writeRestart="0"			# this does nothing for now
plotInterval="10"
restartInterval="1e9"

# specify an array with all objects to be included, order matters. for overlapping objects, the first one listed will be assigned at each point.
# "fill" must be last on the list.
ball = geom.sphere('disk1',[0.,0.,0.],0.4,vel=[-1.0,-1.0,1.0],mat=0,group=0)
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
	defaultDensity="1000"
	defaultBulkModulus="1000"
	defaultShearModulus="1000"/>
"""

mpmSolverParameterString="""
timeIntegrationOption="ExplicitDynamic"
cflFactor="0.5"
initialDt="1e-16"
solverProfiling="2"
boundaryConditionTypes="{ 0, 0, 0, 0, 0, 0 }"
"""
