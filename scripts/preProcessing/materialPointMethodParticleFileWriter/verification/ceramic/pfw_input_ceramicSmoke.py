# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# Material model parameters.  Will be read from input file for plotting.

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8


# Prescribed deformation table: pure shear load and unload.
fTableString="""0	1	1	1
30	1	0.8	1
"""

domainScale=1.0
sampleHeight = 2.0*domainScale # mm
sampleWidth = 2.1*domainScale  # mm

# Domain ---------------------------------------------------------------------------------
cpp=24   # cells per partition in each direction
refine=6
xpar=refine  # grid partitions
ypar=refine
zpar=1

nI=xpar*cpp  # grid cells in the x-direction
nJ=ypar*cpp  # grid cells in the y-direction
nK=3  # grid cells in the z-direction
ppc=2   # particles per cell in each direction

# Define all the geometric objects -------------------------------------------------------
domainHeight = sampleHeight
domainWidth = sampleWidth  # This would be increased for unconfined compression.


xmin =-0.5*domainWidth # mm
xmax = 0.5*domainWidth # mm
ymin =-0.0*domainHeight  # mm
ymax = 1.0*domainHeight  # mm

planeStrain = 1

particleRefinement = [1,1]

dx = domainWidth/(nI-2)/ppc
dy = domainHeight/(nJ-2)/ppc
dz = 0.5*(dx+dy)

zmin =-0.5*dz
zmax = 0.5*dz

lx = xmax - xmin
ly = ymax - ymin
lz = zmax - zmin

# Batch parameters for GEOS runs.  --------------------------------------------------------
# An error will result if there are too many cores for
# a low resolution simulation.  If there is insufficient run-time to obtain a signal
# for a given run, that run will have its results ommited from the Hugoniot analysis.
mBatch=True
# read in the default bank:
import sys
import importlib
import getpass
username = getpass.getuser()
userDefsFile = str('userDefs_'+str(username))
userDefs = importlib.import_module(userDefsFile)
mBank = userDefs.defaultBank

mWallTime="00:01:00"
mCores=xpar*ypar*zpar
mNodes=int(np.ceil(float(mCores)/36.)) 
mSubmitJobs=False

# GEOS MPM input parameters ---------------------------------------------------------------
endTime="30.0"
useDamageAsSurfaceFlag=1  # this is needed to determine whether to flag particles on creation.
writePlot="1"
writeRestart="1"
plotInterval="0.25"
restartInterval="10.0"

# specify an array with all objects to be included, order matters. for overlapping objects, the first one listed will be assigned at each point.
# "fill" must be last on the list.

target = geom.cylinder('target',[xmin+0.5*lx,ymin+0.5*ly,zmin-lz],[xmin+0.5*lx,ymin+0.5*ly,zmax+lz],0.5*ly,vel=[0.0,0.0,0.0],mat=0,group=0)
objects=[target]

# List of materials:
materials = ["sand","plasticSand"]
materialPropertyString="""
<CeramicDamage
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
	compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
	maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
	crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
	/>
<PerfectlyPlastic
	name="plasticSand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	defaultYieldStress="""+'"'+str(tensileStrength)+'"'+"""
	/>
"""

mpmSolverParameterString="""
timeIntegrationOption="ExplicitDynamic"
cflFactor="0.7"    
initialDt="1e-16"

prescribedBcTable="0"
prescribedBoundaryFTable="1"
fTableInterpType="2"
solverProfiling="0"

planeStrain="""+"\""+str(planeStrain)+"\""+"""

neighborRadius="-1.01"
needsNeighborList="1"
useDamageAsSurfaceFlag="""+"\""+str(useDamageAsSurfaceFlag)+"\""+"""

cpdiDomainScaling="1"

surfaceDetection="0"
damageFieldPartitioning="1"
treatFullyDamagedAsSingleField="1"
separabilityMinDamage="0.5"
contactGapCorrection="1"
frictionCoefficient="0.0"

boundaryConditionTypes="{ 0, 0, 2, 2, 2, 2 }"    
"""

# New pfw dictionary interface -----------------------------------------------------------
pfw = {
    "runDebug": True,
    "xpar": xpar,
    "ypar": ypar,
    "zpar": zpar,
    "nI": nI,
    "nJ": nJ,
    "nK": nK,
    "ppc": ppc,
    "xmin": xmin,
    "xmax": xmax,
    "ymin": ymin,
    "ymax": ymax,
    "zmin": zmin,
    "zmax": zmax,
    "planeStrain": planeStrain,
    "particleRefinement": particleRefinement,
    "mBatch": mBatch,
    "mBank": mBank,
    "mWallTime": mWallTime,
    "mCores": mCores,
    "mNodes": mNodes,
    "mSubmitJobs": mSubmitJobs,
    "endTime": float(endTime),
    "plotInterval": float(plotInterval),
    "restartInterval": float(restartInterval),
    "objects": objects,
    "materials": materials,
    "materialPropertyString": materialPropertyString,
    "timeIntegrationOption": "ExplicitDynamic",
    "cflFactor": 0.7,
    "initialDt": 1e-16,
    "prescribedBcTable": 0,
    "prescribedBoundaryFTable": 1,
    "fTableInterpType": 2,
    "fTable": [[0.0, 1.0, 1.0, 1.0], [30.0, 1.0, 0.8, 1.0]],
    "solverProfiling": 0,
    "neighborRadius": -1.01,
    "needsNeighborList": 1,
    "useDamageAsSurfaceFlag": useDamageAsSurfaceFlag,
    "cpdiDomainScaling": 1,
    "damageFieldPartitioning": 1,
    "treatFullyDamagedAsSingleField": 1,
    "separabilityMinDamage": 0.5,
    "contactGapCorrection": 1,
    "frictionCoefficient": 0.0,
    "boundaryConditionTypes": [0, 0, 2, 2, 2, 2],
}
