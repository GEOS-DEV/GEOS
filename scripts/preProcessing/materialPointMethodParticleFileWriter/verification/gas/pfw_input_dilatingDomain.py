# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 20.0

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 1.0  # mm
sampleHeight = 1.0 # mm
sampleLength = 1.0 # mm

domainWidth = sampleWidth  # This would be increased for unconfined compression.
domainHeight = sampleHeight
domainLength = sampleLength

pfw["xmin"] = -0.5*domainWidth  # mm
pfw["xmax"] =  0.5*domainWidth  # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

refine=2  # partitions in each direction
cpp=4     # cells per partition in each direction

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
pfw["cflFactor"]=0.25  
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["maxParticleJacobian"]=1e12
pfw["subdivideParticles"]=1

# DEFORMATION -----------------------------------------------------------------------------

#Need to use surface tractions to dilate a material (keeps grid size constant to subdivide particles)

pfw["fTableInterpType"] = "Smoothstep"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,	         1,	    1,	    1],
               [stopTime,	 1.25,	1.25,	1.25]]

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2.7"
	defaultBulkModulus="70.0"
	defaultShearModulus="24.0"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]
