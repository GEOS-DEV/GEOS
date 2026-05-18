# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 1000.0
maxStretch = 3.0

# Domain ---------------------------------------------------------------------------------

refine=1

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

# We define this to be a single cell in the loading direction
# so we aren't affected by localization or shear banding.

pfw["nI"]=5     # grid cells in the x-direction
pfw["nJ"]=3     # grid cells in the y-direction
pfw["nK"]=5     # grid cells in the z-direction
pfw["ppcx"]=1   # particles per cell in each direction
pfw["ppcy"]=1   # particles per cell in each direction
pfw["ppcz"]=1   # particles per cell in each direction

# Define all the geometric objects -------------------------------------------------------
sampleWidth = 0.1  # mm
sampleHeight = 0.1 # mm
sampleLength = 0.1 # mm

domainWidth = 2.0*sampleWidth  # This would be increased for unconfined compression.
domainHeight = sampleHeight
domainLength = 2.0*sampleLength

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] =  -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:05:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/1000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/1000
pfw["useInteralForceAsFaceReaction"]=0

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1   # Would need to be zo for large stretch tension test if not single element or high ppc
pfw["damageFieldPartitioning"]=0

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["useEvents"]=1
pfw["frictionCoefficient"]=0.0

pfw["updateMethod"]="PIC"


# Deformation ---------------------------------------------------------------------------------
# Prescribed stretch in the loading direction, outflow boundary laterally to allow free 
# contraction or expansion.  There is a significant gap so the material shouldn't
# interact with lateral boundary.

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 0, 0 ]

pfw["fTableInterpType"]="Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,            1.0,   1.0,   1.0],
			   [stopTime,     1.0,   maxStretch, 1.0]]

# Define all the geometric objects -------------------------------------------------------

block = geom.box('block',[-sampleWidth/2, pfw["ymin"] ,-sampleLength/2],[sampleWidth/2,pfw["ymax"] ,sampleLength/2], vel=[0.0,0.0,0.0], mat=0, group=0)
pfw["objects"]=[block]

pfw["materials"] = [ "polymer" ]
pfw["materialPropertyString"]="""
<StrainHardeningPolymer
	name="polymer"
	defaultDensity="2.648"
	defaultBulkModulus="""+'"'+str(2.8)+'"'+"""
  bulkModulusA="0."
  bulkModulusB="4.56"
  bulkModulusT0="300."
  
	defaultShearModulus="""+'"'+str(0.2125)+'"'+"""                                         
	shearModulusA="0.0"
  shearModulusB="1.0"
  shearModulusT0="300."
  
  defaultYieldStrength="""+'"'+str(0.75*0.0215)+'"'+"""
  yieldStrengthA="0.0"
  yieldStrengthB="1.0"
  yieldStrengthT0="300.0"             
	
  strainHardeningSlope="""+'"'+str(0.0083)+'"'+"""                                         
  strainHardeningSlopeA="0.0" 
  strainHardeningSlopeB="1.0" 
  strainHardeningSlopeT0="300." 	
  
  shearSofteningMagnitude="""+'"'+str(0.00925)+'"'+"""                                                           
  shearSofteningMagnitudeA="0.0"                      
  shearSofteningMagnitudeB="1.0"                      
  shearSofteningMagnitudeT0="300."                      
  
	shearSofteningShapeParameter1="0.1"                      
	shearSofteningShapeParameter2="1.0"
 
  maximumStretch="1000."
  maximumStretchA="0.0"     
  maximumStretchB="1.0"     
  maximumStretchT0="2.95"    
  />
"""
