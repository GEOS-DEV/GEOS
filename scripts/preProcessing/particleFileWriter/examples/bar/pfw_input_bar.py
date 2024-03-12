# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
# from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 0.01

domainLength = 8.0 # m
domainWidth = 1.0  # m

# Domain ---------------------------------------------------------------------------------
pfw["xpar"]=1  # grid partitions
pfw["ypar"]=1
pfw["zpar"]=1

pfw["nI"]=10  	# grid cells in the x-direction
pfw["nJ"]=3  	# grid cells in the y-direction
pfw["nK"]=3  	# grid cells in the z-direction
pfw["ppc"]=2    # particles per cell in each direction

# Define all the geometric objects -------------------------------------------------------
pfw["xmin"] = 0.0 			
pfw["xmax"] = domainLength
pfw["ymin"] =-0.5*domainWidth
pfw["ymax"] = 0.5*domainWidth
pfw["zmin"] =-0.5*domainWidth
pfw["zmax"] = 0.5*domainWidth

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

bar = geom.box('bar',[0.0,-domainWidth/2,-domainWidth/2],[domainLength,domainWidth/2,domainWidth/2],[0.0,0.0,0.0],0,0,0)
pfw["objects"]=[bar]

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e9"
	defaultShearModulus="24.0e9"/>
"""
