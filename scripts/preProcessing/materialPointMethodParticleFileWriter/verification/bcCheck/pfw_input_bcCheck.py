# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True

# DOMAIN ---------------------------------------------------------------------------------

pfw["planeSrain"]=1

refine=2
cpp=18

pfw["xpar"]=refine 
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=3  	            # grid cells in the z-direction
pfw["ppc"]=2                # particles per cell in each direction

domainLength = 1.0 # m
domainThickness = domainLength*(pfw["nK"]-2)/(pfw["nI"]-2)  # m, to get cubic cells

pfw["xmin"] =-0.5*domainLength	  # m
pfw["xmax"] = 0.5*domainLength	  # m
pfw["ymin"] =-0.5*domainLength	  # m
pfw["ymax"] = 0.5*domainLength	  # m
pfw["zmin"] =-0.5*domainThickness # m
pfw["zmax"] = 0.5*domainThickness # m


# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=0.01
pfw["plotInterval"]=0.0001
pfw["restartInterval"]=0.0025

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25   
pfw["initialDt"]=1e-16   

pfw["solverProfiling"]=0        

pfw["contactGapCorrection"]=1
pfw["frictionCoefficient"]=0.25

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 1, 1, 1, 1, 1, 1 ]    

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e9"
	defaultShearModulus="24.0e9"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

rad = 0.147
disk1 = geom.cylinder( 'disk1', [ pfw["xmax"] - 1.5*rad, pfw["ymax"] - 1.5*rad, pfw["zmin"] ], [ pfw["xmax"] - 1.5*rad, pfw["ymax"] - 1.5*rad, pfw["zmax"] ], rad, vel=[ -120, -120.0, 0.0 ], mat=0, group=0)
disk2 = geom.cylinder( 'disk1', [ pfw["xmin"] + 1.5*rad, pfw["ymax"] - 1.5*rad, pfw["zmin"] ], [ pfw["xmin"] + 1.5*rad, pfw["ymax"] - 1.5*rad, pfw["zmax"] ], rad, vel=[ 120, -120.0, 0.0 ], mat=0, group=0)
pfw["objects"]=[disk1, disk2]
