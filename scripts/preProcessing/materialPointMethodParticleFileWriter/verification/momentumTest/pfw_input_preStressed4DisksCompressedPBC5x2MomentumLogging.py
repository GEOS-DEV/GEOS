# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree


# 2-D plane strain borhole collapse simulation,
# inner diameter will have pressure BC
# outer boundary will have rigid stress-controlled motion.

pfw = {} 
pfw["runDebug"] = True

stopTime = 500.0
initialVelocity = 0.01
initialPressure = 1.0

# units: mm, mg, us (stres is GPa)

# Domain ---------------------------------------------------------------------------------
pfw["planeStrain"] = 1

domainX = 5.0
domainY = 1.0

# Typically we run 12-16 in 3D, 30-60 in 2D (30,000 particles per partition)
cppx=24   # cells per partition in each direction
cppy=12
cppz=3   # plane strain needs 3, 

pfw["xpar"]=5  # grid partitions
pfw["ypar"]=2
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] = -2.5 # mm
pfw["xmax"] =  2.5 # mm
pfw["ymin"] = -0.5*domainY # mm
pfw["ymax"] = 0.5*domainY # mm

# for plane strain we compute depth to equal grid cell size.
DZ = ( pfw["xmax"] - pfw["xmin"] ) / ( pfw["nI"] - 2 )
domainZ = DZ
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True    # run in background vs. in salloc interactive
pfw["mBank"]= "mahem"
pfw["mWallTime"]="01:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# END BATCH PARAMETERS ---------------------------------------------------------------
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 60
pfw["restartInterval"] = stopTime*5.0

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------
pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1

# pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0
pfw["FSubcycles"]=10

# pfw["plotUnscaledParticles"]=0

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]="2"

pfw["useEvents"]=1

pfw["particleFileFields"] = ["Velocity",
                            "MaterialType",
                            "ContactGroup",
                            "SurfaceFlag",
                            "RVector"]

pfw["logMomentum"]=1

# END GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------

# GEOMETRY OBJECTS -------------------------------------------------------------------------------------------
d1 = geom.cylinder('d1',
    x1 = [-1.0, 0.0, pfw["zmin"] ],
    x2 = [-1.0, 0.0, pfw["zmax"] ],
    r=0.4,
    vel=np.array( [ -initialVelocity, 0.5*initialVelocity, 0. ] ),
    mat = 0,
    group = 0 )
d2 = geom.cylinder('d1',
    x1 = [0.0, 0.0, pfw["zmin"] ],
    x2 = [0.0, 0.0, pfw["zmax"] ],
    r=0.4,
    vel=np.array( [ +initialVelocity, -0.25*0.5*initialVelocity, 0. ] ),
    mat = 0,
    group = 0 )
d3 = geom.cylinder('d1',
    x1 = [ 1.0, 0.0, pfw["zmin"] ],
    x2 = [ 1.0, 0.0, pfw["zmax"] ],
    r=0.4,
    vel=np.array( [ -initialVelocity, 0., 0. ] ),
    mat = 0,
    group = 0 )
d4 = geom.cylinder('d1',
    x1 = [ 2.0, 0.0, pfw["zmin"] ],
    x2 = [ 2.0, 0.0, pfw["zmax"] ],
    r=0.4,
    vel=np.array( [ initialVelocity, -0.25*initialVelocity, 0. ] ),
    mat = 0,
    group = 0 )



pfw["objects"]=[d1,d2,d3,d4]

# MATERIAL PROPERTIES ----------------------------------------------------------------
# Import material properties from validation materials file --------------------------
# [pfw_dependency] /pfw_materials.py
import importlib
matFile = importlib.import_module('pfw_materials')
aluminum = matFile.aluminum

# Assign imported model to pfw dictionary.
pfw["materials"] = [aluminum["name"]]
pfw["materialPropertyString"]=aluminum["materialString"]
# ------------------------------------------------------------------------------------

# DEFORMATION ---------------------------------------------------------------------------------
pfw["periodic"] = [True, False, False]
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]
pfw["fTableInterpType"]='Smoothstep'
pfw["prescribedFTable"]=0
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0, 1.00, 1.00, 1.00],
[stopTime,         1.00, 0.50, 1.00]]

# MPM EVENTS -------------------------------------------------------------------------------
# attempting stress boundary condition with confiningPressure Event.
# the box should be interior to the domain in the dimensions where
# the pressure is to be applied (XY, in this case) but the box
# shouldn't overlap with the borholePressure region
pfw["mpmEventsString"]="""
<InitializeStress 
    startTime="0.0"
    targetRegion="all"
    pressure=""" + '"' + str(initialPressure) + '"' + """
/>
"""
