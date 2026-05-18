# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree


# 2-D plane strain borhole collapse simulation,
# inner diameter will have pressure BC
# outer boundary will have rigid stress-controlled motion.

pfw = {} 
pfw["runDebug"] = True

stopTime = 1.0
initialVelocity = 0.001
initialPressure = 1.0

# units: mm, mg, us (stres is GPa)

# Domain ---------------------------------------------------------------------------------
pfw["planeStrain"] = 1

domainX = 5.0
domainY = 1.0

# Typically we run 12-16 in 3D, 30-60 in 2D (30,000 particles per partition)
cppx=5   # cells per partition in each direction
cppy=5
cppz=3   # plane strain needs 3, 

pfw["xpar"]=5  # grid partitions
pfw["ypar"]=3
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] = -2.5 # mm
pfw["xmax"] =  2.5 # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = 3.0 # mm

# for plane strain we compute depth to equal grid cell size.
DZ = ( pfw["xmax"] - pfw["xmin"] ) / ( pfw["nI"] - 2 )
domainZ = DZ
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True    # run in background vs. in salloc interactive
pfw["mWallTime"]="00:10:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# END BATCH PARAMETERS ---------------------------------------------------------------
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 10
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
sample = geom.twoMaterialBox('block',
    x0=[pfw["xmin"],pfw["ymin"],pfw["zmin"]],
    x1=[pfw["xmax"],pfw["ymax"],pfw["zmax"]],
    vel=np.array([initialVelocity,0.,0.]),
    mat1=0,
    mat2=1,
    group=0)

pfw["objects"]=[sample]

# MATERIAL PROPERTIES ----------------------------------------------------------------
# Import material properties from validation materials file --------------------------
# [pfw_dependency] /pfw_materials.py
import importlib
matFile = importlib.import_module('pfw_materials')
aluminum = matFile.aluminum
steel = matFile.steel

# Assign imported model to pfw dictionary.
pfw["materials"] = [aluminum["name"], steel["name"]]
pfw["materialPropertyString"]=aluminum["materialString"]+steel["materialString"]
# ------------------------------------------------------------------------------------

# DEFORMATION ---------------------------------------------------------------------------------
# This could be [0,0,0,0,1,1], but then the material might drift out of the domain and get 
# deleted.  We'll see if the confiningPressure event works with material against the 
# bounday (it should).  contact BC might be better here.
pfw["boundaryConditionTypes"]=[0,0,1,1,1,1]   
pfw["periodic"]=[True,False,False]

# MPM EVENTS -------------------------------------------------------------------------------
# attempting stress boundary condition with confiningPressure Event.
# the box should be interior to the domain in the dimensions where
# the pressure is to be applied (XY, in this case) but the box
# shouldn't overlap with the borholePressure region
pfw["mpmEventsString"]="""
<InitializeStress 
    startTime="0.0"
    targetRegion="all"
    pressure=
""" + '"' + str(initialPressure) + '"' + """
/>
"""
