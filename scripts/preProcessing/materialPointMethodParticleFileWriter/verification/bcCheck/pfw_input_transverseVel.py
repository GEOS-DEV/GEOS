# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# crushing of 4 disks in 2D where each uses the graphite
# model with a different preferred direction

pfw = {} 
pfw["runDebug"] = True
stopTime = 20.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.46 # mass density mg/mm^3
K = 38.67 # bulk modulus (GPa)
G = 29.0 # shear modulus (GPa)

# Domain ---------------------------------------------------------------------------------

domainWidth = 1.0
domainHeight = 1.0

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [False, False, False]

refine = 5 # grid partitions
cpp = 8  # cells per partition in each direction

pfw["xpar"]=refine  
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3          # grid cells in the z-direction
pfw["ppc"]=2         # particles per cell in each direction

domainLength = domainHeight/(pfw["nJ"]-2)

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageWriteInterval"]=stopTime/2000
pfw["frictionCoefficient"]=0.25

# pfw["plotUnscaledParticles"]=1
# pfw["overlapCorrection"]=2
# pfw["overlapThreshold1"]=1.05
# pfw["overlapThreshold2"]=1.10

pfw["maxParticleVelocity"]=10.0               
pfw["minParticleJacobian"]=0.01                       
pfw["maxParticleJacobian"]=10.0  

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

pfw["useEvents"]=1
pfw["plotGridFields"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Define all the geometric objects -------------------------------------------------------

box_side = domainWidth*0.2

# -X boundary
box1=geom.box('box1',[-domainWidth/2, -box_side/2, pfw["zmin"]],[-domainWidth/2+box_side, box_side/2, pfw["zmax"]], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, False])

# +X boundary
box2=geom.box('box2',[domainWidth/2-box_side, -box_side/2, pfw["zmin"]],[domainWidth/2, box_side/2, pfw["zmax"]], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, False])

# -Y boundary
box3=geom.box('box3',[-box_side/2, -domainHeight/2, pfw["zmin"]],[box_side/2, -domainHeight/2 + box_side, pfw["zmax"]], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, False])

# +Y boundary
box4=geom.box('box4',[-box_side/2, domainHeight/2-box_side, pfw["zmin"]],[box_side/2, domainHeight/2, pfw["zmax"]], vel=[0.0, 0.0, 0.0], mat=0, group=0, dim=2, flaggedSurfaces=[False, False, False, False])

pfw["objects"]=[box1,box2, box3, box4]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 1, 1 ]

bc_vel = 0.01
pfw["enablePrescribedBoundaryTransverseVelocities"]=[1,1,1,1,0,0]
pfw["prescribedBoundaryTransverseVelocities"]=[[bc_vel, 0.0],
                                               [-bc_vel, 0.0],
                                               [0.0, bc_vel],
                                               [0.0, -bc_vel],
                                               [0.0, 0.0],
                                               [0.0, 0.0]]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0.0,          1.00,  1.00,  1.00],
               [stopTime,     1.10,  1.10,  1.00]]

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:15:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["materials"] = ["elasticIsotropic"]
pfw["materialPropertyString"] = f"""
<ElasticIsotropic
    name="elasticIsotropic"
    defaultDensity="{density}"
    defaultBulkModulus="{K}"
    defaultShearModulus="{G}"/>"""

# # --- PFW VERIFICATION FAST DEBUG OVERRIDES BEGIN ---
# # Debug-only runtime caps.  Keep this block below all source-file pfw assignments.
# def _vv_fast_int(_value, _default):
#     try:
#         return int(float(str(_value).strip().strip('"').strip("'")))
#     except Exception:
#         return int(_default)

# def _vv_fast_bool(_value):
#     if isinstance(_value, bool):
#         return _value
#     return str(_value).strip().strip('"').strip("'").lower() in ("1", "true", "yes", "on")

# try:
#     refine = 1
# except Exception:
#     pass

# # Fix common legacy typo before GEOS XML is written.
# if "planeStrain" in pfw and "planeStrain" not in pfw:
#     pfw["planeStrain"] = pfw.pop("planeStrain")

# _vv_fast_plane = _vv_fast_bool(pfw.get("planeStrain", False))
# # Treat thin 2D/plane-strain legacy cases as plane-like even when planeStrain was omitted.
# try:
#     if _vv_fast_int(pfw.get("zpar", 1), 1) == 1 and "nK" not in pfw:
#         _vv_fast_plane = True
# except Exception:
#     pass

# _vv_fast_cpp_cap = 24 if _vv_fast_plane else 8
# _vv_fast_max_partitions = 2
# pfw["mWallTime"] = "00:05:00"

# for _vv_key in ("xpar", "ypar", "zpar"):
#     pfw[_vv_key] = max(1, min(_vv_fast_int(pfw.get(_vv_key, 1), 1), _vv_fast_max_partitions))
# if _vv_fast_plane:
#     pfw["zpar"] = 1

# # Preserve already coarser grids, but cap high cells-per-partition values.
# def _vv_fast_cap_cells(_nkey, _pkey, _default_cells=1):
#     _p = max(1, _vv_fast_int(pfw.get(_pkey, 1), 1))
#     _n = _vv_fast_int(pfw.get(_nkey, 0), 0)
#     if _n <= 0:
#         return max(1, _p * min(_default_cells, _vv_fast_cpp_cap))
#     _cpp = max(1, (_n + _p - 1) // _p)
#     return max(1, _p * min(_cpp, _vv_fast_cpp_cap))

# pfw["nI"] = _vv_fast_cap_cells("nI", "xpar", _vv_fast_cpp_cap)
# pfw["nJ"] = _vv_fast_cap_cells("nJ", "ypar", _vv_fast_cpp_cap)
# if _vv_fast_plane:
#     if "nK" in pfw:
#         pfw["nK"] = max(1, min(_vv_fast_int(pfw.get("nK", 1), 1), 8))
# else:
#     pfw["nK"] = _vv_fast_cap_cells("nK", "zpar", 8)

# pfw["mCores"] = max(1, _vv_fast_int(pfw.get("xpar", 1), 1) * _vv_fast_int(pfw.get("ypar", 1), 1) * _vv_fast_int(pfw.get("zpar", 1), 1))
# # --- PFW VERIFICATION FAST DEBUG OVERRIDES END ---
