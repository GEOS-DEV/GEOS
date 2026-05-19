# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree


# 2-D plane strain borhole collapse simulation,
# inner diameter will have pressure BC
# outer boundary will have rigid stress-controlled motion.

pfw = {} 
pfw["runDebug"] = True

stopTime = 3000.0
initialVelocity = 0.01
initialPressure = 1.0

# units: mm, mg, us (stres is GPa)

# Domain ---------------------------------------------------------------------------------
pfw["planeStrain"] = 1

domainX = 3.0
domainY = 1.0

# Typically we run 12-16 in 3D, 30-60 in 2D (30,000 particles per partition)
cppx = 5   # cells per partition in each direction
cppy = 5
cppz = 3   # plane strain needs 3, 

pfw["xpar"]=3  # grid partitions
pfw["ypar"]=1
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] = -1.5 # mm
pfw["xmax"] =  1.5 # mm
pfw["ymin"] = -0.5 # mm
pfw["ymax"] = 0.5 # mm

# for plane strain we compute depth to equal grid cell size.
DZ = ( pfw["xmax"] - pfw["xmin"] ) / ( pfw["nI"] - 2 )
domainZ = DZ
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True    # run in background vs. in salloc interactive
pfw["mWallTime"] = "00:05:00"
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
sample = geom.box('block',
    x0 = [-0.5*DZ, -0.5*DZ, pfw["zmin"] ],
    x1 = [ 0.0, 0.0, pfw["zmax"] ],
    vel=np.array( [ initialVelocity, 0., 0. ] ),
    mat = 0,
    group = 0 )

pfw["objects"]=[sample]

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
# --- PFW VERIFICATION FAST DEBUG OVERRIDES BEGIN ---
# Debug-only runtime caps.  Keep this block below all source-file pfw assignments.
def _vv_fast_int(_value, _default):
    try:
        return int(float(str(_value).strip().strip('"').strip("'")))
    except Exception:
        return int(_default)

def _vv_fast_bool(_value):
    if isinstance(_value, bool):
        return _value
    return str(_value).strip().strip('"').strip("'").lower() in ("1", "true", "yes", "on")

try:
    refine = 1
except Exception:
    pass

# Fix common legacy typo before GEOS XML is written.
if "planeStrain" in pfw and "planeStrain" not in pfw:
    pfw["planeStrain"] = pfw.pop("planeStrain")

_vv_fast_plane = _vv_fast_bool(pfw.get("planeStrain", False))
# Treat thin 2D/plane-strain legacy cases as plane-like even when planeStrain was omitted.
try:
    if _vv_fast_int(pfw.get("zpar", 1), 1) == 1 and "nK" not in pfw:
        _vv_fast_plane = True
except Exception:
    pass

_vv_fast_cpp_cap = 24 if _vv_fast_plane else 8
_vv_fast_max_partitions = 2
pfw["mWallTime"] = "00:05:00"

for _vv_key in ("xpar", "ypar", "zpar"):
    pfw[_vv_key] = max(1, min(_vv_fast_int(pfw.get(_vv_key, 1), 1), _vv_fast_max_partitions))
if _vv_fast_plane:
    pfw["zpar"] = 1

# Preserve already coarser grids, but cap high cells-per-partition values.
def _vv_fast_cap_cells(_nkey, _pkey, _default_cells=1):
    _p = max(1, _vv_fast_int(pfw.get(_pkey, 1), 1))
    _n = _vv_fast_int(pfw.get(_nkey, 0), 0)
    if _n <= 0:
        return max(1, _p * min(_default_cells, _vv_fast_cpp_cap))
    _cpp = max(1, (_n + _p - 1) // _p)
    return max(1, _p * min(_cpp, _vv_fast_cpp_cap))

pfw["nI"] = _vv_fast_cap_cells("nI", "xpar", _vv_fast_cpp_cap)
pfw["nJ"] = _vv_fast_cap_cells("nJ", "ypar", _vv_fast_cpp_cap)
if _vv_fast_plane:
    if "nK" in pfw:
        pfw["nK"] = max(1, min(_vv_fast_int(pfw.get("nK", 1), 1), 8))
else:
    pfw["nK"] = _vv_fast_cap_cells("nK", "zpar", 8)

pfw["mCores"] = max(1, _vv_fast_int(pfw.get("xpar", 1), 1) * _vv_fast_int(pfw.get("ypar", 1), 1) * _vv_fast_int(pfw.get("zpar", 1), 1))
pfw["mNodes"] = max(1, (pfw["mCores"] + 111) // 112)
# --- PFW VERIFICATION FAST DEBUG OVERRIDES END ---
