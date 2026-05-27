# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = True
stopTime = 0.01

# DOMAIN ---------------------------------------------------------------------------------

refine = 1 # partitions in each direction
cpp = 8 # grid cells per direction for each partition

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  	# grid cells in the y-direction
pfw["nK"]=pfw["ypar"]*cpp  	# grid cells in the z-direction
pfw["ppc"]=2   		        # particles per cell in each direction

domainLength = 1.0 # m
domainThickness = domainLength*(pfw["nK"]-2)/(pfw["nI"]-2)  # m, to get cubic cells

pfw["xmin"] =-0.5*domainLength	   # m
pfw["xmax"] = 0.5*domainLength	   # m
pfw["ymin"] =-0.5*domainLength	   # m
pfw["ymax"] = 0.5*domainLength	   # m
pfw["zmin"] =-0.5*domainThickness  # m
pfw["zmax"] = 0.5*domainThickness  # m

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS ------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/1000.0001
pfw["restartInterval"]=stopTime/4

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25    
pfw["initialDt"]=1e-16

pfw["contactGapCorrection"] = "Implicit"
pfw["frictionCoefficient"]=0.0 

# MATERIAL PROPERTIES ---------------------------------------------------------

pfw["materials"] = [ "aluminum", "aluminum10x" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2700"
	defaultBulkModulus="70.0e8"
	defaultShearModulus="24.0e8"/>
<ElasticIsotropic
	name="aluminum10x"
	defaultDensity="27000"
	defaultBulkModulus="70.0e9"
	defaultShearModulus="24.0e9"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

ball1 = geom.sphere('ball1',[pfw["xmax"]/2,pfw["ymax"]/2,pfw["zmax"]/2],0.15,vel=[-60.0,-60.0,-60.0],mat=0,group=0)
ball2 = geom.sphere('ball2',[-pfw["xmax"]/2,-pfw["ymax"]/2,-pfw["zmax"]/2],0.15,vel=[6.0,6.0,6.0],mat=1,group=1)
pfw["objects"]=[ball1,ball2]
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
# --- PFW VERIFICATION FAST DEBUG OVERRIDES END ---
