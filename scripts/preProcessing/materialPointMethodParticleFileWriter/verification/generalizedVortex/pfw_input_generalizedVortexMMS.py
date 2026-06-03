# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 4

# DOMAIN ---------------------------------------------------------------------------------

pfw["planeStrain"] = 1

refine = 1 # partitions in each direction
cpp = 20   # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp
pfw["nJ"]=pfw["ypar"]*cpp
pfw["nK"]=3
pfw["ppc"]=2

domainHeight = 4
domainWidth = 4 

pfw["xmin"] = -0.5*domainHeight # mm
pfw["xmax"] = 0.5*domainHeight # mm
pfw["ymin"] =-0.5*domainWidth # mm
pfw["ymax"] = 0.5*domainWidth # mm

dx = domainWidth/(pfw["nI"]-2)
dy = domainHeight/(pfw["nJ"]-2)

domainLength = 3*dx

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# PATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 100
pfw["restartInterval"] = stopTime

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.05   
pfw["initialDt"]=1e-16

pfw["cpdiDomainScaling"] = 0
pfw["bodyForce"] = [ 0, 0, 0 ]
pfw["generalizedVortexMMS"] = 1
pfw["reactionHistory"] = 0
pfw["boxAverageHistory"] = 0

pfw["boundaryConditionTypes"] = [ 0, 0, 0, 0, 0, 0 ]
pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0  

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1000
firstLameConstant = 577
shearModulus = 384.615384615384585

pfw["materials"] = [ "hyperelasticMMS" ]
pfw["materialPropertyString"]="""
<HyperelasticMMS
	name="hyperelasticMMS"
	defaultDensity=""" + '"' + str(density) + '"' + """
	defaultLambda=""" + '"' + str(firstLameConstant) + '"' + """
	defaultShearModulus=""" + '"' + str(shearModulus) + '"' + """/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('test_specimen',[pfw["xmin"]+2*dx,pfw["ymin"]+2*dx,pfw["zmin"]],[pfw["xmax"]-2*dx,pfw["ymax"]-2*dx,pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]
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
