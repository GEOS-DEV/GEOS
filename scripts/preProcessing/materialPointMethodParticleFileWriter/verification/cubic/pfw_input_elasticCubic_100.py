# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 10.0

# DOMAIN ---------------------------------------------------------------------------------

sampleHeight = 1.0 # mm
sampleWidth = 1.0  # mm
sampleLength = 1.0 # mm 

domainHeight = sampleHeight
domainWidth = sampleWidth
domainLength = sampleLength

refine = 1 # partitions in each direction
cpp = 5   # cells per partition in each direction

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=1   # particles per cell in each direction

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] =  0.5*domainWidth # mm
pfw["ymin"] =  0.0 # mm
pfw["ymax"] =  domainHeight # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"] = stopTime
pfw["plotInterval"] = 1.0
pfw["restartInterval"] = stopTime*5.0

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.05    
pfw["initialDt"]=1e-16

pfw["solverProfiling"]=0
pfw["planeStrain"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["cpdiDomainScaling"]=0
pfw["damageFieldPartitioning"]=0

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
							 "MaterialDirection"]

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-sampleWidth/2, 0.0, -sampleLength/2],[sampleWidth/2, sampleHeight, sampleLength/2],vel=[0.0,0.0,0.0],mat=0,group=0)
blockWDir = geom.materialDirectionWrapper('blockWDir',block,[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
pfw["objects"]=[blockWDir]


# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.329
c11 = 165.64
c12 = 63.94
c44 = 79.51

pfw["materials"] = ["silicon"]
pfw["materialPropertyString"]="""
<ElasticCubic
	name="silicon"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultC11="""+'"'+str(c11)+'"'+"""
	defaultC12="""+'"'+str(c12)+'"'+"""
    defaultC44="""+'"'+str(c44)+'"'+"""/>
"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[0, 0, 2, 2, 0, 0]   

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[
[0,  1.0, 1.0,	1.0],
[10, 1.0, 1.10, 1.0]]
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
