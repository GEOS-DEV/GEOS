# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 10.0

# DOMAIN ---------------------------------------------------------------------------------

domainWidth = 1.0
domainHeight = 1.0
domainLength = 1.0

pfw["xmin"] = -0.5*domainWidth  # mm
pfw["xmax"] =  0.5*domainWidth  # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] =  0.5*domainHeight # mm
pfw["zmin"] = -0.5*domainLength # mm
pfw["zmax"] =  0.5*domainLength # mm

refine = 1 # partitions in each direction
cpp = 8   # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

dx = domainWidth / (pfw["nI"]-1)
dy = domainHeight / (pfw["nJ"]-1)
dz = domainLength / (pfw["nK"]-1)

# BATCH PARAMETERS --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1

pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1

pfw["minParticleJacobian"]=1e-8
pfw["maxParticleJacobian"]=1e12
pfw["subdivideParticles"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "SurfaceFlag",
                             "RVector",
                             "SurfaceNormal",
                             "SurfacePosition",
							 "SurfaceTraction"]

pfw["plottableFields"]=["particleID",
                        "particleMass",
                        "particleVolume",
                        "particleDamage",
                        "particleInitialVolume",
                        "particleDensity",
                        "particleMaterialType",
                        "particleGroup",
                        "particleSurfaceFlag",
                        "particleStrengthScale",
                        "particleCenter",
                        "particleReferencePosition",
                        "particleVelocity",
                        "particleStress",
                        "particleSurfaceNormal",
                        "particleSurfacePosition",
                        "particleMaterialDirection",
                        "particleDeformationGradient",
                        "particlePlasticStrain",
                        "particleCohesiveZoneFlag",
                        "particleCohesiveForce",
                        "particleRVectors",
                        "particleReferenceRVectors",
                        "gridMass",
                        "gridDamage",
                        "gridPosition",
                        "gridVelocity",
                        "gridAcceleration",
                        "gridSurfaceNormal",
                        "gridSurfacePosition",
                        "gridContactForce",
                        "gridCohesiveNode"]

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "gas" ]
pfw["materialPropertyString"] ="""
<Gas
	name="gas"
	defaultDensity="0.1"
	referencePressure="0.101325"
    referenceTemperature="325"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-dx/2,-dy/2,-dz/2],[dx/2,dy/2,dz/2],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 1, 1, 1, 1, 1, 1 ]
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
