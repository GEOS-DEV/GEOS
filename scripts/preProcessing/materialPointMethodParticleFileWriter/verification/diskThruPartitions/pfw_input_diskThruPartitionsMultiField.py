# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# START!

# Domain ---------------------------------------------------------------------------------
refine = 1
xpar=refine  # grid partitions
ypar=refine
zpar=1

cpp = 8
nI=xpar*cpp  	# grid cells in the x-direction
nJ=ypar*cpp  	# grid cells in the y-direction
nK=3  			# grid cells in the z-direction
ppc=2   		# particles per cell in each direction

domainLength = 1.0 # m
domainThickness = domainLength*(nK-2)/(nI-2)  # m, to get cubic cells

# Define all the geometric objects -------------------------------------------------------
xmin =-0.5*domainLength	# m
xmax = 0.5*domainLength	# m
ymin =-0.5*domainLength	# m
ymax = 0.5*domainLength	# m
zmin =-0.5*domainThickness 	# m
zmax = 0.5*domainThickness 	# m

lx = xmax - xmin
ly = ymax - ymin
lz = zmax - zmin

# pfw flags
planeStrain = 1
useDamageAsSurfaceFlag = 1

# Batch parameters for GEOS runs.  --------------------------------------------------------
# read in the default bank:
import sys
import importlib
import getpass
username = getpass.getuser()
userDefsFile = str('userDefs_'+str(username))
userDefs = importlib.import_module(userDefsFile)
mBank = userDefs.defaultBank

mWallTime = "00:05:00"
mCores=xpar*ypar*zpar
mNodes=int(np.ceil(float(mCores)/36.)) 
mSubmitJobs=False

# GEOSX MPM input parameters ---------------------------------------------------------------
endTime="0.01"				# seconds
writePlot="1"				# this does nothing for now
writeRestart="1"			# this does nothing for now
plotInterval="0.00005"
restartInterval="0.0025"

# specify an array with all objects to be included, order matters. for overlapping objects, the first one listed will be assigned at each point.
# "fill" must be last on the list.

disk1 = geom.cylinder('disk1',[xmax/2,ymax/2,zmin-lz],[xmax/2,ymax/2,zmax+lz],r=0.15,vel=[-120.0,-120.0,0.0],mat=0,group=0)
disk2 = geom.cylinder('disk2',[-xmax/2,-ymax/2,zmin-lz],[-xmax/2,-ymax/2,zmax+lz],r=0.15,vel=[12.0,12.0,0.0],mat=1,group=1)

objects=[disk1,disk2]

# Material Properties
# Notes:
# 1. It DOES NOT matter if the 'materials' array matches the order of 'materialPropertyString'.
# 2. It's okay if there are more materials specified than used.
# 3. The 'materials' array is a map from material name to material ID in the particle file.
#    i.e. if the first material in 'materials' is 'sand', then 'sand' is material ID 0.
#    This matters when using integers to flag objects with materials.
materials = [ "aluminum", "aluminum10x" ]
materialPropertyString="""
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

mpmSolverParameterString="""
timeIntegrationOption="ExplicitDynamic"
cflFactor="0.5"    
initialDt="1e-16"

prescribedBcTable="0"
prescribedBoundaryFTable="0"
fTableInterpType="0"    

solverProfiling="0"

planeStrain="""+"\""+str(planeStrain)+"\""+"""

neighborRadius="-1.01"
needsNeighborList="0"

cpdiDomainScaling="1"

damageFieldPartitioning="0"

contactGapCorrection="Simple"
frictionCoefficient="0.25"

boundaryConditionTypes="{ 0, 0, 0, 0, 0, 0 }"    
"""

# New pfw dictionary interface -----------------------------------------------------------
pfw = {
    "runDebug": True,
    "xpar": xpar,
    "ypar": ypar,
    "zpar": zpar,
    "nI": nI,
    "nJ": nJ,
    "nK": nK,
    "ppc": ppc,
    "xmin": xmin,
    "xmax": xmax,
    "ymin": ymin,
    "ymax": ymax,
    "zmin": zmin,
    "zmax": zmax,
    "planeStrain": planeStrain,
    "mBatch": True,
    "mBank": mBank,
    "mWallTime": mWallTime,
    "mCores": mCores,
    "mNodes": mNodes,
    "mSubmitJobs": mSubmitJobs,
    "endTime": float(endTime),
    "plotInterval": float(plotInterval),
    "restartInterval": float(restartInterval),
    "objects": objects,
    "materials": materials,
    "materialPropertyString": materialPropertyString,
    "timeIntegrationOption": "ExplicitDynamic",
    "cflFactor": 0.5,
    "initialDt": 1e-16,
    "prescribedBcTable": 0,
    "prescribedBoundaryFTable": 0,
    "fTableInterpType": 0,
    "solverProfiling": 1,
    "neighborRadius": -1.01,
    "needsNeighborList": 0,
    "cpdiDomainScaling": 1,
    "damageFieldPartitioning": 0,
    "contactGapCorrection": "Simple",
    "frictionCoefficient": 0.25,
    "boundaryConditionTypes": [0, 0, 0, 0, 0, 0],
}
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
