# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = True
stopTime = 30.0

# Domain ---------------------------------------------------------------------------------

sampleHeight = 1.0 # mm
sampleWidth = 1.0  # mm
sampleLength = 1.0 # mm 

domainHeight = sampleHeight
domainWidth = sampleWidth  # This would be increased for unconfined compression.
domainLength = sampleLength

cpp = 3   # cells per partition in each direction
refine = 1
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=1   # particles per cell in each direction

pfw["xmin"] = 0.0 # mm
pfw["xmax"] = domainHeight # mm
pfw["ymin"] =-0.5*domainWidth # mm
pfw["ymax"] = 0.5*domainWidth # mm
pfw["zmin"] =-0.5*domainWidth # mm
pfw["zmax"] = 0.5*domainWidth # mm

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"]=False

# END BATCH PARAMETERS ---------------------------------------------------------------


pfw["endTime"] = stopTime
pfw["plotInterval"] = 1.0
pfw["restartInterval"] = stopTime*5.0

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

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
							 "Damage",
							 "StrengthScale",
                             "SurfaceFlag",
                             "RVector"]

# END GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------

# GEOMETRIC OBJECTS ----------------------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
strengthScale = geom.strengthScaleWrapper('strengthScale', block, 1.0)
pfw["objects"]=[strengthScale]

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8e-32
thirdInvariantDependence=1

pfw["materials"] = ["sand"]
pfw["materialPropertyString"]="""
<CeramicDamage
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
	compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
	maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
	crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
	thirdInvariantDependence=""" + '"' + str(thirdInvariantDependence) + '"' + """/>
"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[2, 2, 2, 2, 2, 2]   

# For these parameters, the correct solution is for the stress to the yield surface
# yield with increasing pressure, then unload past zero shear stress to half the yield
# stress as the pressure goes to 0
shearStrain = maximumStrength/shear

p0=10.0
volStrain = -p0/bulk
J0 = np.exp(volStrain)
volStretch = J0**(1/3)

p02=-1.0
volStrain2 = -p02/bulk
J02 = np.exp(volStrain2)
volStretch2 = J02**(1/3)

pfw["fTableInterpType"] = "Linear"
pfw["prescribedBoundaryFTable"]=1

print("J0 = ",J0,", vol. strain = ",volStrain,", shearStrain=",shearStrain)

# Prescribed deformation table: pure shear load and unload.
pfw["fTable"]=[
[0,  1,	                                      1,	                                            1],
[10, volStretch,                              volStretch,                                       volStretch],
[11, volStretch*np.exp(shearStrain*1.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*1.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*1.0/10.0))],
[12, volStretch*np.exp(shearStrain*2.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*2.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*2.0/10.0))],
[13, volStretch*np.exp(shearStrain*3.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*3.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*3.0/10.0))],
[14, volStretch*np.exp(shearStrain*4.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*4.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*4.0/10.0))],
[15, volStretch*np.exp(shearStrain*5.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*5.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*5.0/10.0))],
[16, volStretch*np.exp(shearStrain*6.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*6.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*6.0/10.0))],
[17, volStretch*np.exp(shearStrain*7.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*7.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*7.0/10.0))],
[18, volStretch*np.exp(shearStrain*8.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*8.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*8.0/10.0))],
[19, volStretch*np.exp(shearStrain*9.0/10.0), volStretch/np.sqrt(np.exp(shearStrain*9.0/10.0)), volStretch/np.sqrt(np.exp(shearStrain*9.0/10.0))],
[20, volStretch*np.exp(shearStrain),          volStretch/np.sqrt(np.exp(shearStrain)),          volStretch/np.sqrt(np.exp(shearStrain))],
[30, volStretch2*np.exp(shearStrain),         volStretch2/np.sqrt(np.exp(shearStrain)),         volStretch2/np.sqrt(np.exp(shearStrain))]
]
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
