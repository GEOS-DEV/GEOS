import pfw_geometryObjects as geom   
import numpy as np                   
from sklearn.neighbors import KDTree 

pfw = {} 
pfw["runDebug"] = False
stopTime = 10

# DOMAIN ---------------------------------------------------------------------------------

refine = 1 # partitions per direction
cpp = 8   # cells per partition in each direction

domainHeight = 2.0
domainWidth = 2.0
domainLength = 2.0

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth  # mm
pfw["ymin"] = 0.0              # mm
pfw["ymax"] = domainHeight     # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

pfw["planeStrain"] = 0
pfw["periodic"] = [False, False, False]

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

ppvc=10               # average number of particles across voronoi cell.

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8
thirdInvariantDependence = 1
damagedMaterialFrictionalSlope = 0.10

pfw["materials"] = ["sand"] # Only include constitutive models for the particles in this list in the order they appear below
# Include both particle and cohesive zone models here
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
<PolymerCohesiveZone
    name="cz1"
    thickness="0.01"
	bulkModulus="1.0"
	shearModulus="0.25"
	yieldStrength0="0.005"
	r0="0.025"
	r1="0.01"
	r2="0.15"
	Gr="0.002"
	maxStretch="3.0"/>
"""

# BATCH PARAMETERS  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime             
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20
pfw["solverProfiling"]=0

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["needsNeighborList"]=1
pfw["damageFieldPartitioning"]=1
pfw["treatFullyDamagedAsSingleField"]=1
pfw["resetDefGradForFullyDamagedParticles"]=1
pfw["plotUnscaledParticles"]=1
pfw["FSubcycles"]=10

pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["frictionCoefficient"]=0.01
pfw["contactGapCorrection"] = "Implicit"
pfw["explicitSurfaceNormalInfluence"]= 1000.0
pfw["useSurfacePositionForContact"]= 1
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"]=1

pfw["maxParticleVelocity"]=10.0               
pfw["minParticleJacobian"]=0.01                       
pfw["maxParticleJacobian"]=10.0  

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

pfw["useEvents"]=1

pfw["preventCZInterpenetration"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
							 "StrengthScale",
                             "SurfaceFlag",
                             "RVector",
							 "SurfaceNormal",
							 "SurfacePosition"]

# GEOMETRY OBJECTS -------------------------------------------------------
grainDiameter = 0.225
radialBias = 1.0
seed = 1436359
porosity = 0.0

dx = domainHeight/(pfw["nI"]-2)/pfw["ppc"]
dy = domainWidth/(pfw["nJ"]-2)/pfw["ppc"]
dz = domainWidth/(pfw["nK"]-2)/pfw["ppc"]

# Weibull variability is applied per each granule
vMin=dx*dy*dz
flawSize = dx*ppvc     # average diameter of voronoi cell.

weibullModulus = 5
weibullVolume = 100

radius = 0.75
height = domainHeight

def make_objects():    
    tesselation=geom.czPrill("tesselation", [-radius, 0.0, -radius], [radius, height, radius], grainDiameter=grainDiameter, porosity=porosity, radialBias=radialBias, seed=seed, vel=[0.0,0.0,0.0], mat=0, group=0, dim=3)
    cylinder=geom.cylinder("cylinder", [0.0, 0.0, 0.0], [0.0, height, 0.0], radius, vel=[0.0,0.0,0.0], mat=0, group=0)
    czCylinder=geom.intersection("czCylinder", tesselation, cylinder)
    czCylinderWStrength = geom.voronoiWeibullBoxWrapper("prill1WStrength", czCylinder, [-radius, 0.0, -radius], [radius, height, radius], flawSize, weibullModulus, weibullVolume, seed, vMin, vpts=tesselation.vpts, dim=3)

    return [czCylinderWStrength]

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 0, 0 ]

pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,            1.0,  1.0,    1.0],  
               [stopTime,     1.0,  0.9,    1.0]]

pfw["mpmEventsString"]="""
<CohesiveZone 
    name="cz1" 
    startTime=
""" + '"' + str(0.0) + '"' + """
    constitutiveModel="cz1"
    czVolumeNormalization="1"/>
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
# --- PFW VERIFICATION FAST DEBUG OVERRIDES END ---
