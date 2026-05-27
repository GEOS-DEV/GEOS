# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {} 
pfw["runDebug"] = False
stopTime = 5.

innerDiameter = 2.*(7.96) - 0.25
outerDiameter = 2*(7.96) + 0.25
initialRadialVelocity = 0.016 #0.01592

density = 2.75						# mg/mm^3
bulk = 275./( 3.*( 1. - 2.*0.3) )	# GPa
shear = 275./( 2.*( 1. + 0.3) )		# GPa
failureStrength = 0.3 #0.30				# GPa
energyReleaseRate = 1e-05 #1.0e-5 			# mg/us^2

weibullModulus = 6.



# Domain ---------------------------------------------------------------------------------
domainX = 1.1*outerDiameter
domainY = 1.1*outerDiameter

cppx = 8   # cells per partition in each direction
cppy = 8
cppz = 5

refine = 1

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cppx  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cppy  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cppz  # grid cells in the z-direction
pfw["ppc"]=2   # particles per cell in each direction

pfw["xmin"] =-0.5*domainX # mm
pfw["xmax"] = 0.5*domainX # mm
pfw["ymin"] =-0.5*domainY # mm
pfw["ymax"] = 0.5*domainY # mm

DZ = ( pfw["xmax"] - pfw["xmin"] ) / ( pfw["nI"] - 2 )
domainZ = 3.*DZ
pfw["zmin"] =-0.5*domainZ # mm
pfw["zmax"] = 0.5*domainZ # mm

# BATCH PARAMETERS --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"] = "00:05:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=False
pfw["autoRestart"]=False

# END BATCH PARAMETERS ---------------------------------------------------------------
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 100
pfw["restartInterval"] = stopTime*5.0

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------
pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1
pfw["separabilityMinDamage"]=0.5               
pfw["treatFullyDamagedAsSingleField"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=0
pfw["boxAverageHistory"]=0

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0
pfw["FSubcycles"]=10

pfw["plotUnscaledParticles"]=0

pfw["frictionCoefficient"]=0.05

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]="2"

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
							 "StrengthScale",
							 "SurfaceFlag",
                             "RVector"]

# END GEOSX MPM SOLVER PARAMETERS ---------------------------------------------------------------

# GEOMETRY OBJECTS -------------------------------------------------------------------------------------------

# Radially varying velocity for cylinder:
def getVelocity(self,pt):
    x = np.array(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point
    if ( r > 1.e-12):
      er = xr/r
      vel= initialRadialVelocity*er
    else:
      vel= np.array([0.,0.,0.])
    return vel
def make_objects():
    ring = geom.cylinder('ring',
        x1=[0.0,0.0,-0.5*DZ],
        x2=[0.0,0.0,0.5*DZ],
        ri=0.5*innerDiameter,
        r=0.5*outerDiameter,
        vel=getVelocity,
        mat=0,
        group=0)

    # weibullVolume = DZ*0.25*np.pi*( outerDiameter**2 - innerDiameter**2 )
    # weibullring = geom.voronoiWeibullBoxWrapper('weibullring',
    #     subObject=ring,
    #     x0=[pfw["xmin"],pfw["ymin"],pfw["zmin"]],
    #     x1=[pfw["xmax"],pfw["ymax"],pfw["zmax"]],
    #     flawSize=5*DZ,
    #     weibullVolume=weibullVolume,
    #     weibullModulus=weibullModulus,
    #     weibullSeed=1,
    #     vMin=(DZ)**3.
    #     )

    objects=[ring]
    return objects

# MATERIAL PROPERTIES --------------------------------------------------------------------



pfw["materials"] = ["sand"]
pfw["materialPropertyString"]="""
<Chiumenti
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	failureStrength="""+'"'+str(failureStrength)+'"'+"""
	energyReleaseRate="""+'"'+str(energyReleaseRate)+'"'+"""/>
"""

# DEFORMATION ---------------------------------------------------------------------------------

pfw["boundaryConditionTypes"]=[0, 0, 0, 0, 0, 0]
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
