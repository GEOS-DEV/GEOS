# -*- coding: utf-8 -*-
"""Full-3D planar weak-interface smoke test.

This is intentionally smaller and simpler than the 2D verification problems.
Its purpose is to exercise the prescribed particle-surface C0 mortar path in
full 3D after applying the 3D fallback patch.  The two elastic blocks fill the
same Cartesian domain and are split by an oblique analytic plane.  In the trace
variant, particles near the split plane carry WeakDiscontinuity surface flags,
opposing contact groups, and explicit surface normals/positions.
"""

from __future__ import annotations

import math
import os

import numpy as np

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers


pfw = {}

CASE_ID = "weakTracePlanarInterface3D"
VARIANT = os.environ.get("WEAK_TRACE_VARIANT", "singleField")
VARIANT_LABEL = os.environ.get("VV_VARIANT_LABEL", VARIANT)
CASE_NAME = os.environ.get("VV_CASE_NAME", f"{CASE_ID}_{VARIANT}")


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


def elastic_xml(name: str, density: float, young: float, poisson: float) -> str:
    return f"""
<ElasticIsotropic
    name="{name}"
    defaultDensity="{density}"
    defaultYoungModulus="{young}"
    defaultPoissonRatio="{poisson}"/>"""


class HalfSpaceBox(geom.Geometry):
    """Cartesian box clipped by one analytic plane.

    PFW asks each geometry object whether a particle center is inside and, for
    surface-flagged particles, asks for a surface normal and position.  Keeping
    the signed-plane calculation here ensures refined child particles recompute
    their own distance to the interface instead of inheriting a parent offset.
    """

    def __init__(
        self,
        name: str,
        bounds_min,
        bounds_max,
        plane_normal,
        plane_offset: float,
        negative_side: bool,
        flag_interface: bool,
        vel,
        mat: int,
        group: int,
        particle_type: int,
    ):
        super().__init__(name, vel=vel, mat=mat, group=group, particleType=particle_type)
        self.bounds_min = np.asarray(bounds_min, dtype=float)
        self.bounds_max = np.asarray(bounds_max, dtype=float)
        normal = np.asarray(plane_normal, dtype=float)
        normal_norm = float(np.linalg.norm(normal))
        if normal_norm <= 0.0:
            raise ValueError("plane_normal must be nonzero")
        self.normal = normal / normal_norm
        self.offset = float(plane_offset)
        self.negative_side = bool(negative_side)
        self.flag_interface = bool(flag_interface)

    def _signed_distance(self, pt) -> float:
        return float(np.dot(self.normal, np.asarray(pt, dtype=float)) - self.offset)

    def _inside_box(self, pt) -> bool:
        x = np.asarray(pt, dtype=float)
        return bool(np.all(x >= self.bounds_min) and np.all(x < self.bounds_max))

    def isInterior(self, pt, skinDepth):
        if not self._inside_box(pt):
            return -1
        signed = self._signed_distance(pt)
        if self.negative_side:
            if signed > 0.0:
                return -1
        elif signed < 0.0:
            return -1
        if self.flag_interface and abs(signed) <= float(skinDepth):
            return int(geom.SurfaceFlag.WeakDiscontinuity)
        return 0

    def getSurfaceNormal(self, pt):
        # Outward normal from each material subdomain across the internal split.
        return self.normal.copy() if self.negative_side else -self.normal.copy()

    def getSurfacePosition(self, pt):
        # Vector from particle center to its closest point on the split plane.
        return -self._signed_distance(pt) * self.normal

    def xMin(self):
        return float(self.bounds_min[0])

    def xMax(self):
        return float(self.bounds_max[0])

    def getSubregions(self):
        return [(self.mat, self.particleType)]


pfw["caseName"] = CASE_NAME
pfw["runDebug"] = env_bool("WEAK_TRACE_RUN_DEBUG", True)
pfw["mPartition"] = os.environ.get("WEAK_TRACE_PARTITION", "pdebug" if pfw["runDebug"] else "pbatch")
pfw["mBatch"] = True
pfw["mWallTime"] = os.environ.get("WEAK_TRACE_WALLTIME", "00:10:00")
pfw["mSubmitJobs"] = False

# Domain and discretization -------------------------------------------------------------
refine = int(os.environ.get("WEAK_TRACE_REFINE", "1"))
cpp = int(os.environ.get("WEAK_TRACE_CPP", "14"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", "1"))
surface_refinement_level = int(os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_LEVEL", "1"))

pfw["planeStrain"] = 0
pfw["periodic"] = [False, False, False]
pfw["xpar"] = int(os.environ.get("WEAK_TRACE_XPAR", str(refine)))
pfw["ypar"] = int(os.environ.get("WEAK_TRACE_YPAR", str(refine)))
pfw["zpar"] = int(os.environ.get("WEAK_TRACE_ZPAR", str(refine)))
pfw["nI"] = int(os.environ.get("WEAK_TRACE_NI", str(pfw["xpar"] * cpp)))
pfw["nJ"] = int(os.environ.get("WEAK_TRACE_NJ", str(pfw["ypar"] * cpp)))
pfw["nK"] = int(os.environ.get("WEAK_TRACE_NK", str(pfw["zpar"] * cpp)))
if pfw["nI"] % pfw["xpar"] != 0 or pfw["nJ"] % pfw["ypar"] != 0 or pfw["nK"] % pfw["zpar"] != 0:
    raise ValueError("WEAK_TRACE_NI/NJ/NK must be divisible by WEAK_TRACE_XPAR/YPAR/ZPAR")
pfw["ppc"] = ppc
pfw["surfaceRefinementLevel"] = surface_refinement_level

pfw["xmin"] = -0.5
pfw["xmax"] = 0.5
pfw["ymin"] = -0.35
pfw["ymax"] = 0.35
pfw["zmin"] = -0.30
pfw["zmax"] = 0.30

# Time integration and output -----------------------------------------------------------
stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "0.60"))
final_strain = float(os.environ.get("WEAK_TRACE_FINAL_STRAIN", "0.002"))
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 20.0
pfw["restartInterval"] = 10.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = float(os.environ.get("WEAK_TRACE_CFL", "0.18"))
pfw["initialDt"] = 1.0e-16
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 200.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 200.0
pfw["frictionCoefficient"] = 0.0
pfw["contactNormalType"] = "Aligned"
pfw["useSurfacePositionForContact"] = 1
pfw["damageFieldPartitioning"] = 0
pfw["weakDiscontinuityInsertion"] = "off"
pfw["weakInterfaceTraceGeometry"] = "surface"
pfw["cpdiDomainScaling"] = 0
pfw["updateMethod"] = os.environ.get("WEAK_TRACE_UPDATE_METHOD", "XPIC")
default_update_order = "2" if str(pfw["updateMethod"]).upper() == "XPIC" else "1"
pfw["updateOrder"] = int(os.environ.get("WEAK_TRACE_UPDATE_ORDER", default_update_order))
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = [
    "gridMass",
    "gridVelocity",
    "gridActive",
    "gridSurfaceFieldMass",
    "gridSurfaceNormal",
    "gridExplicitSurfaceNormal",
    "gridPrincipalExplicitSurfaceNormal",
    "gridSurfacePosition",
    "gridDisplacement",
    "gridWeakInterfaceTraceActive",
    "gridWeakInterfaceTraceContactSuppressed",
    "gridWeakInterfaceTraceSkipReason",
    "gridWeakInterfaceTraceAnchorWeight",
    "gridWeakInterfaceTraceSupportWeight",
    "gridWeakInterfaceTraceSurfaceJump",
    "gridWeakInterfaceTraceVelocityJump",
    "gridWeakInterfaceTraceVelocityJumpPost",
    "gridWeakInterfaceTraceForce",
    "gridInterfaceMechanism",
    "gridInterfaceNormalForce",
    "gridInterfaceTangentialForce",
    "gridInterfaceNormalVelocityJump",
    "gridInterfaceTangentialVelocityJump",
    "gridInterfaceNormalDisplacementJump",
    "gridInterfaceTangentialDisplacementJump",
]
pfw["plottableFields"] = [
    "particleID",
    "particleRank",
    "ghostRank",
    "particleMass",
    "particleVolume",
    "particleDensity",
    "particleMaterialType",
    "particleGroup",
    "particleSurfaceFlag",
    "particleCenter",
    "particleReferencePosition",
    "particleVelocity",
    "particleDisplacement",
    "particleStress",
    "particleSurfaceNormal",
    "particleSurfacePosition",
    "particleDeformationGradient",
    "particleRVectors",
    "particleReferenceRVectors",
]
pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "SurfaceFlag",
    "RVector",
    "SurfaceNormal",
    "SurfacePosition",
]

# Move only the x faces.  The oblique internal plane then cuts the active 3D
# cells while the boundary condition remains simple and robust for a smoke test.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 0, 0, 0, 0]
pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0 + final_strain, 1.0, 1.0]]

if VARIANT == "traceContactGroups":
    pfw["enableContact"] = 1
    pfw["enableWeakInterfaceTraceProjection"] = 1
    pfw["explicitSurfaceNormalInfluence"] = float(os.environ.get("WEAK_TRACE_EXPLICIT_NORMAL_INFLUENCE", "25.0"))
    pfw["weakInterfaceTraceProjectionIterations"] = int(os.environ.get("WEAK_TRACE_PROJECTION_ITERATIONS", "100"))
    pfw["weakInterfaceTraceProjectionScale"] = float(os.environ.get("WEAK_TRACE_PROJECTION_SCALE", "1.0"))
    pfw["weakInterfaceTraceGapStabilization"] = float(os.environ.get("WEAK_TRACE_GAP_STABILIZATION", "0.0"))
    pfw["weakInterfaceTraceMaxGapCorrectionVelocity"] = float(os.environ.get("WEAK_TRACE_MAX_GAP_CORRECTION_VELOCITY", "0.0"))
    pfw["weakInterfaceTraceMinWeight"] = float(os.environ.get("WEAK_TRACE_MIN_WEIGHT", "1.0e-12"))
    pfw["weakInterfaceTraceCutTolerance"] = float(os.environ.get("WEAK_TRACE_CUT_TOLERANCE", "1.0e-8"))
    pfw["weakInterfaceTraceMinSegmentFraction"] = float(os.environ.get("WEAK_TRACE_MIN_SEGMENT_FRACTION", "0.02"))
    pfw["weakInterfaceTraceMortarQuadratureOrder"] = int(os.environ.get("WEAK_TRACE_MORTAR_QUADRATURE_ORDER", "2"))
    pfw["weakInterfaceTraceSurfaceFitNeighbors"] = int(os.environ.get("WEAK_TRACE_SURFACE_FIT_NEIGHBORS", "8"))
    pfw["weakInterfaceTraceSurfaceFitRadius"] = float(os.environ.get("WEAK_TRACE_SURFACE_FIT_RADIUS", "10.0"))
    pfw["weakInterfaceTraceMinNormalCoherence"] = float(os.environ.get("WEAK_TRACE_MIN_NORMAL_COHERENCE", "0.0"))
    pfw["weakInterfaceTraceMaxFitResidual"] = float(os.environ.get("WEAK_TRACE_MAX_FIT_RESIDUAL", "1.0e6"))
    pfw["weakInterfaceTraceMinFacetAreaFraction"] = float(os.environ.get("WEAK_TRACE_MIN_FACET_AREA_FRACTION", "0.0"))
    pfw["weakInterfaceTraceMinTriangleQuality"] = float(os.environ.get("WEAK_TRACE_MIN_TRIANGLE_QUALITY", "0.0"))
    pfw["weakInterfaceTraceFallbackToSingleField"] = int(os.environ.get("WEAK_TRACE_FALLBACK_TO_SINGLE_FIELD", "1"))
    pfw["weakInterfaceTraceLinearSolverTolerance"] = float(os.environ.get("WEAK_TRACE_LINEAR_SOLVER_TOLERANCE", "1.0e-10"))
    pfw["weakInterfaceTraceLinearSolverRegularization"] = float(os.environ.get("WEAK_TRACE_LINEAR_SOLVER_REGULARIZATION", "1.0e-12"))
    pfw["weakInterfaceTraceSuppressNodalContact"] = 1
    pfw["weakInterfaceTracePairs"] = [[0, 1]]
else:
    pfw["enableContact"] = 0
    pfw["enableWeakInterfaceTraceProjection"] = 0

# Materials -----------------------------------------------------------------------------
left_name = "weakTrace3DCompliant"
right_name = "weakTrace3DStiff"
density = float(os.environ.get("WEAK_TRACE_DENSITY", "1.0"))
poisson = float(os.environ.get("WEAK_TRACE_POISSON", "0.25"))
left_E = float(os.environ.get("WEAK_TRACE_LEFT_E", "2.0"))
right_E = float(os.environ.get("WEAK_TRACE_RIGHT_E", "8.0"))
pfw["materials"] = [left_name, right_name]
pfw["materialPropertyString"] = elastic_xml(left_name, density, left_E, poisson) + "\n" + elastic_xml(right_name, density, right_E, poisson)

# Geometry -------------------------------------------------------------------------------
plane_normal = np.array([1.0, 0.35, 0.25], dtype=float)
plane_normal /= float(np.linalg.norm(plane_normal))
plane_offset = float(os.environ.get("WEAK_TRACE_INTERFACE_OFFSET", "0.025"))
bounds_min = [pfw["xmin"], pfw["ymin"], pfw["zmin"]]
bounds_max = [pfw["xmax"], pfw["ymax"], pfw["zmax"]]
trace_groups = VARIANT == "traceContactGroups"

pfw["objects"] = [
    HalfSpaceBox(
        "negative_halfspace",
        bounds_min,
        bounds_max,
        plane_normal,
        plane_offset,
        negative_side=True,
        flag_interface=trace_groups,
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        particle_type=particle_type,
    ),
    HalfSpaceBox(
        "positive_halfspace",
        bounds_min,
        bounds_max,
        plane_normal,
        plane_offset,
        negative_side=False,
        flag_interface=trace_groups,
        vel=[0.0, 0.0, 0.0],
        mat=1,
        group=1 if trace_groups else 0,
        particle_type=particle_type,
    ),
]

# Tracers sample both materials close to the interface plus two far-field points.
tracer_points = []
for signed_offset in [-0.08, -0.025, 0.025, 0.08]:
    point = plane_normal * (plane_offset + signed_offset)
    tracer_points.append((float(point[0]), float(point[1]), float(point[2])))
tracer_points.extend([(-0.35, -0.20, -0.15), (0.35, 0.20, 0.15)])

tracers.set_tracers(
    pfw,
    tracer_points,
    variables=[
        "particleID",
        "materialType",
        "density",
        "stressXX",
        "stressYY",
        "stressZZ",
        "stressYZ",
        "stressXZ",
        "stressXY",
        "velocityX",
        "velocityY",
        "velocityZ",
    ],
    write_interval=stop_time / 200.0,
    output_prefix="weakTracePlanarInterface3DTracer",
)

pfw_expected = {
    "case_id": CASE_ID,
    "variant": VARIANT,
    "variant_label": VARIANT_LABEL,
    "stop_time": stop_time,
    "final_strain": final_strain,
    "plane_strain": pfw["planeStrain"],
    "interface_normal": [float(v) for v in plane_normal],
    "interface_offset": plane_offset,
    "particle_type": particle_type,
    "update_method": pfw["updateMethod"],
    "update_order": pfw["updateOrder"],
    "surface_refinement_level": surface_refinement_level,
    "xpar": pfw["xpar"],
    "ypar": pfw["ypar"],
    "zpar": pfw["zpar"],
    "nI": pfw["nI"],
    "nJ": pfw["nJ"],
    "nK": pfw["nK"],
}
