# -*- coding: utf-8 -*-
"""Full-3D diagonal-bar weak-interface test.

The sample is a square prism running along the long diagonal of a Cartesian
box.  Equal F-table stretch is applied on the x, y, and z moving faces, so the
opposite corner tabs separate along the diagonal.  The internal material split
is normal to the diagonal and located near the box center.  The prism
cross-section is small enough, and the box large enough, that the mortar
interface and its free side edges remain several cells away from the
domain-boundary nodes.
"""

from __future__ import annotations

import math
import os

import numpy as np

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers


pfw = {}

CASE_ID = "weakTraceRotatedBar3D"
VARIANT = os.environ.get("WEAK_TRACE_VARIANT", "singleField")
VARIANT_LABEL = os.environ.get("VV_VARIANT_LABEL", VARIANT)
CASE_NAME = os.environ.get("VV_CASE_NAME", f"{CASE_ID}_{VARIANT}")


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


def unit(v) -> np.ndarray:
    arr = np.asarray(v, dtype=float)
    norm = float(np.linalg.norm(arr))
    if norm <= 0.0:
        raise ValueError("cannot normalize a zero vector")
    return arr / norm


def elastic_xml(name: str, density: float, young: float, poisson: float) -> str:
    return f"""
<ElasticIsotropic
    name="{name}"
    defaultDensity="{density}"
    defaultYoungModulus="{young}"
    defaultPoissonRatio="{poisson}"/>"""


class DiagonalPrismHalfSpace(geom.Geometry):
    """Square prism along a diagonal axis, split by one normal plane."""

    def __init__(
        self,
        name: str,
        axis,
        transverse_1,
        transverse_2,
        half_width_1: float,
        half_width_2: float,
        plane_offset: float,
        negative_side: bool,
        flag_interface: bool,
        bounds_min,
        bounds_max,
        vel,
        mat: int,
        group: int,
        particle_type: int,
    ):
        super().__init__(name, vel=vel, mat=mat, group=group, particleType=particle_type)
        self.axis = unit(axis)
        self.t1 = unit(transverse_1)
        self.t2 = unit(transverse_2)
        self.half_width_1 = float(half_width_1)
        self.half_width_2 = float(half_width_2)
        self.offset = float(plane_offset)
        self.negative_side = bool(negative_side)
        self.flag_interface = bool(flag_interface)
        self.bounds_min = np.asarray(bounds_min, dtype=float)
        self.bounds_max = np.asarray(bounds_max, dtype=float)
        if self.half_width_1 <= 0.0 or self.half_width_2 <= 0.0:
            raise ValueError("diagonal prism half widths must be positive")

    def _local(self, pt) -> tuple[float, float, float]:
        x = np.asarray(pt, dtype=float)
        return (
            float(np.dot(self.axis, x)),
            float(np.dot(self.t1, x)),
            float(np.dot(self.t2, x)),
        )

    def _inside_grid_box(self, pt) -> bool:
        x = np.asarray(pt, dtype=float)
        return bool(np.all(x >= self.bounds_min) and np.all(x < self.bounds_max))

    def isInterior(self, pt, skinDepth):
        if not self._inside_grid_box(pt):
            return -1
        s, t1, t2 = self._local(pt)
        if abs(t1) > self.half_width_1 or abs(t2) > self.half_width_2:
            return -1
        signed = s - self.offset
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
        return self.axis.copy() if self.negative_side else -self.axis.copy()

    def getSurfacePosition(self, pt):
        s, _, _ = self._local(pt)
        return -(s - self.offset) * self.axis

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
pfw["mWallTime"] = os.environ.get("WEAK_TRACE_WALLTIME", "00:15:00")
pfw["mSubmitJobs"] = False

# Domain and discretization -------------------------------------------------------------
refine = int(os.environ.get("WEAK_TRACE_REFINE", "2"))
cpp = int(os.environ.get("WEAK_TRACE_CPP", "30"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", "1"))
surface_refinement_level = int(os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_LEVEL", "4"))
surface_refinement_method = os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_METHOD", "simplex")
surface_refinement_fallback = os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_FALLBACK", "cartesian")
surface_refinement_plane_tolerance = float(os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_PLANE_TOLERANCE", "0.05"))
particle_reference_geometry = os.environ.get("WEAK_TRACE_PARTICLE_REFERENCE_GEOMETRY", "fromRefinement")
particle_domain_output = os.environ.get("WEAK_TRACE_PARTICLE_DOMAIN_OUTPUT", "native")

pfw["planeStrain"] = 0
pfw["periodic"] = [False, False, False]
pfw["xpar"] = int(os.environ.get("WEAK_TRACE_XPAR", str(refine)))
pfw["ypar"] = int(os.environ.get("WEAK_TRACE_YPAR", str(refine)))
pfw["zpar"] = int(os.environ.get("WEAK_TRACE_ZPAR", str(refine)))
pfw["nI"] = int(os.environ.get("WEAK_TRACE_NI", str(cpp)))
pfw["nJ"] = int(os.environ.get("WEAK_TRACE_NJ", str(cpp)))
pfw["nK"] = int(os.environ.get("WEAK_TRACE_NK", str(cpp)))
if pfw["nI"] <= 2 or pfw["nJ"] <= 2 or pfw["nK"] <= 2:
    raise ValueError("WEAK_TRACE_NI/NJ/NK must leave at least one interior cell")
if pfw["nI"] % pfw["xpar"] != 0 or pfw["nJ"] % pfw["ypar"] != 0 or pfw["nK"] % pfw["zpar"] != 0:
    raise ValueError("WEAK_TRACE_NI/NJ/NK must be divisible by WEAK_TRACE_XPAR/YPAR/ZPAR")
pfw["ppc"] = ppc
pfw["surfaceRefinementLevel"] = surface_refinement_level
pfw["surfaceRefinementMethod"] = surface_refinement_method
pfw["surfaceRefinementFallback"] = surface_refinement_fallback
pfw["surfaceRefinementPlaneTolerance"] = surface_refinement_plane_tolerance
pfw["particleReferenceGeometry"] = particle_reference_geometry
pfw["particleDomainOutput"] = particle_domain_output

axis = unit([1.0, 1.0, 1.0])
transverse_1 = unit([1.0, -1.0, 0.0])
transverse_2 = unit(np.cross(axis, transverse_1))
bar_half_width = float(os.environ.get("WEAK_TRACE_BAR_HALF_WIDTH", "0.18"))
interface_offset = float(os.environ.get("WEAK_TRACE_INTERFACE_OFFSET", "0.018"))
side_buffer_cells = int(os.environ.get("WEAK_TRACE_SIDE_BUFFER_CELLS", "4"))
if side_buffer_cells < 0:
    raise ValueError("WEAK_TRACE_SIDE_BUFFER_CELLS must be nonnegative")

requested_half_extent = float(os.environ.get("WEAK_TRACE_BOX_HALF_EXTENT", "0.75"))
interior_cells = min(pfw["nI"] - 2, pfw["nJ"] - 2, pfw["nK"] - 2)
if interior_cells <= 2 * side_buffer_cells:
    raise ValueError("grid is too coarse for the requested side-buffer cells")
transverse_coordinate_extent = float(max(
    bar_half_width * (abs(transverse_1[i]) + abs(transverse_2[i])) for i in range(3)
))
interface_coordinate_extent = float(abs(interface_offset) * max(abs(v) for v in axis) + transverse_coordinate_extent)
minimum_half_extent = float(interface_coordinate_extent / (1.0 - 2.0 * side_buffer_cells / float(interior_cells)))
box_half_extent = float(max(requested_half_extent, 1.05 * minimum_half_extent))
pfw["xmin"] = -box_half_extent
pfw["xmax"] = box_half_extent
pfw["ymin"] = -box_half_extent
pfw["ymax"] = box_half_extent
pfw["zmin"] = -box_half_extent
pfw["zmax"] = box_half_extent
cell_size = float(2.0 * box_half_extent / float(interior_cells))
interface_boundary_clearance = float(box_half_extent - interface_coordinate_extent)
interface_boundary_clearance_cells = float(interface_boundary_clearance / cell_size)

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

# Equal stretch on all box faces moves the two occupied corner tabs apart along
# the long diagonal axis.  The sample sides near the central interface are free.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 2, 2]
pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [
    [0.0, 1.0, 1.0, 1.0],
    [stop_time, 1.0 + final_strain, 1.0 + final_strain, 1.0 + final_strain],
]

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
left_name = "weakTrace3DRotatedCompliant"
right_name = "weakTrace3DRotatedStiff"
density = float(os.environ.get("WEAK_TRACE_DENSITY", "1.0"))
poisson = float(os.environ.get("WEAK_TRACE_POISSON", "0.25"))
left_E = float(os.environ.get("WEAK_TRACE_LEFT_E", "2.0"))
right_E = float(os.environ.get("WEAK_TRACE_RIGHT_E", "8.0"))
pfw["materials"] = [left_name, right_name]
pfw["materialPropertyString"] = elastic_xml(left_name, density, left_E, poisson) + "\n" + elastic_xml(right_name, density, right_E, poisson)

# Geometry -------------------------------------------------------------------------------
trace_groups = VARIANT == "traceContactGroups"
refinement_uses_surface_flags = surface_refinement_method.lower().replace("-conforming", "") == "simplex"
flag_interface = trace_groups or refinement_uses_surface_flags
grid_bounds_min = [pfw["xmin"], pfw["ymin"], pfw["zmin"]]
grid_bounds_max = [pfw["xmax"], pfw["ymax"], pfw["zmax"]]

pfw["objects"] = [
    DiagonalPrismHalfSpace(
        "negative_diagonal_prism",
        axis,
        transverse_1,
        transverse_2,
        bar_half_width,
        bar_half_width,
        interface_offset,
        negative_side=True,
        flag_interface=flag_interface,
        bounds_min=grid_bounds_min,
        bounds_max=grid_bounds_max,
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        particle_type=particle_type,
    ),
    DiagonalPrismHalfSpace(
        "positive_diagonal_prism",
        axis,
        transverse_1,
        transverse_2,
        bar_half_width,
        bar_half_width,
        interface_offset,
        negative_side=False,
        flag_interface=flag_interface,
        bounds_min=grid_bounds_min,
        bounds_max=grid_bounds_max,
        vel=[0.0, 0.0, 0.0],
        mat=1,
        group=1 if trace_groups else 0,
        particle_type=particle_type,
    ),
]

# Tracers sample both materials near the interface and away from it along the bar.
tracer_points = []
for s in [interface_offset - 0.08, interface_offset - 0.025, interface_offset + 0.025, interface_offset + 0.08]:
    for t1, t2 in [(0.0, 0.0), (0.45 * bar_half_width, 0.0), (0.0, -0.45 * bar_half_width)]:
        p = s * axis + t1 * transverse_1 + t2 * transverse_2
        tracer_points.append((float(p[0]), float(p[1]), float(p[2])))
for s in [-0.55 * box_half_extent * math.sqrt(3.0), 0.55 * box_half_extent * math.sqrt(3.0)]:
    p = s * axis
    tracer_points.append((float(p[0]), float(p[1]), float(p[2])))

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
    output_prefix="weakTraceRotatedBar3DTracer",
)

pfw_expected = {
    "case_id": CASE_ID,
    "variant": VARIANT,
    "variant_label": VARIANT_LABEL,
    "stop_time": stop_time,
    "final_strain": final_strain,
    "plane_strain": pfw["planeStrain"],
    "bar_axis": [float(v) for v in axis],
    "transverse_1": [float(v) for v in transverse_1],
    "transverse_2": [float(v) for v in transverse_2],
    "interface_normal": [float(v) for v in axis],
    "interface_offset": interface_offset,
    "bar_half_width": bar_half_width,
    "box_half_extent": box_half_extent,
    "requested_box_half_extent": requested_half_extent,
    "side_buffer_cells_requested": side_buffer_cells,
    "interface_boundary_clearance": interface_boundary_clearance,
    "interface_boundary_clearance_cells": interface_boundary_clearance_cells,
    "cell_size": cell_size,
    "grid_bounds_min": grid_bounds_min,
    "grid_bounds_max": grid_bounds_max,
    "particle_type": particle_type,
    "update_method": pfw["updateMethod"],
    "update_order": pfw["updateOrder"],
    "surface_refinement_level": surface_refinement_level,
    "surface_refinement_method": surface_refinement_method,
    "surface_refinement_fallback": surface_refinement_fallback,
    "surface_refinement_plane_tolerance": surface_refinement_plane_tolerance,
    "particle_reference_geometry": particle_reference_geometry,
    "particle_domain_output": particle_domain_output,
    "surface_flags_for_refinement": flag_interface,
    "xpar": pfw["xpar"],
    "ypar": pfw["ypar"],
    "zpar": pfw["zpar"],
    "nI": pfw["nI"],
    "nJ": pfw["nJ"],
    "nK": pfw["nK"],
}
