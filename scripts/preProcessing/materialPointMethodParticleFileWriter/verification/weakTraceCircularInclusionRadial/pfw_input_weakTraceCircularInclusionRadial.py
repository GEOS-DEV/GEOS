# -*- coding: utf-8 -*-
"""Circular inclusion under biaxial radial-type loading for weak trace tests."""

from __future__ import annotations

import math
import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db


pfw = {}

CASE_ID = "weakTraceCircularInclusionRadial"
VARIANT = os.environ.get("WEAK_TRACE_VARIANT", "singleField")
VARIANT_LABEL = os.environ.get("VV_VARIANT_LABEL", VARIANT)
CASE_NAME = os.environ.get("VV_CASE_NAME", f"{CASE_ID}_{VARIANT}")


def elastic_xml(name: str, density: float, young: float, poisson: float) -> str:
    return f"""
<ElasticIsotropic
    name="{name}"
    defaultDensity="{density}"
    defaultYoungModulus="{young}"
    defaultPoissonRatio="{poisson}"/>"""


pfw["caseName"] = CASE_NAME
pfw["runDebug"] = True
pfw["mBatch"] = True
pfw["mWallTime"] = os.environ.get("WEAK_TRACE_WALLTIME", "00:10:00")
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("WEAK_TRACE_REFINE", "1"))
cpp = int(os.environ.get("WEAK_TRACE_CPP", "64"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
inclusion_radius = float(os.environ.get("WEAK_TRACE_INCLUSION_RADIUS", "0.22"))
outer_radius = float(os.environ.get("WEAK_TRACE_OUTER_RADIUS", "0.72"))
domain_half_width = float(os.environ.get("WEAK_TRACE_DOMAIN_HALF_WIDTH", str(1.15 * outer_radius)))
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = ppc

pfw["xmin"] = -domain_half_width
pfw["xmax"] = domain_half_width
pfw["ymin"] = -domain_half_width
pfw["ymax"] = domain_half_width
cell = min(
    (pfw["xmax"] - pfw["xmin"]) / float(max(1, pfw["nI"] - 2)),
    (pfw["ymax"] - pfw["ymin"]) / float(max(1, pfw["nJ"] - 2)),
)
pfw["zmin"] = -0.5 * cell
pfw["zmax"] = 0.5 * cell

stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "8.0"))
external_pressure = float(os.environ.get("WEAK_TRACE_EXTERNAL_PRESSURE", "0.004"))
pressure_ramp_time = float(os.environ.get("WEAK_TRACE_PRESSURE_RAMP_TIME", str(stop_time)))
confining_half_width = float(os.environ.get("WEAK_TRACE_CONFINING_BOX_HALF_WIDTH", str(0.40 * outer_radius)))
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 80.0
pfw["restartInterval"] = 10.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = float(os.environ.get("WEAK_TRACE_CFL", "0.20"))
pfw["initialDt"] = 1.0e-16
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 400.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 400.0
pfw["profileHistory"] = 1
pfw["profileDirection"] = "x"
pfw["profileNumSlices"] = 0
pfw["profileWriteInterval"] = pfw["plotInterval"]
pfw["profileCycleInterval"] = 0
pfw["profileVariables"] = "{ velocityX, velocityY, stressXX, stressYY, stressXY }"
pfw["frictionCoefficient"] = 0.0
pfw["contactNormalType"] = "Aligned"
pfw["useSurfacePositionForContact"] = 1
pfw["damageFieldPartitioning"] = 0
pfw["weakDiscontinuityInsertion"] = "off"
pfw["weakInterfaceTraceGeometry"] = "surface"
pfw["cpdiDomainScaling"] = 0
pfw["updateMethod"] = os.environ.get("WEAK_TRACE_UPDATE_METHOD", "PIC")
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = [
    "gridMass",
    "gridVelocity",
    "gridInternalForce",
    "gridActive",
    "gridSurfaceFieldMass",
    "gridSurfaceNormal",
    "gridExplicitSurfaceNormal",
    "gridPrincipalExplicitSurfaceNormal",
    "gridSurfacePosition",
    "gridContactForce",
    "gridCohesiveNode",
    "gridCohesiveFieldFlag",
    "gridCohesiveArea",
    "gridCohesiveForce",
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
pfw["plottableFields"] = ["particleStress", "particleVelocity", "particleDisplacement"]
pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "SurfaceFlag",
    "RVector",
    "SurfaceNormal",
    "SurfacePosition",
]
if VARIANT == "falseElasticCZ":
    pfw["particleFileFields"].insert(pfw["particleFileFields"].index("RVector"), "CZTag")

# The first radial-inclusion runs use the same virtual-fluid pressure machinery
# as the borehole example.  The physical specimen is a circular disk; the
# confining-pressure event supplies a background stress in the surrounding
# virtual region so the load is applied as a traction instead of an affine square
# boundary stretch.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [1, 1, 1, 1, 1, 1]
pfw["prescribedBoundaryFTable"] = 0
pfw["useEvents"] = 1

if VARIANT == "traceContactGroups":
    pfw["enableContact"] = 1
    pfw["enableWeakInterfaceTraceProjection"] = 1
    pfw["explicitSurfaceNormalInfluence"] = float(os.environ.get("WEAK_TRACE_EXPLICIT_NORMAL_INFLUENCE", "25.0"))
    pfw["weakInterfaceTraceProjectionIterations"] = int(os.environ.get("WEAK_TRACE_PROJECTION_ITERATIONS", "8"))
    pfw["weakInterfaceTraceProjectionScale"] = float(os.environ.get("WEAK_TRACE_PROJECTION_SCALE", "1.0"))
    pfw["weakInterfaceTraceGapStabilization"] = float(os.environ.get("WEAK_TRACE_GAP_STABILIZATION", "0.01"))
    pfw["weakInterfaceTraceMaxGapCorrectionVelocity"] = float(os.environ.get("WEAK_TRACE_MAX_GAP_CORRECTION_VELOCITY", "1.0e-5"))
    pfw["weakInterfaceTraceMinWeight"] = float(os.environ.get("WEAK_TRACE_MIN_WEIGHT", "1.0e-12"))
    pfw["weakInterfaceTraceSuppressNodalContact"] = 1
    pfw["weakInterfaceTracePairs"] = [[0, 1]]
elif VARIANT == "falseElasticCZ":
    pfw["enableContact"] = 1
    pfw["enableWeakInterfaceTraceProjection"] = 0
    pfw["explicitSurfaceNormalInfluence"] = float(os.environ.get("WEAK_TRACE_EXPLICIT_NORMAL_INFLUENCE", "25.0"))
    pfw["damageFieldPartitioning"] = 1
else:
    pfw["enableContact"] = 0
    pfw["enableWeakInterfaceTraceProjection"] = 0

matrix_name = "weakTraceMatrixElastic"
inclusion_name = "weakTraceInclusionElastic"
matrix_E = 1.0
inclusion_E = 4.0
poisson = 0.25
density = 1.0
pfw["materials"] = [matrix_name, inclusion_name]
pfw["materialPropertyString"] = elastic_xml(matrix_name, density, matrix_E, poisson) + "\n" + elastic_xml(inclusion_name, density, inclusion_E, poisson)
if VARIANT == "falseElasticCZ":
    cz = material_db.weakTraceFalseElasticCohesiveZone.copy()
    pfw["materialPropertyString"] += "\n" + cz["materialString"]
    pfw["cohesiveZoneRegions"] = f"""
<CohesiveZoneRegion
    name="weakTraceCircularInclusionElasticTie"
    constitutiveModel="{cz['name']}"
    tag="0"/>"""

event_blocks = [
    f"""
<ConfiningPressure
    name="weakTraceCircularInclusionOuterPressure"
    startTime="0.0"
    endTime="{pressure_ramp_time}"
    confiningPressureBoxMin="{{{-confining_half_width},{-confining_half_width},{-2.0 * cell}}}"
    confiningPressureBoxMax="{{{confining_half_width},{confining_half_width},{2.0 * cell}}}"
    startPressure="0.0"
    endPressure="{external_pressure}"
    interpType="1"/>"""
]
if VARIANT == "falseElasticCZ":
    event_blocks.append(
        """
<ReferenceCohesiveZones
    name="weakTraceCircularInclusionElasticTie"
    startTime="0.0"
    regionNames="{ weakTraceCircularInclusionElasticTie }"
    czVolumeNormalization="1"/>"""
    )
pfw["mpmEventsString"] = "\n".join(event_blocks)

particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", "2"))
geometry_z_pad = float(os.environ.get("WEAK_TRACE_GEOMETRY_Z_PAD", str(2.0 * cell)))
cylinder_zmin = pfw["zmin"] - geometry_z_pad
cylinder_zmax = pfw["zmax"] + geometry_z_pad


class radialInterfaceWrapper(geom.BaseWrapper):
    def __init__(self, name: str, sub_object, radius: float, normal_sign: float):
        super().__init__(name, sub_object)
        self.radius = float(radius)
        self.normal_sign = float(normal_sign)

    def _radial_unit(self, pt):
        x = float(pt[0])
        y = float(pt[1])
        r = math.hypot(x, y)
        if r <= 1.0e-14:
            return 1.0, 0.0, 0.0
        return x / r, y / r, r

    def getSurfaceNormal(self, pt):
        ux, uy, _ = self._radial_unit(pt)
        return [self.normal_sign * ux, self.normal_sign * uy, 0.0]

    def getSurfacePosition(self, pt):
        ux, uy, r = self._radial_unit(pt)
        distance = self.radius - r
        return [distance * ux, distance * uy, 0.0]


def interface_group(side: int) -> int:
    return side if VARIANT == "traceContactGroups" else 0


def disk(name: str, mat: int, group: int):
    return geom.cylinder(
        name,
        [0.0, 0.0, cylinder_zmin],
        [0.0, 0.0, cylinder_zmax],
        inclusion_radius,
        vel=[0.0, 0.0, 0.0],
        mat=mat,
        group=group,
        particleType=particle_type,
        flaggedSurfaces=[False, False, True, False],
    )


def annulus(name: str, mat: int, group: int):
    return geom.cylinder(
        name,
        [0.0, 0.0, cylinder_zmin],
        [0.0, 0.0, cylinder_zmax],
        outer_radius,
        ri=inclusion_radius,
        vel=[0.0, 0.0, 0.0],
        mat=mat,
        group=group,
        particleType=particle_type,
        flaggedSurfaces=[False, False, False, True],
    )


def mark_interface(name: str, obj):
    if VARIANT == "traceContactGroups":
        return geom.surfaceFlagWrapper(name + "_weak", obj, int(geom.SurfaceFlag.WeakDiscontinuity))
    if VARIANT == "falseElasticCZ":
        return geom.czTagWrapper(name + "_cz_tag", geom.surfaceFlagWrapper(name + "_cz", obj, int(geom.SurfaceFlag.Cohesive)), 0)
    return obj


matrix = radialInterfaceWrapper(
    "matrix_annulus_radial_interface",
    annulus("matrix_annulus", 0, interface_group(0)),
    inclusion_radius,
    -1.0,
)
inclusion = radialInterfaceWrapper(
    "circular_inclusion_radial_interface",
    disk("circular_inclusion", 1, interface_group(1)),
    inclusion_radius,
    1.0,
)
pfw["objects"] = [mark_interface("matrix_interface", matrix), mark_interface("inclusion_interface", inclusion)]

tracer_points = []
axis_coordinates = [
    -0.90 * outer_radius,
    -0.70 * outer_radius,
    -1.08 * inclusion_radius,
    -0.96 * inclusion_radius,
    -0.50 * inclusion_radius,
    0.50 * inclusion_radius,
    0.96 * inclusion_radius,
    1.08 * inclusion_radius,
    0.70 * outer_radius,
    0.90 * outer_radius,
]
for x in axis_coordinates:
    tracer_points.append((x, 0.0, 0.0))
for y in axis_coordinates:
    tracer_points.append((0.0, y, 0.0))
for r in [0.12, inclusion_radius * 0.95, inclusion_radius * 1.05, 0.40]:
    tracer_points.extend([(r, 0.0, 0.0), (0.0, r, 0.0), (-r, 0.0, 0.0), (0.0, -r, 0.0)])

surface_tracer_count = int(os.environ.get("WEAK_TRACE_SURFACE_TRACER_COUNT", "32"))
surface_tracer_offset = float(os.environ.get("WEAK_TRACE_SURFACE_TRACER_OFFSET", str(0.5 * cell)))


def add_tracer_ring(radius: float, count: int):
    if count <= 0 or radius <= 0.0:
        return
    for i in range(count):
        theta = 2.0 * math.pi * (float(i) + 0.5) / float(count)
        tracer_points.append((radius * math.cos(theta), radius * math.sin(theta), 0.0))


add_tracer_ring(max(0.0, inclusion_radius - surface_tracer_offset), surface_tracer_count)
add_tracer_ring(min(outer_radius, inclusion_radius + surface_tracer_offset), surface_tracer_count)
add_tracer_ring(max(0.0, outer_radius - surface_tracer_offset), surface_tracer_count)
tracers.set_tracers(
    pfw,
    tracer_points,
    ["velocityX", "velocityY", "stressXX", "stressYY", "stressXY", "materialType"],
    write_interval=stop_time / 400.0,
    output_prefix="weakTraceCircularInclusionRadialTracer",
)

pfw_expected = {
    "case_id": CASE_ID,
    "variant_label": VARIANT_LABEL,
    "inclusion_radius": inclusion_radius,
    "outer_radius": outer_radius,
    "external_pressure": external_pressure,
    "pressure_ramp_time": pressure_ramp_time,
    "confining_box_half_width": confining_half_width,
    "geometry_z_pad": geometry_z_pad,
    "surface_tracer_count": surface_tracer_count,
    "surface_tracer_offset": surface_tracer_offset,
    "matrix_E": matrix_E,
    "inclusion_E": inclusion_E,
    "poisson": poisson,
}
