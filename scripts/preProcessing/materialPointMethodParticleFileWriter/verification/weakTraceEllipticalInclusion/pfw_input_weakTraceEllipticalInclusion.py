# -*- coding: utf-8 -*-
"""Elliptical inclusion under remote uniaxial loading for weak trace tests."""

from __future__ import annotations

import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db


pfw = {}

CASE_ID = "weakTraceEllipticalInclusion"
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
cpp = int(os.environ.get("WEAK_TRACE_CPP", "56"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = ppc

pfw["xmin"] = -0.75
pfw["xmax"] = 0.75
pfw["ymin"] = -0.75
pfw["ymax"] = 0.75
cell = min(
    (pfw["xmax"] - pfw["xmin"]) / float(max(1, pfw["nI"] - 2)),
    (pfw["ymax"] - pfw["ymin"]) / float(max(1, pfw["nJ"] - 2)),
)
pfw["zmin"] = -0.5 * cell
pfw["zmax"] = 0.5 * cell

stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "10.0"))
final_strain = float(os.environ.get("WEAK_TRACE_FINAL_STRAIN", "0.004"))
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

pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 0, 0, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0 + final_strain, 1.0, 1.0]]

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
    pfw["useEvents"] = 1
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
    name="weakTraceEllipticalInclusionElasticTie"
    constitutiveModel="{cz['name']}"
    tag="0"/>"""
    pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="weakTraceEllipticalInclusionElasticTie"
    startTime="0.0"
    regionNames="{ weakTraceEllipticalInclusionElasticTie }"
    czVolumeNormalization="1"/>"""

axis_a = float(os.environ.get("WEAK_TRACE_ELLIPSE_A", "0.30"))
axis_b = float(os.environ.get("WEAK_TRACE_ELLIPSE_B", "0.13"))
particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", "2"))


def interface_group(side: int) -> int:
    return side if VARIANT == "traceContactGroups" else 0


def ellipse(name: str, mat: int, group: int):
    return geom.ellipsoid(
        name,
        [0.0, 0.0, 0.0],
        axis_a,
        axis_b,
        10.0 * cell,
        vel=[0.0, 0.0, 0.0],
        mat=mat,
        group=group,
        particleType=particle_type,
    )


def mark_interface(name: str, obj):
    if VARIANT == "traceContactGroups":
        return geom.surfaceFlagWrapper(name + "_weak", obj, int(geom.SurfaceFlag.WeakDiscontinuity))
    if VARIANT == "falseElasticCZ":
        return geom.czTagWrapper(name + "_cz_tag", geom.surfaceFlagWrapper(name + "_cz", obj, int(geom.SurfaceFlag.Cohesive)), 0)
    return obj


matrix_box = geom.box(
    "matrix_box",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=interface_group(0),
    dim=2,
    flaggedSurfaces=[False, False, False, False],
)
inclusion_cut = ellipse("elliptical_inclusion_cut", 1, interface_group(1))
matrix = geom.difference("matrix_with_elliptical_hole", matrix_box, inclusion_cut)
inclusion = ellipse("elliptical_inclusion", 1, interface_group(1))
pfw["objects"] = [mark_interface("matrix_interface", matrix), mark_interface("inclusion_interface", inclusion)]

tracer_points = []
for x in [-0.45, -axis_a * 1.05, -axis_a * 0.95, 0.0, axis_a * 0.95, axis_a * 1.05, 0.45]:
    tracer_points.append((x, 0.0, 0.0))
for y in [-0.35, -axis_b * 1.05, -axis_b * 0.95, axis_b * 0.95, axis_b * 1.05, 0.35]:
    tracer_points.append((0.0, y, 0.0))
tracers.set_tracers(
    pfw,
    tracer_points,
    ["velocityX", "velocityY", "stressXX", "stressYY", "stressXY", "materialType"],
    write_interval=stop_time / 400.0,
    output_prefix="weakTraceEllipticalInclusionTracer",
)

pfw_expected = {
    "case_id": CASE_ID,
    "variant_label": VARIANT_LABEL,
    "ellipse_a": axis_a,
    "ellipse_b": axis_b,
    "final_remote_x_strain": final_strain,
    "matrix_E": matrix_E,
    "inclusion_E": inclusion_E,
    "poisson": poisson,
}
