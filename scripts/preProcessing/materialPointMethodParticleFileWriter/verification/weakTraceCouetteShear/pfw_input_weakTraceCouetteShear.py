# -*- coding: utf-8 -*-
"""Parallel-plate Couette shear for weak trace no-slip testing."""

from __future__ import annotations

import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db


pfw = {}

CASE_ID = "weakTraceCouetteShear"
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
cpp = int(os.environ.get("WEAK_TRACE_CPP", "40"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
pfw["planeStrain"] = 1
pfw["periodic"] = [True, False, False]
pfw["xpar"] = max(3, 3 * refine)
pfw["ypar"] = max(1, refine)
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = ppc

pfw["xmin"] = -0.5
pfw["xmax"] = 0.5
pfw["ymin"] = -0.25
pfw["ymax"] = 0.25
cell = (pfw["ymax"] - pfw["ymin"]) / float(max(1, pfw["nJ"] - 2))
pfw["zmin"] = -0.5 * cell
pfw["zmax"] = 0.5 * cell

stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "5.0"))
plate_velocity = float(os.environ.get("WEAK_TRACE_COUETTE_VELOCITY", "0.002"))
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

# Boundary face order is x-, x+, y-, y+, z-, z+.  With planeStrain enabled,
# the second transverse entry on y-faces controls the x velocity.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Linear"
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0, 1.0, 1.0]]
pfw["enablePrescribedBoundaryTransverseVelocities"] = [0, 0, 1, 1, 0, 0]
pfw["prescribedBoundaryTransverseVelocities"] = [
    [0.0, 0.0],
    [0.0, 0.0],
    [0.0, -plate_velocity],
    [0.0, plate_velocity],
    [0.0, 0.0],
    [0.0, 0.0],
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
    pfw["weakInterfaceTraceLinearSolverTolerance"] = float(os.environ.get("WEAK_TRACE_LINEAR_SOLVER_TOLERANCE", "1.0e-10"))
    pfw["weakInterfaceTraceLinearSolverRegularization"] = float(os.environ.get("WEAK_TRACE_LINEAR_SOLVER_REGULARIZATION", "1.0e-12"))
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

material_name = "weakTraceCouetteElastic"
young = 1.0
poisson = 0.25
density = 1.0
pfw["materials"] = [material_name]
pfw["materialPropertyString"] = elastic_xml(material_name, density, young, poisson)
if VARIANT == "falseElasticCZ":
    cz = material_db.weakTraceFalseElasticCohesiveZone.copy()
    pfw["materialPropertyString"] += "\n" + cz["materialString"]
    pfw["cohesiveZoneRegions"] = f"""
<CohesiveZoneRegion
    name="weakTraceCouetteElasticTie"
    constitutiveModel="{cz['name']}"
    tag="0"/>"""
    pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="weakTraceCouetteElasticTie"
    startTime="0.0"
    regionNames="{ weakTraceCouetteElasticTie }"
    czVolumeNormalization="1"/>"""

particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", "2"))


def interface_group(side: int) -> int:
    return side if VARIANT == "traceContactGroups" else 0


def mark_interface(name: str, obj):
    if VARIANT == "traceContactGroups":
        return geom.surfaceFlagWrapper(name + "_weak", obj, int(geom.SurfaceFlag.WeakDiscontinuity))
    if VARIANT == "falseElasticCZ":
        return geom.czTagWrapper(name + "_cz_tag", geom.surfaceFlagWrapper(name + "_cz", obj, int(geom.SurfaceFlag.Cohesive)), 0)
    return obj


lower = geom.box(
    "lower_layer",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], 0.0, pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=interface_group(0),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, False, False, True],
)
upper = geom.box(
    "upper_layer",
    [pfw["xmin"], 0.0, pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=interface_group(1),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, False, True, False],
)
pfw["objects"] = [mark_interface("lower_interface", lower), mark_interface("upper_interface", upper)]

tracer_y = [-0.20, -0.10, -0.01, 0.01, 0.10, 0.20]
tracer_points = [(0.0, y, 0.0) for y in tracer_y]
tracers.set_tracers(
    pfw,
    tracer_points,
    ["velocityX", "velocityY", "stressXY", "stressXX", "stressYY", "materialType"],
    write_interval=stop_time / 400.0,
    output_prefix="weakTraceCouetteShearTracer",
)

pfw_expected = {
    "case_id": CASE_ID,
    "variant_label": VARIANT_LABEL,
    "plate_velocity_x": plate_velocity,
    "gap": pfw["ymax"] - pfw["ymin"],
    "young": young,
    "poisson": poisson,
}
