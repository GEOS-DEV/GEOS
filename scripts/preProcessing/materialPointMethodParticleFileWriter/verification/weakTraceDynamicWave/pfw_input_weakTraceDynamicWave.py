# -*- coding: utf-8 -*-
"""Dynamic wave transmission across an explicit weak-trace interface."""

from __future__ import annotations

import math
import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db


pfw = {}

CASE_ID = "weakTraceDynamicWave"
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

# Plane-strain quasi-1D mesh.  nK=3 is required so plane strain has one
# interior grid-cell layer through thickness.
refine = int(os.environ.get("WEAK_TRACE_REFINE", "1"))
cpp = int(os.environ.get("WEAK_TRACE_CPP", "48"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
default_particle_type = "1" if VARIANT == "primalEnrichment" else "2"
particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", default_particle_type))
left_length = float(os.environ.get("WEAK_TRACE_LEFT_LENGTH", "1.0"))
right_length = float(os.environ.get("WEAK_TRACE_RIGHT_LENGTH", "1.5"))
bar_length = left_length + right_length
interface_x = -right_length
bar_height = 0.12
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
cells_x_blocks = max(1, int(math.ceil(refine * bar_length)))
cells_y_blocks = max(1, refine)
if VARIANT == "primalEnrichment":
    # The primal enrichment prototype is serial-only. Keep the same cell count
    # used by the other variants, but do not split the background grid over MPI
    # partitions.
    pfw["xpar"] = 1
    pfw["ypar"] = 1
else:
    pfw["xpar"] = cells_x_blocks
    pfw["ypar"] = cells_y_blocks
pfw["zpar"] = 1
pfw["nI"] = cells_x_blocks * cpp
pfw["nJ"] = cells_y_blocks * max(8, cpp // 6)
pfw["nK"] = 3
pfw["ppc"] = ppc

pfw["xmin"] = -bar_length
pfw["xmax"] = 0.0
pfw["ymin"] = -0.5 * bar_height
pfw["ymax"] = 0.5 * bar_height
dx = (pfw["xmax"] - pfw["xmin"]) / float(max(1, pfw["nI"] - 2))
dy = (pfw["ymax"] - pfw["ymin"]) / float(max(1, pfw["nJ"] - 2))
thickness = min(dx, dy)
pfw["zmin"] = -0.5 * thickness
pfw["zmax"] = 0.5 * thickness

stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "1.6"))
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
pfw["profileVariables"] = "{ velocityX, stressXX }"
pfw["frictionCoefficient"] = 0.0
pfw["contactNormalType"] = "Aligned"
pfw["useSurfacePositionForContact"] = 1
pfw["damageFieldPartitioning"] = 0
pfw["weakDiscontinuityInsertion"] = "off"
pfw["weakInterfaceTraceGeometry"] = "surface"
pfw["cpdiDomainScaling"] = 0
default_update_method = "PIC" if VARIANT == "primalEnrichment" else "FLIP"
pfw["updateMethod"] = os.environ.get("WEAK_TRACE_UPDATE_METHOD", default_update_method)
pfw["updateOrder"] = int(os.environ.get("WEAK_TRACE_UPDATE_ORDER", "1"))
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
    "gridWeakPrimalCrossMass",
    "gridWeakPrimalRawMass",
    "gridWeakPrimalRawMomentum",
    "gridWeakPrimalRawInternalForce",
    "gridWeakPrimalMean",
    "gridWeakPrimalMass",
    "gridWeakPrimalMassRatio",
    "gridWeakPrimalVelocity",
    "gridWeakPrimalAcceleration",
    "gridWeakPrimalActive",
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

pulse_time_center = float(os.environ.get("WEAK_TRACE_PULSE_TIME_CENTER", "0.20"))
pulse_time_width = float(os.environ.get("WEAK_TRACE_PULSE_TIME_WIDTH", "0.05"))
pulse_velocity = float(os.environ.get("WEAK_TRACE_PULSE_VELOCITY", "0.02"))
pulse_table_samples = int(os.environ.get("WEAK_TRACE_PULSE_TABLE_SAMPLES", "161"))


def boundary_velocity(time):
    return pulse_velocity * math.exp(-((time - pulse_time_center) / pulse_time_width) ** 2)


def make_boundary_velocity_ftable():
    sample_count = max(2, pulse_table_samples)
    times = [stop_time * i / float(sample_count - 1) for i in range(sample_count)]
    log_f = 0.0
    rows = [[times[0], 1.0, 1.0, 1.0]]
    previous_time = times[0]
    previous_rate = boundary_velocity(previous_time) / pfw["xmin"]
    for time in times[1:]:
        rate = boundary_velocity(time) / pfw["xmin"]
        log_f += 0.5 * (previous_rate + rate) * (time - previous_time)
        rows.append([time, math.exp(log_f), 1.0, 1.0])
        previous_time = time
        previous_rate = rate
    return rows


f_table_rows = make_boundary_velocity_ftable()

# Constrain transverse strain and drive the left boundary with an F-table that
# approximates a prescribed Gaussian boundary velocity.  The right boundary is
# placed at x=0; the same F-table therefore prescribes zero normal velocity there.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Linear"
pfw["fTable"] = f_table_rows

pfw["enableWeakInterfacePrimalEnrichment"] = 0
if VARIANT == "primalEnrichment":
    pfw["enableContact"] = 0
    pfw["enableWeakInterfaceTraceProjection"] = 0
    pfw["enableWeakInterfacePrimalEnrichment"] = 1
    pfw["weakInterfacePrimalGeometryType"] = 1
    pfw["weakInterfacePrimalPointX"] = float(os.environ.get("PRIMAL_WEAK_POINT_X", str(interface_x)))
    pfw["weakInterfacePrimalPointY"] = float(os.environ.get("PRIMAL_WEAK_POINT_Y", "0.0"))
    pfw["weakInterfacePrimalNormalX"] = float(os.environ.get("PRIMAL_WEAK_NORMAL_X", "1.0"))
    pfw["weakInterfacePrimalNormalY"] = float(os.environ.get("PRIMAL_WEAK_NORMAL_Y", "0.0"))
    pfw["weakInterfacePrimalMinMassRatio"] = float(os.environ.get("PRIMAL_WEAK_MIN_MASS_RATIO", "1.0e-8"))
    pfw["weakInterfacePrimalMassCFL"] = float(os.environ.get("PRIMAL_WEAK_MASS_CFL", "0.35"))
elif VARIANT == "traceContactGroups":
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

left_name = "weakTraceWaveLeftElastic"
right_name = "weakTraceWaveRightElastic"
left_density = float(os.environ.get("WEAK_TRACE_LEFT_DENSITY", "1.0"))
right_density = float(os.environ.get("WEAK_TRACE_RIGHT_DENSITY", "1.0"))
left_E = float(os.environ.get("WEAK_TRACE_LEFT_E", "1.0"))
right_E = float(os.environ.get("WEAK_TRACE_RIGHT_E", "4.0"))
poisson = float(os.environ.get("WEAK_TRACE_POISSON", "0.25"))
pfw["materials"] = [left_name, right_name]
pfw["materialPropertyString"] = elastic_xml(left_name, left_density, left_E, poisson) + "\n" + elastic_xml(right_name, right_density, right_E, poisson)
if VARIANT == "falseElasticCZ":
    cz = material_db.weakTraceFalseElasticCohesiveZone.copy()
    cz_stiffness = os.environ.get("WEAK_TRACE_CZ_STIFFNESS", "")
    cz_normal = float(os.environ.get("WEAK_TRACE_CZ_NORMAL_FORCE_CONSTANT", cz_stiffness or str(cz["normalForceConstant"])))
    cz_shear = float(os.environ.get("WEAK_TRACE_CZ_SHEAR_FORCE_CONSTANT", cz_stiffness or str(cz["shearForceConstant"])))
    cz["normalForceConstant"] = cz_normal
    cz["shearForceConstant"] = cz_shear
    cz["materialString"] = material_db.generateMaterialString(cz)
    pfw["materialPropertyString"] += "\n" + cz["materialString"]
    pfw["cohesiveZoneRegions"] = f"""
<CohesiveZoneRegion
    name="weakTraceDynamicWaveElasticTie"
    constitutiveModel="{cz['name']}"
    tag="0"/>"""
    pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="weakTraceDynamicWaveElasticTie"
    startTime="0.0"
    regionNames="{ weakTraceDynamicWaveElasticTie }"
    czVolumeNormalization="1"/>"""

def interface_group(side: int) -> int:
    return side if VARIANT == "traceContactGroups" else 0


def mark_interface(name: str, obj):
    if VARIANT == "traceContactGroups":
        return geom.surfaceFlagWrapper(name + "_weak", obj, int(geom.SurfaceFlag.WeakDiscontinuity))
    if VARIANT == "falseElasticCZ":
        return geom.czTagWrapper(name + "_cz_tag", geom.surfaceFlagWrapper(name + "_cz", obj, int(geom.SurfaceFlag.Cohesive)), 0)
    return obj


left = geom.box(
    "left_bar",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [-right_length, pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=interface_group(0),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, False, True, False],
)
right = geom.box(
    "right_bar",
    [-right_length, pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=1,
    group=interface_group(1),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[True, False, False, False],
)
pfw["objects"] = [mark_interface("left_interface", left), mark_interface("right_interface", right)]

tracer_x = [interface_x + offset for offset in [-0.75, -0.25, 0.25, 0.75]]
tracer_points = [(x, 0.0, 0.0) for x in tracer_x]
tracers.set_tracers(
    pfw,
    tracer_points,
    ["velocityX", "stressXX", "stressYY", "stressXY", "materialType"],
    write_interval=stop_time / 500.0,
    output_prefix="weakTraceDynamicWaveTracer",
)

plane_strain_modulus_left = left_E * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
plane_strain_modulus_right = right_E * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
pfw_expected = {
    "case_id": CASE_ID,
    "variant_label": VARIANT_LABEL,
    "xmin": pfw["xmin"],
    "xmax": pfw["xmax"],
    "interface_x": interface_x,
    "left_density": left_density,
    "right_density": right_density,
    "left_E": left_E,
    "right_E": right_E,
    "poisson": poisson,
    "left_length": left_length,
    "right_length": right_length,
    "bar_length": bar_length,
    "left_plane_strain_modulus": plane_strain_modulus_left,
    "right_plane_strain_modulus": plane_strain_modulus_right,
    "left_wave_speed": math.sqrt(plane_strain_modulus_left / left_density),
    "right_wave_speed": math.sqrt(plane_strain_modulus_right / right_density),
    "left_impedance": math.sqrt(left_density * plane_strain_modulus_left),
    "right_impedance": math.sqrt(right_density * plane_strain_modulus_right),
    "source": "left_boundary_gaussian_velocity",
    "source_x": pfw["xmin"],
    "stop_time": stop_time,
    "pulse_time_center": pulse_time_center,
    "pulse_time_width": pulse_time_width,
    "pulse_velocity": pulse_velocity,
    "pulse_table_samples": pulse_table_samples,
    "update_method": pfw["updateMethod"],
    "update_order": pfw["updateOrder"],
    "particle_type": particle_type,
    "primal_enrichment": pfw.get("enableWeakInterfacePrimalEnrichment", 0),
    "primal_geometry_type": pfw.get("weakInterfacePrimalGeometryType", 0),
    "primal_min_mass_ratio": pfw.get("weakInterfacePrimalMinMassRatio", None),
}
