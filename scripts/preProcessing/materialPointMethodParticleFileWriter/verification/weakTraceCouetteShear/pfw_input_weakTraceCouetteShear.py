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


def liquid_xml(name: str, density: float, bulk: float, viscosity: float) -> str:
    return f"""
<Liquid
    name="{name}"
    defaultDensity="{density}"
    defaultBulkModulus="{bulk}"
    defaultViscosity="{viscosity}"/>"""


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


pfw["caseName"] = CASE_NAME
pfw["runDebug"] = env_bool("WEAK_TRACE_RUN_DEBUG", True)
pfw["mPartition"] = os.environ.get("WEAK_TRACE_PARTITION", "pdebug" if pfw["runDebug"] else "pbatch")
pfw["mBatch"] = True
pfw["mWallTime"] = os.environ.get("WEAK_TRACE_WALLTIME", "00:10:00")
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("WEAK_TRACE_REFINE", "1"))
cpp = int(os.environ.get("WEAK_TRACE_CPP", "40"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))
surface_refinement_level = int(os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_LEVEL", "1"))
pfw["planeStrain"] = 1
pfw["periodic"] = [True, False, False]
pfw["xpar"] = int(os.environ.get("WEAK_TRACE_XPAR", str(max(3, 3 * refine))))
pfw["ypar"] = int(os.environ.get("WEAK_TRACE_YPAR", str(max(3, 3 * refine))))
pfw["zpar"] = 1
pfw["nI"] = int(os.environ.get("WEAK_TRACE_NI", str(pfw["xpar"] * cpp)))
pfw["nJ"] = int(os.environ.get("WEAK_TRACE_NJ", str(pfw["ypar"] * cpp)))
pfw["nK"] = 3
if pfw["nI"] % pfw["xpar"] != 0 or pfw["nJ"] % pfw["ypar"] != 0:
    raise ValueError("WEAK_TRACE_NI/NJ must be divisible by WEAK_TRACE_XPAR/YPAR")
pfw["ppc"] = ppc
pfw["surfaceRefinementLevel"] = surface_refinement_level

domain_width = float(os.environ.get("WEAK_TRACE_DOMAIN_WIDTH", "1.0"))
domain_height = float(os.environ.get("WEAK_TRACE_DOMAIN_HEIGHT", "0.5"))
platen_thickness = float(os.environ.get("WEAK_TRACE_PLATEN_THICKNESS", "0.05"))
if 2.0 * platen_thickness >= domain_height:
    raise ValueError("WEAK_TRACE_PLATEN_THICKNESS must leave a positive fluid gap")

pfw["xmin"] = -0.5 * domain_width
pfw["xmax"] = 0.5 * domain_width
pfw["ymin"] = -0.5 * domain_height
pfw["ymax"] = 0.5 * domain_height
cell = (pfw["ymax"] - pfw["ymin"]) / float(max(1, pfw["nJ"] - 2))
pfw["zmin"] = -0.5 * cell
pfw["zmax"] = 0.5 * cell
fluid_ymin = pfw["ymin"] + platen_thickness
fluid_ymax = pfw["ymax"] - platen_thickness
fluid_gap = fluid_ymax - fluid_ymin

stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "5.0"))
top_plate_velocity = float(os.environ.get("WEAK_TRACE_COUETTE_VELOCITY", "0.002"))
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
pfw["profileDirection"] = "y"
pfw["profileNumSlices"] = 0
pfw["profileWriteInterval"] = pfw["plotInterval"]
pfw["profileCycleInterval"] = 0
pfw["profileVariables"] = "{ velocityX, velocityY, stressXY, stressXX, stressYY }"
pfw["frictionCoefficient"] = 0.0
pfw["contactNormalType"] = "Aligned"
pfw["useSurfacePositionForContact"] = 1
pfw["damageFieldPartitioning"] = 0
pfw["weakDiscontinuityInsertion"] = "off"
pfw["weakInterfaceTraceGeometry"] = "surface"
pfw["cpdiDomainScaling"] = 0
pfw["updateMethod"] = os.environ.get("WEAK_TRACE_UPDATE_METHOD", "PIC")
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
if VARIANT == "falseElasticCZ":
    pfw["particleFileFields"].insert(pfw["particleFileFields"].index("RVector"), "CZTag")

# Boundary face order is x-, x+, y-, y+, z-, z+.  With planeStrain enabled,
# the second transverse entry on y-faces controls the x velocity.  This is a
# standard Couette setup: fixed lower platen and moving upper platen.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Linear"
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0, 1.0, 1.0]]
pfw["enablePrescribedBoundaryTransverseVelocities"] = [0, 0, 1, 1, 0, 0]
pfw["prescribedBoundaryTransverseVelocities"] = [
    [0.0, 0.0],
    [0.0, 0.0],
    [0.0, 0.0],
    [0.0, top_plate_velocity],
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

platen_material_name = "weakTraceCouettePlaten"
fluid_material_name = "weakTraceCouetteFluid"
platen_young = float(os.environ.get("WEAK_TRACE_PLATEN_YOUNG", "20.0"))
platen_poisson = float(os.environ.get("WEAK_TRACE_PLATEN_POISSON", "0.25"))
platen_density = float(os.environ.get("WEAK_TRACE_PLATEN_DENSITY", "2.0"))
fluid_density = float(os.environ.get("WEAK_TRACE_FLUID_DENSITY", "1.0"))
fluid_bulk = float(os.environ.get("WEAK_TRACE_FLUID_BULK", "2.2"))
fluid_viscosity = float(os.environ.get("WEAK_TRACE_FLUID_VISCOSITY", "0.1"))
pfw["materials"] = [platen_material_name, fluid_material_name]
pfw["materialPropertyString"] = elastic_xml(platen_material_name, platen_density, platen_young, platen_poisson)
pfw["materialPropertyString"] += "\n" + liquid_xml(fluid_material_name, fluid_density, fluid_bulk, fluid_viscosity)
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

# Cubic B-spline particles are the primary C0 mortar path. Periodic B-spline
# stencils require at least two MPI partitions in x; the default layout uses
# three. CPDI remains selectable with WEAK_TRACE_PARTICLE_TYPE=2.
particle_type = int(os.environ.get("WEAK_TRACE_PARTICLE_TYPE", "1"))
if pfw["xpar"] < 2 and (particle_type == 1 or VARIANT == "traceContactGroups"):
    raise ValueError(
        "Periodic SinglePointBSpline and C0-mortar Couette runs require "
        "WEAK_TRACE_XPAR >= 2"
    )


def interface_group(material_group: int) -> int:
    return material_group if VARIANT == "traceContactGroups" else 0


def mark_interface(name: str, obj):
    if VARIANT in ("singleField", "traceContactGroups"):
        return geom.surfaceFlagWrapper(name + "_weak", obj, int(geom.SurfaceFlag.WeakDiscontinuity))
    if VARIANT == "falseElasticCZ":
        return geom.czTagWrapper(name + "_cz_tag", geom.surfaceFlagWrapper(name + "_cz", obj, int(geom.SurfaceFlag.Cohesive)), 0)
    return obj


lower = geom.box(
    "lower_platen",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], fluid_ymin, pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=interface_group(0),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, False, False, True],
)
fluid = geom.box(
    "fluid_layer",
    [pfw["xmin"], fluid_ymin, pfw["zmin"]],
    [pfw["xmax"], fluid_ymax, pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=1,
    group=interface_group(1),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, False, True, True],
)
upper = geom.box(
    "upper_platen",
    [pfw["xmin"], fluid_ymax, pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=interface_group(0),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, False, True, False],
)
pfw["objects"] = [
    mark_interface("lower_platen_interface", lower),
    mark_interface("fluid_interface", fluid),
    mark_interface("upper_platen_interface", upper),
]

tracer_y = [
    pfw["ymin"] + 0.5 * platen_thickness,
    fluid_ymin + 0.10 * fluid_gap,
    fluid_ymin + 0.25 * fluid_gap,
    fluid_ymin + 0.50 * fluid_gap,
    fluid_ymin + 0.75 * fluid_gap,
    fluid_ymin + 0.90 * fluid_gap,
    pfw["ymax"] - 0.5 * platen_thickness,
]
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
    "top_plate_velocity_x": top_plate_velocity,
    "domain_width": domain_width,
    "domain_height": domain_height,
    "platen_thickness": platen_thickness,
    "fluid_gap": fluid_gap,
    "expected_fluid_shear_rate": top_plate_velocity / fluid_gap,
    "expected_fluid_shear_stress_xy": fluid_viscosity * top_plate_velocity / fluid_gap,
    "platen_young": platen_young,
    "platen_poisson": platen_poisson,
    "platen_density": platen_density,
    "fluid_density": fluid_density,
    "fluid_bulk": fluid_bulk,
    "fluid_viscosity": fluid_viscosity,
    "particle_type": particle_type,
    "update_method": pfw["updateMethod"],
    "update_order": pfw["updateOrder"],
    "xpar": pfw["xpar"],
    "ypar": pfw["ypar"],
    "nI": pfw["nI"],
    "nJ": pfw["nJ"],
    "surface_refinement_level": surface_refinement_level,
}
