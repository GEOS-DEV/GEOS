# -*- coding: utf-8 -*-
"""Large-deformation compaction benchmark for C0 weak-interface reconstruction.

The specimen is a 2-D plane-strain copper compact with two aluminum inclusions
between tool-steel platens.  The copper matrix is one contact group, while the
platens and aluminum inclusions are a second group.  Particle-carried
SurfaceFlag::WeakDiscontinuity data seed a C0 vertex-mortar tie between those
two groups at all bonded interfaces.
"""

from __future__ import annotations

import math
import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db


pfw = {}

CASE_ID = "weakTraceLargeDeformationCompaction"
VARIANT = os.environ.get("WEAK_TRACE_VARIANT", "traceContactGroups")
VARIANT_LABEL = os.environ.get("VV_VARIANT_LABEL", VARIANT)
CASE_NAME = os.environ.get("VV_CASE_NAME", f"{CASE_ID}_{VARIANT}")


def env_int(name: str, default: int) -> int:
    return int(os.environ.get(name, str(default)))


def env_float(name: str, default: float) -> float:
    return float(os.environ.get(name, str(default)))


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


def material_copy(base: dict, name: str, yield_strength: float | None = None) -> dict:
    material = base.copy()
    material["name"] = name
    if yield_strength is not None:
        material["defaultYieldStrength"] = yield_strength
    material["materialString"] = material_db.generateMaterialString(material)
    return material


pfw["caseName"] = CASE_NAME
pfw["runDebug"] = env_bool("WEAK_TRACE_RUN_DEBUG", True)
pfw["mPartition"] = os.environ.get("WEAK_TRACE_PARTITION", "pdebug" if pfw["runDebug"] else "pbatch")
pfw["mBatch"] = True
pfw["mWallTime"] = os.environ.get("WEAK_TRACE_WALLTIME", "00:20:00")
pfw["mSubmitJobs"] = False
pfw["autoRestart"] = False
pfw["lastRestartBufferInSeconds"] = 300

# Mesh and geometry defaults are intentionally small enough for a first
# diagnostic run but expose the same refinement knobs as the other weak-trace
# verification cases.
refine = env_int("WEAK_TRACE_REFINE", 1)
cpp = env_int("WEAK_TRACE_CPP", 56)
ppc = env_int("WEAK_TRACE_PPC", 4)
surface_refinement_level = env_int("WEAK_TRACE_SURFACE_REFINEMENT_LEVEL", 2)
particle_type = env_int("WEAK_TRACE_PARTICLE_TYPE", 1)

domain_width = env_float("WEAK_TRACE_DOMAIN_WIDTH", 1.12)
matrix_width = env_float("WEAK_TRACE_MATRIX_WIDTH", 0.78)
matrix_height = env_float("WEAK_TRACE_MATRIX_HEIGHT", 0.70)
platen_thickness = env_float("WEAK_TRACE_PLATEN_THICKNESS", 0.12)
side_clearance = 0.5 * max(domain_width - matrix_width, 0.0)
total_height = matrix_height + 2.0 * platen_thickness

pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
pfw["xpar"] = 1
pfw["ypar"] = 1
pfw["zpar"] = 1
pfw["nI"] = refine * cpp
pfw["nJ"] = refine * cpp
pfw["nK"] = 3
pfw["ppc"] = ppc
pfw["surfaceRefinementLevel"] = surface_refinement_level

pfw["xmin"] = -0.5 * domain_width
pfw["xmax"] = 0.5 * domain_width
pfw["ymin"] = -0.5 * total_height
pfw["ymax"] = 0.5 * total_height
dx = (pfw["xmax"] - pfw["xmin"]) / float(max(1, pfw["nI"] - 2))
dy = (pfw["ymax"] - pfw["ymin"]) / float(max(1, pfw["nJ"] - 2))
cell_size = min(dx, dy)
pfw["zmin"] = -0.5 * cell_size
pfw["zmax"] = 0.5 * cell_size

stop_time = env_float("WEAK_TRACE_STOP_TIME", 12.0)
final_y_stretch = env_float("WEAK_TRACE_FINAL_Y_STRETCH", 0.62)
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 60.0
pfw["restartInterval"] = 10.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = env_float("WEAK_TRACE_CFL", 0.16)
pfw["initialDt"] = 1.0e-16
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 600.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 600.0
pfw["profileHistory"] = 1
pfw["profileDirection"] = "y"
pfw["profileNumSlices"] = 0
pfw["profileWriteInterval"] = pfw["plotInterval"]
pfw["profileCycleInterval"] = 0
pfw["profileVariables"] = "{ velocityY, stressXX, stressYY, stressXY }"

pfw["solverProfiling"] = 0
pfw["needsNeighborList"] = 0
pfw["frictionCoefficient"] = 0.0
pfw["contactNormalType"] = "Aligned"
pfw["useSurfacePositionForContact"] = 1
pfw["damageFieldPartitioning"] = 0
pfw["weakDiscontinuityInsertion"] = "off"
pfw["weakInterfaceTraceGeometry"] = "surface"
pfw["cpdiDomainScaling"] = 0
pfw["flagParticlesWithBadMappingArraysForDeletion"] = 1
pfw["plotUnscaledParticles"] = 1
pfw["maxParticleVelocity"] = env_float("WEAK_TRACE_MAX_PARTICLE_VELOCITY", 10.0)
pfw["minParticleJacobian"] = env_float("WEAK_TRACE_MIN_PARTICLE_JACOBIAN", 0.01)
pfw["maxParticleJacobian"] = env_float("WEAK_TRACE_MAX_PARTICLE_JACOBIAN", 10.0)

pfw["updateMethod"] = os.environ.get("WEAK_TRACE_UPDATE_METHOD", "PIC")
pfw["enableWeakInterfacePrimalEnrichment"] = 0

if VARIANT == "singleField":
    pfw["enableContact"] = 0
    pfw["enableWeakInterfaceTraceProjection"] = 0
elif VARIANT == "traceContactGroups":
    # This benchmark exercises the duplicated-field C0 vertex-mortar path.
    # Copper is contact group 1; tool-steel platens and aluminum inclusions are
    # group 0.  The single-field variant below collapses these groups back to 0.
    pfw["enableContact"] = 1
    pfw["enableWeakInterfaceTraceProjection"] = 1
    pfw["explicitSurfaceNormalInfluence"] = env_float("WEAK_TRACE_EXPLICIT_NORMAL_INFLUENCE", 25.0)
    pfw["weakInterfaceTraceProjectionIterations"] = env_int("WEAK_TRACE_PROJECTION_ITERATIONS", 30)
    pfw["weakInterfaceTraceProjectionScale"] = env_float("WEAK_TRACE_PROJECTION_SCALE", 0.6)
    pfw["weakInterfaceTraceSurfaceFitNeighbors"] = env_int("WEAK_TRACE_SURFACE_FIT_NEIGHBORS", 4)
    pfw["weakInterfaceTraceMortarQuadratureOrder"] = env_int("WEAK_TRACE_MORTAR_QUADRATURE_ORDER", 2)
    pfw["weakInterfaceTraceMassCutoffRatio"] = env_float("WEAK_TRACE_MASS_CUTOFF_RATIO", 1.0e-12)
    pfw["weakInterfaceTraceGapStabilization"] = env_float("WEAK_TRACE_GAP_STABILIZATION", 0.0)
    pfw["weakInterfaceTraceMaxGapCorrectionVelocity"] = env_float("WEAK_TRACE_MAX_GAP_CORRECTION_VELOCITY", 0.0)
    pfw["weakInterfaceTraceMinWeight"] = env_float("WEAK_TRACE_MIN_WEIGHT", 1.0e-12)
    pfw["weakInterfaceTraceCutTolerance"] = env_float("WEAK_TRACE_CUT_TOLERANCE", 1.0e-8)
    pfw["weakInterfaceTraceMinSegmentFraction"] = env_float("WEAK_TRACE_MIN_SEGMENT_FRACTION", 1.0e-4)
    pfw["weakInterfaceTraceLinearSolverTolerance"] = env_float("WEAK_TRACE_LINEAR_SOLVER_TOLERANCE", 1.0e-10)
    pfw["weakInterfaceTraceLinearSolverRegularization"] = env_float("WEAK_TRACE_LINEAR_SOLVER_REGULARIZATION", 1.0e-12)
    pfw["weakInterfaceTraceSuppressNodalContact"] = 1
    pfw["weakInterfaceTracePairs"] = [[0, 1]]
else:
    raise ValueError(
        "weakTraceLargeDeformationCompaction supports WEAK_TRACE_VARIANT=singleField or traceContactGroups"
    )

pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = [
    "gridMass",
    "gridVelocity",
    "gridMomentum",
    "gridInternalForce",
    "gridActive",
    "gridSurfaceFieldMass",
    "gridSurfaceNormal",
    "gridExplicitSurfaceNormal",
    "gridPrincipalExplicitSurfaceNormal",
    "gridSurfacePosition",
    "gridContactForce",
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
    "particleStress",
    "particleVelocity",
    "particleReferencePosition",
    "particleMaterialType",
    "particleSurfaceFlag",
    "particleSurfaceNormal",
    "particleSurfacePosition",
    "particleReferenceSurfaceNormal",
    "particleReferenceSurfacePosition",
    "particlePlasticStrain",
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

# Symmetric displacement-controlled compaction.  The top and bottom y faces move
# inward, while lateral faces remain free so copper can flow into the side void.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [
    [0.0, 1.0, 1.0, 1.0],
    [stop_time, 1.0, final_y_stretch, 1.0],
]

# Material properties ------------------------------------------------------------------
platen = material_copy(
    material_db.toolSteelA2,
    "weakTraceCompactionToolSteel",
    env_float("WEAK_TRACE_PLATEN_YIELD", 20.0),
)
copper = material_copy(
    material_db.copper,
    "weakTraceCompactionAnnealedCopper",
    env_float("WEAK_TRACE_COPPER_YIELD", material_db.copper["defaultYieldStrength"]),
)
aluminum = material_copy(
    material_db.aluminum,
    "weakTraceCompactionAnnealedAluminum",
    env_float("WEAK_TRACE_ALUMINUM_YIELD", material_db.aluminum["defaultYieldStrength"]),
)
pfw["materials"] = [platen["name"], copper["name"], aluminum["name"]]
pfw["materialPropertyString"] = "\n".join(
    [platen["materialString"], copper["materialString"], aluminum["materialString"]]
)

# Geometry -----------------------------------------------------------------------------
matrix_xmin = -0.5 * matrix_width
matrix_xmax = 0.5 * matrix_width
matrix_ymin = -0.5 * matrix_height
matrix_ymax = 0.5 * matrix_height
platen_xmin = pfw["xmin"]
platen_xmax = pfw["xmax"]
z_pad = env_float("WEAK_TRACE_GEOMETRY_Z_PAD", 2.0 * cell_size)
cyl_zmin = pfw["zmin"] - z_pad
cyl_zmax = pfw["zmax"] + z_pad
inclusion_radius = env_float("WEAK_TRACE_INCLUSION_RADIUS", 0.105)
inclusion_centers = [
    (
        env_float("WEAK_TRACE_INCLUSION_0_X", -0.16),
        env_float("WEAK_TRACE_INCLUSION_0_Y", 0.15),
    ),
    (
        env_float("WEAK_TRACE_INCLUSION_1_X", 0.18),
        env_float("WEAK_TRACE_INCLUSION_1_Y", -0.14),
    ),
]


def weak_interface(name: str, obj):
    return geom.surfaceFlagWrapper(name + "_weak", obj, int(geom.SurfaceFlag.WeakDiscontinuity))


def interface_group(group: int) -> int:
    return group if VARIANT == "traceContactGroups" else 0


def inclusion_cut(name: str, center: tuple[float, float]):
    return geom.cylinder(
        name,
        [center[0], center[1], cyl_zmin],
        [center[0], center[1], cyl_zmax],
        inclusion_radius,
        vel=[0.0, 0.0, 0.0],
        mat=2,
        group=interface_group(0),
        particleType=particle_type,
        flaggedSurfaces=[False, False, True, False],
    )


matrix_box = geom.box(
    "copper_matrix_box",
    [matrix_xmin, matrix_ymin, pfw["zmin"]],
    [matrix_xmax, matrix_ymax, pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=1,
    group=interface_group(1),
    particleType=particle_type,
    dim=2,
    flaggedSurfaces=[False, True, False, True],
)
matrix_with_first_hole = geom.difference(
    "copper_matrix_minus_inclusion_0",
    matrix_box,
    inclusion_cut("inclusion_0_cut", inclusion_centers[0]),
)
matrix = weak_interface(
    "copper_matrix_interfaces",
    geom.difference(
        "copper_matrix_with_inclusion_holes",
        matrix_with_first_hole,
        inclusion_cut("inclusion_1_cut", inclusion_centers[1]),
    ),
)

aluminum_inclusions = []
for index, center in enumerate(inclusion_centers):
    aluminum_inclusions.append(
        weak_interface(
            f"aluminum_inclusion_{index}_interface",
            geom.cylinder(
                f"aluminum_inclusion_{index}",
                [center[0], center[1], cyl_zmin],
                [center[0], center[1], cyl_zmax],
                inclusion_radius,
                vel=[0.0, 0.0, 0.0],
                mat=2,
                group=interface_group(0),
                particleType=particle_type,
                flaggedSurfaces=[False, False, True, False],
            ),
        )
    )

def platen_box(name: str, xmin: float, ymin: float, xmax: float, ymax: float, flagged_surfaces: list[bool]):
    return geom.box(
        name,
        [xmin, ymin, pfw["zmin"]],
        [xmax, ymax, pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=interface_group(0),
        particleType=particle_type,
        dim=2,
        flaggedSurfaces=flagged_surfaces,
    )


# Only the platen faces directly opposite copper should seed weak-interface
# geometry.  The overhanging platen faces above the side void are ordinary free
# surfaces; flagging them creates one-sided weak interfaces and corner artifacts.
unflagged_faces = [False, False, False, False]
bottom_platen_left = platen_box(
    "bottom_platen_left_free",
    platen_xmin,
    pfw["ymin"],
    matrix_xmin,
    matrix_ymin,
    unflagged_faces,
)
bottom_platen_center = weak_interface(
    "bottom_platen_center_top_interface",
    platen_box(
        "bottom_platen_center",
        matrix_xmin,
        pfw["ymin"],
        matrix_xmax,
        matrix_ymin,
        [False, False, False, True],
    ),
)
bottom_platen_right = platen_box(
    "bottom_platen_right_free",
    matrix_xmax,
    pfw["ymin"],
    platen_xmax,
    matrix_ymin,
    unflagged_faces,
)
top_platen_left = platen_box(
    "top_platen_left_free",
    platen_xmin,
    matrix_ymax,
    matrix_xmin,
    pfw["ymax"],
    unflagged_faces,
)
top_platen_center = weak_interface(
    "top_platen_center_bottom_interface",
    platen_box(
        "top_platen_center",
        matrix_xmin,
        matrix_ymax,
        matrix_xmax,
        pfw["ymax"],
        [False, True, False, False],
    ),
)
top_platen_right = platen_box(
    "top_platen_right_free",
    matrix_xmax,
    matrix_ymax,
    platen_xmax,
    pfw["ymax"],
    unflagged_faces,
)

pfw["objects"] = [
    bottom_platen_left,
    bottom_platen_center,
    bottom_platen_right,
    matrix,
] + aluminum_inclusions + [
    top_platen_left,
    top_platen_center,
    top_platen_right,
]

# Tracers are sparse and only meant to catch gross large-strain pathologies in
# the aggregate report.  Visit frames remain the primary diagnostic for this
# qualitative benchmark.
tracer_points = [
    (0.0, matrix_ymin + 0.04, 0.0),
    (0.0, matrix_ymax - 0.04, 0.0),
    (matrix_xmin + 0.04, 0.0, 0.0),
    (matrix_xmax - 0.04, 0.0, 0.0),
    (inclusion_centers[0][0], inclusion_centers[0][1], 0.0),
    (inclusion_centers[1][0], inclusion_centers[1][1], 0.0),
]
for cx, cy in inclusion_centers:
    for theta in (0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi):
        tracer_points.append(
            (
                cx + (inclusion_radius + 0.5 * cell_size) * math.cos(theta),
                cy + (inclusion_radius + 0.5 * cell_size) * math.sin(theta),
                0.0,
            )
        )
for x in (-0.5 * (matrix_width + side_clearance), 0.5 * (matrix_width + side_clearance)):
    tracer_points.append((x, 0.0, 0.0))

tracers.set_tracers(
    pfw,
    tracer_points,
    [
        "particleID",
        "materialType",
        "velocityX",
        "velocityY",
        "stressXX",
        "stressYY",
        "stressXY",
        "plasticStrainMagnitude",
    ],
    write_interval=stop_time / 600.0,
    output_prefix="weakTraceLargeDeformationCompactionTracer",
)

pfw_expected = {
    "case_id": CASE_ID,
    "variant": VARIANT,
    "variant_label": VARIANT_LABEL,
    "domain_width": domain_width,
    "matrix_width": matrix_width,
    "matrix_height": matrix_height,
    "platen_thickness": platen_thickness,
    "side_clearance_each_side": side_clearance,
    "initial_sample_area": matrix_width * matrix_height,
    "initial_side_void_area": max(domain_width - matrix_width, 0.0) * matrix_height,
    "final_y_stretch": final_y_stretch,
    "nominal_height_reduction": 1.0 - final_y_stretch,
    "inclusion_radius": inclusion_radius,
    "inclusion_centers": inclusion_centers,
    "particle_type": particle_type,
    "ppc": ppc,
    "surface_refinement_level": surface_refinement_level,
    "run_debug": pfw["runDebug"],
    "partition": pfw["mPartition"],
    "wall_time": pfw["mWallTime"],
    "update_method": pfw["updateMethod"],
    "trace_projection": pfw.get("enableWeakInterfaceTraceProjection", 0),
    "trace_pairs": pfw.get("weakInterfaceTracePairs", []),
    "trace_projection_iterations": pfw.get("weakInterfaceTraceProjectionIterations", None),
    "trace_projection_scale": pfw.get("weakInterfaceTraceProjectionScale", None),
    "trace_surface_fit_neighbors": pfw.get("weakInterfaceTraceSurfaceFitNeighbors", None),
    "trace_mortar_quadrature_order": pfw.get("weakInterfaceTraceMortarQuadratureOrder", None),
    "trace_suppress_nodal_contact": pfw.get("weakInterfaceTraceSuppressNodalContact", None),
    "materials": pfw["materials"],
}
