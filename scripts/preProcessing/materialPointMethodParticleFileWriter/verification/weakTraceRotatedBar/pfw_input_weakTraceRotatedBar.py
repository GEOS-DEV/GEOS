# -*- coding: utf-8 -*-
"""Rotated two-material bar for weak-interface trace projection.

A diagonal two-phase coupon is loaded in plane-strain tension by moving the
x- and y-F-table faces.  The coupon uses two VonMisesJ materials: a compliant,
strong phase and a stiff, weak phase.  The intended series solution remains
elastic at the prescribed strain, but a single-field mixed-cell kinematic error
can over-stress and plastify the stiff/weak phase near the inclined interface.

The geometry uses standard pfw_geometry polygon objects.  The lower and upper
bar segments are parallelograms in global coordinates that extend into the
lower-left and upper-right domain corners so the moving boundaries load finite
boundary segments rather than corner points.  Only the internal A/B interface
edge is flagged in the CZ and trace variants.
"""

from __future__ import annotations

import math
import os

import numpy as np

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}

VARIANT = os.environ.get("WEAK_TRACE_VARIANT", "singleField")
VARIANT_LABEL = os.environ.get("WEAK_TRACE_VARIANT_LABEL", VARIANT)
CASE_NAME = os.environ.get("WEAK_TRACE_CASE_NAME", f"weakTraceRotatedBar_{VARIANT}")

pfw["caseName"] = CASE_NAME
pfw["runDebug"] = True
pfw["mBatch"] = True
pfw["mWallTime"] = os.environ.get("WEAK_TRACE_WALLTIME", "00:10:00")
pfw["mSubmitJobs"] = False

# Mesh and domain -----------------------------------------------------------------------
refine = int(os.environ.get("WEAK_TRACE_REFINE", "1"))
cpp = int(os.environ.get("WEAK_TRACE_CPP", "40"))
ppc = int(os.environ.get("WEAK_TRACE_PPC", "2"))

pfw["planeStrain"] = 1
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = ppc

pfw["xmin"] = -0.55
pfw["xmax"] = 0.55
pfw["ymin"] = -0.55
pfw["ymax"] = 0.55
thickness = (pfw["xmax"] - pfw["xmin"]) / float(max(1, pfw["nI"] - 2))
pfw["zmin"] = -0.5 * thickness
pfw["zmax"] = 0.5 * thickness

# Time integration and diagnostics ------------------------------------------------------
stop_time = float(os.environ.get("WEAK_TRACE_STOP_TIME", "20.0"))
final_strain = float(os.environ.get("WEAK_TRACE_FINAL_STRAIN", "0.006"))
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 80.0
pfw["restartInterval"] = 10.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = float(os.environ.get("WEAK_TRACE_CFL", "0.18"))
pfw["initialDt"] = 1.0e-16
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 1000.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 1000.0
pfw["frictionCoefficient"] = 0.0
pfw["contactNormalType"] = "Aligned"
pfw["useSurfacePositionForContact"] = 1
pfw["damageFieldPartitioning"] = 0
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
    "gridSurfacePosition",
    "gridWeakInterfaceTraceActive",
    "gridWeakInterfaceTraceContactSuppressed",
    "gridWeakInterfaceTraceSkipReason",
    "gridWeakInterfaceTraceAnchorWeight",
    "gridWeakInterfaceTraceSupportWeight",
    "gridWeakInterfaceTraceSurfaceJump",
    "gridWeakInterfaceTraceVelocityJump",
    "gridWeakInterfaceTraceVelocityJumpPost",
    "gridWeakInterfaceTraceForce",
]

pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "SurfaceFlag",
    "RVector",
    "SurfaceNormal",
    "SurfacePosition",
    "PlasticStrainMagnitude",
]

# Moving F-table faces.  Stretching x and y moves the lower-left and upper-right
# tabs apart along the coupon axis a=(1,1)/sqrt(2).  z is plane strain.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 1, 1]
pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0 + final_strain, 1.0 + final_strain, 1.0]]

# Trace controls ------------------------------------------------------------------------
if VARIANT == "traceContactGroups":
    pfw["enableContact"] = 1
    pfw["enableWeakInterfaceTraceProjection"] = 1
    pfw["weakInterfaceTraceProjectionIterations"] = int(os.environ.get("WEAK_TRACE_PROJECTION_ITERATIONS", "1"))
    pfw["weakInterfaceTraceProjectionScale"] = float(os.environ.get("WEAK_TRACE_PROJECTION_SCALE", "0.25"))
    pfw["weakInterfaceTraceGapStabilization"] = float(os.environ.get("WEAK_TRACE_GAP_STABILIZATION", "0.0"))
    pfw["weakInterfaceTraceMinWeight"] = float(os.environ.get("WEAK_TRACE_MIN_WEIGHT", "1.0e-12"))
    pfw["weakInterfaceTraceSuppressNodalContact"] = 1
    pfw["weakInterfaceTracePairs"] = [[0, 1]]
elif VARIANT == "falseElasticCZ":
    pfw["enableContact"] = 1
    pfw["enableWeakInterfaceTraceProjection"] = 0
else:
    pfw["enableContact"] = 0
    pfw["enableWeakInterfaceTraceProjection"] = 0

# Materials -----------------------------------------------------------------------------
matA = material_db.weakTraceCompliantVonMises.copy()
matB = material_db.weakTraceStiffWeakVonMises.copy()
cz = material_db.weakTraceFalseElasticCohesiveZone.copy()
pfw["materials"] = [matA["name"], matB["name"]]
pfw["materialPropertyString"] = matA["materialString"] + "\n" + matB["materialString"]
if VARIANT == "falseElasticCZ":
    pfw["materialPropertyString"] += "\n" + cz["materialString"]
    pfw["cohesiveZoneRegions"] = """
<CohesiveZoneRegion
    name="weakTraceFalseElasticCZ"
    constitutiveModel="weakTraceFalseElasticCohesiveZone"
    tag="0"/>"""
    pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="weakTraceFalseElasticCZ"
    startTime="0.0"
    regionNames="{ weakTraceFalseElasticCZ }"
    czVolumeNormalization="1"/>"""

# Geometry ------------------------------------------------------------------------------
a = np.array([1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)])
b = np.array([-1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)])
bar_half_width = float(os.environ.get("WEAK_TRACE_BAR_HALF_WIDTH", "0.11"))
s0 = -0.82
s1 = 0.82
# A small nonzero split offset avoids placing the material interface exactly on
# a diagonal row of particle centers for the default structured seeding.  That
# gives both sides a well-populated layer of WeakDiscontinuity particles while
# preserving the non-grid-aligned interface.
s_interface = float(os.environ.get("WEAK_TRACE_INTERFACE_OFFSET", "0.006875"))

def xy(s, t):
    return (s * a + t * b).tolist()

# Vertices are ordered counter-clockwise.  The lower segment interface edge is
# edge 1; the upper segment interface edge is edge 3.  Non-interface faces are
# deliberately unflagged.
lower_vertices = [xy(s0, -bar_half_width), xy(s_interface, -bar_half_width), xy(s_interface, bar_half_width), xy(s0, bar_half_width)]
upper_vertices = [xy(s_interface, -bar_half_width), xy(s1, -bar_half_width), xy(s1, bar_half_width), xy(s_interface, bar_half_width)]

def make_objects(surface_flag=None, trace_groups=False):
    lower_flags = [False, surface_flag is not None, False, False]
    upper_flags = [False, False, False, surface_flag is not None]
    lower_group = 0
    upper_group = 1 if trace_groups else 0
    lower = geom.polygon(
        "weakTraceCompliantLower",
        lower_vertices,
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=lower_group,
        particleType=2,
        flaggedSurfaces=lower_flags,
    )
    upper = geom.polygon(
        "weakTraceStiffUpper",
        upper_vertices,
        vel=[0.0, 0.0, 0.0],
        mat=1,
        group=upper_group,
        particleType=2,
        flaggedSurfaces=upper_flags,
    )
    if surface_flag is not None:
        lower = geom.surfaceFlagWrapper("weakTraceLowerInterfaceFlag", lower, int(surface_flag))
        upper = geom.surfaceFlagWrapper("weakTraceUpperInterfaceFlag", upper, int(surface_flag))
    return [lower, upper]

if VARIANT == "singleField":
    pfw["objects"] = make_objects(surface_flag=None, trace_groups=False)
elif VARIANT == "falseElasticCZ":
    pfw["objects"] = make_objects(surface_flag=geom.SurfaceFlag.Cohesive, trace_groups=False)
elif VARIANT == "traceContactGroups":
    pfw["objects"] = make_objects(surface_flag=geom.SurfaceFlag.WeakDiscontinuity, trace_groups=True)
else:
    # Placeholders currently run as single-field unless explicitly dispatched for
    # development.  Their rows are still reported by the folder postprocessor.
    pfw["objects"] = make_objects(surface_flag=None, trace_groups=False)

# Local tracer layers -------------------------------------------------------------------
tracer_points = []
tracer_labels = []
layer_specs = [
    ("compliant_near", s_interface - 0.035),
    ("stiff_near", s_interface + 0.035),
    ("compliant_far", s_interface - 0.28),
    ("stiff_far", s_interface + 0.28),
]
for label, s in layer_specs:
    for j, t in enumerate([-0.06, 0.0, 0.06]):
        p = s * a + t * b
        tracer_points.append((float(p[0]), float(p[1]), 0.0))
        tracer_labels.append(f"{label}_{j}")

tracers.set_tracers(
    pfw,
    tracer_points,
    variables=[
        "particleID",
        "materialType",
        "density",
        "plasticStrainMagnitude",
        "stressXX",
        "stressYY",
        "stressZZ",
        "stressYZ",
        "stressXZ",
        "stressXY",
    ],
    write_interval=stop_time / 500.0,
    output_prefix="weakTraceRotatedBarTracer",
)

pfw_expected = {
    "variant": VARIANT,
    "variant_label": VARIANT_LABEL,
    "final_strain": final_strain,
    "stop_time": stop_time,
    "bar_axis": [float(a[0]), float(a[1]), 0.0],
    "interface_normal": [float(a[0]), float(a[1]), 0.0],
    "trace_projection_scale": pfw.get("weakInterfaceTraceProjectionScale", 0.0),
    "trace_gap_stabilization": pfw.get("weakInterfaceTraceGapStabilization", 0.0),
    "trace_projection_iterations": pfw.get("weakInterfaceTraceProjectionIterations", 0),
    "tracer_labels": tracer_labels,
}
