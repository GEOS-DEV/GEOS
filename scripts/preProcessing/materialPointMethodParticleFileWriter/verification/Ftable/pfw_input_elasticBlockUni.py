# -*- coding: utf-8 -*-
"""Elastic F-table boundary-switch verification.

A single elastic cube is stretched in x with a prescribed boundary deformation
history.  The first load segment holds all faces kinematically constrained, so
the expected axial stress follows the constrained modulus.  The second segment
releases the lateral faces with a boundary-condition table, so the incremental
slope should drop to the Young's modulus.  The final segment unloads to zero
strain.  This isolates prescribed F-table interpolation, boundary-condition table
switching, reaction output, and box-average stress output.
"""

import os

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}

case_name = os.environ.get("FTABLE_CASE_NAME", "elasticFTableBoundarySwitch")
variant_label = os.environ.get("FTABLE_VARIANT_LABEL", "Smoothstep")
interp_type = os.environ.get("FTABLE_INTERP_TYPE", "Smoothstep")

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("FTABLE_REFINE", "1"))
cpp = int(os.environ.get("FTABLE_CPP", "8"))
ppc = int(os.environ.get("FTABLE_PPC", "2"))
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = refine
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = pfw["zpar"] * cpp
pfw["ppc"] = ppc

pfw["xmin"] = -0.5
pfw["xmax"] = 0.5
pfw["ymin"] = -0.5
pfw["ymax"] = 0.5
pfw["zmin"] = -0.5
pfw["zmax"] = 0.5

stop_time = 1.25
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 5.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.15
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 250.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 250.0
pfw["frictionCoefficient"] = 0.0
pfw["useInternalForceAsFaceReaction"] = 1
pfw["updateMethod"] = "XPIC"
pfw["updateOrder"] = 2
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity", "gridInternalForce"]

pfw["fTableInterpType"] = interp_type
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [
    [0.00, 1.000, 1.0, 1.0],
    [0.50, 1.010, 1.0, 1.0],
    [1.00, 1.020, 1.0, 1.0],
    [1.25, 1.000, 1.0, 1.0],
]

pfw["prescribedBcTable"] = 1
pfw["bcTable"] = [
    [0.00, 2, 2, 2, 2, 2, 2],
    [0.50, 2, 2, 0, 0, 0, 0],
    [2.00, 2, 2, 0, 0, 0, 0],
]

material = material_db.verificationElastic.copy()
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

block = geom.box(
    "elastic_block",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)
pfw["objects"] = [block]

# Human-readable metadata consumed by the post-processor.
pfw_expected = {
    "variant_label": variant_label,
    "young_modulus": material["defaultYoungModulus"],
    "poisson_ratio": material["defaultPoissonRatio"],
    "release_time": 0.50,
    "stop_time": stop_time,
}
