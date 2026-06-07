# -*- coding: utf-8 -*-
"""Elastic stress-control verification.

A single elastic cube is driven by the MPM stress-control boundary condition.
The stress table first applies a hydrostatic compression, then adds a uniaxial
compressive increment in y, and finally unloads to zero.  The expected solution
is the prescribed stress table itself; the verification metric is the tracking
error between box-average stress and the analytical target history.
"""

import os

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}

case_name = os.environ.get("STRESS_CONTROL_CASE_NAME", "elasticStressControl")
variant_label = os.environ.get("STRESS_CONTROL_VARIANT_LABEL", "P-only")
kp = float(os.environ.get("STRESS_CONTROL_KP", "1.0"))
ki = float(os.environ.get("STRESS_CONTROL_KI", "0.0"))
kd = float(os.environ.get("STRESS_CONTROL_KD", "0.0"))

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("STRESS_CONTROL_REFINE", "1"))
cpp = int(os.environ.get("STRESS_CONTROL_CPP", "8"))
ppc = int(os.environ.get("STRESS_CONTROL_PPC", "2"))
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

stop_time = 1.0
pressure = 0.01
y_increment = -0.03
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 5.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.05
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 250.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 250.0
pfw["frictionCoefficient"] = 0.0
pfw["updateMethod"] = "PIC"
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity", "gridInternalForce"]

pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 2, 2]
pfw["prescribedBoundaryFTable"] = 0
pfw["stressControl"] = [1, 1, 1]
pfw["stressTableInterpType"] = "Cosine"
pfw["stressControlKp"] = kp
pfw["stressControlKi"] = ki
pfw["stressControlKd"] = kd
pfw["stressTable"] = [
    [0.00, 0.0, 0.0, 0.0],
    [0.25, -pressure, -pressure, -pressure],
    [0.75, -pressure, -pressure + y_increment, -pressure],
    [1.00, 0.0, 0.0, 0.0],
]

material = material_db.verificationElastic.copy()
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

block = geom.box(
    "stress_control_block",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)
pfw["objects"] = [block]

pfw_expected = {
    "variant_label": variant_label,
    "stress_table": pfw["stressTable"],
    "kp": kp,
    "ki": ki,
    "kd": kd,
}
