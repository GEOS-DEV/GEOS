# -*- coding: utf-8 -*-
"""Periodic free-advection verification.

A disk translates through a two-dimensional periodic domain with no body force,
contact, or prescribed deformation.  The analytical solution is constant
velocity and position modulo the periodic box.  The test isolates periodic MPM
communication, particle remapping across periodic seams, tracer output, and
momentum conservation for a force-free body.
"""

import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
case_name = os.environ.get("PERIODIC_CASE_NAME", "periodicAdvection")
variant_label = os.environ.get("PERIODIC_VARIANT_LABEL", "2D disk")
velocity = [float(x) for x in os.environ.get("PERIODIC_VELOCITY", "0.42,0.31,0.0").split(",")]

pfw["caseName"] = case_name
pfw["runDebug"] = True
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("PERIODIC_REFINE", "3"))
cpp = int(os.environ.get("PERIODIC_CPP", "8"))
ppc = int(os.environ.get("PERIODIC_PPC", "2"))
pfw["planeStrain"] = 1
pfw["periodic"] = [True, True, False]
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = ppc

pfw["xmin"] = -0.5
pfw["xmax"] = 0.5
pfw["ymin"] = -0.5
pfw["ymax"] = 0.5
thickness = (pfw["xmax"] - pfw["xmin"]) / float(pfw["nI"] - 2)
pfw["zmin"] = -0.5 * thickness
pfw["zmax"] = 0.5 * thickness

stop_time = 2.0
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 8.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 200.0
pfw["logMomentum"] = 1
pfw["frictionCoefficient"] = 0.0
pfw["updateMethod"] = "FLIP"
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity", "gridMomentum"]
pfw["particleFileFields"] = ["Velocity", "MaterialType"]

pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 0, 0, 1, 1]

material = material_db.verificationElastic.copy()
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

center0 = (-0.18, -0.16, 0.0)
radius = 0.12
disk = geom.cylinder(
    "periodic_disk",
    [center0[0], center0[1], pfw["zmin"]],
    [center0[0], center0[1], pfw["zmax"]],
    radius,
    vel=velocity,
    mat=0,
    group=0,
)
pfw["objects"] = [disk]

tracers.set_tracers(
    pfw,
    [center0] + tracers.disk(center0, radius * 0.75, normal_axis="z", radial_count=2, angular_count=12, include_center=False),
    variables=["particleID", "velocityX", "velocityY", "speed"],
    write_interval=stop_time / 80.0,
    output_prefix="periodicAdvectionTracer",
)

pfw_expected = {
    "variant_label": variant_label,
    "center0": center0,
    "velocity": velocity,
    "domain": [pfw["xmin"], pfw["xmax"], pfw["ymin"], pfw["ymax"]],
    "stop_time": stop_time,
}
