# -*- coding: utf-8 -*-
"""Spinning-disk conservation verification.

A two-dimensional elastic disk is initialized in solid-body rotation and then
allowed to evolve without external force or contact.  The ideal solution is
constant angular momentum.  The folder runTest dispatches PIC, FLIP, and FMPM2
variants from this one input so the post-processor can compare their angular
momentum and kinetic-energy drift.
"""

import os
import numpy as np

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
case_name = os.environ.get("SPINNING_DISK_CASE_NAME", "spinningDiskConservation")
variant_label = os.environ.get("SPINNING_DISK_VARIANT_LABEL", "FLIP")
update_method = os.environ.get("SPINNING_DISK_UPDATE_METHOD", "FLIP")
update_order = os.environ.get("SPINNING_DISK_UPDATE_ORDER", "")

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("SPINNING_DISK_REFINE", "1"))
cpp = int(os.environ.get("SPINNING_DISK_CPP", "24"))
ppc = int(os.environ.get("SPINNING_DISK_PPC", "2"))
pfw["planeStrain"] = 1
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

stop_time = 4.0
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 8.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.20
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 200.0
pfw["logMomentum"] = 1
pfw["frictionCoefficient"] = 0.0
pfw["cpdiDomainScaling"] = 1
pfw["maxParticleVelocity"] = 2.0
pfw["updateMethod"] = update_method
if update_order.strip():
    pfw["updateOrder"] = int(update_order)
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

radius = 0.30
omega0 = 0.60

def get_velocity(self, pt):
    x, y, _z = np.array(pt, dtype=float)
    return [-omega0 * y, omega0 * x, 0.0]

disk = geom.cylinder(
    "spinning_disk",
    [0.0, 0.0, pfw["zmin"]],
    [0.0, 0.0, pfw["zmax"]],
    radius,
    vel=get_velocity,
    mat=0,
    group=0,
)
pfw["objects"] = [disk]

tracers.set_tracers(
    pfw,
    tracers.disk((0.0, 0.0, 0.0), radius * 0.90, normal_axis="z", radial_count=3, angular_count=16, include_center=False),
    variables=["particleID", "velocityX", "velocityY", "speed"],
    write_interval=stop_time / 80.0,
    output_prefix="spinningDiskTracer",
)

pfw_expected = {
    "variant_label": variant_label,
    "radius": radius,
    "omega0": omega0,
    "stop_time": stop_time,
}
