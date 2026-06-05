# -*- coding: utf-8 -*-
"""Partition-crossing free-flight verification.

A disk translates diagonally across a two-dimensional 3 x 3 partitioned domain
without external forces.  The analytical solution is constant velocity and a
linear center-of-mass trajectory.  The test isolates MPI partition ownership
transfer, particle migration, tracer output, and global momentum conservation.
"""

import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
case_name = os.environ.get("PARTITION_CASE_NAME", "partitionCrossing2D")
variant_label = os.environ.get("PARTITION_VARIANT_LABEL", "single material")
velocity = [float(x) for x in os.environ.get("PARTITION_VELOCITY", "0.20,0.16,0.0").split(",")]

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:02:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("PARTITION_REFINE", "3"))
cpp = int(os.environ.get("PARTITION_CPP", "8"))
ppc = int(os.environ.get("PARTITION_PPC", "2"))
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

stop_time = 2.0
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 8.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 200.0
pfw["logMomentum"] = 1
pfw["frictionCoefficient"] = 0.0
pfw["updateMethod"] = os.environ.get("PARTITION_UPDATE_METHOD", "FLIP")
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

center0 = (-0.25, -0.22, 0.0)
radius = 0.09
disk = geom.cylinder(
    "partition_disk",
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
    output_prefix="partitionCrossingTracer",
)

pfw_expected = {
    "variant_label": variant_label,
    "center0": center0,
    "velocity": velocity,
    "stop_time": stop_time,
    "partitions": [pfw["xpar"], pfw["ypar"], pfw["zpar"]],
}
