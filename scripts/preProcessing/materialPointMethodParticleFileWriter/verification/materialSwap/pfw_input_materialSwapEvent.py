# -*- coding: utf-8 -*-
"""MaterialSwap event verification.

A small elastic block starts in ParticleRegion1 with material type 0 and zero
velocity.  At a prescribed time the MaterialSwap event transfers the particles
to ParticleRegion2, whose material has the same density but a different elastic
stiffness.  There is no applied load, so the analytical kinematic solution is
zero displacement, zero velocity, and zero stress.  The event-specific expected
solution is a step in the tracer particle material type from 0 to 1 at the swap
time, with conserved particle position and momentum.
"""

import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}

case_name = os.environ.get("MATERIAL_SWAP_CASE_NAME", "materialSwapEvent")
variant_label = os.environ.get("MATERIAL_SWAP_VARIANT_LABEL", "identity kinematics")
swap_time = float(os.environ.get("MATERIAL_SWAP_TIME", "0.50"))
stop_time = float(os.environ.get("MATERIAL_SWAP_STOP_TIME", "1.00"))

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("MATERIAL_SWAP_REFINE", "1"))
cpp = int(os.environ.get("MATERIAL_SWAP_CPP", "6"))
ppc = int(os.environ.get("MATERIAL_SWAP_PPC", "2"))
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

pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 4.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 100.0
pfw["outputType"] = "silo"
pfw["particleFileFields"] = ["Velocity", "MaterialType", "ContactGroup", "SurfaceFlag", "RVector"]
pfw["plottableFields"] = "{ particleMaterialType, particleVelocity, particleDisplacement, particleDensity }"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridMomentum", "gridVelocity"]

soft = material_db.verificationMaterialSwapSoft.copy()
stiff = material_db.verificationMaterialSwapStiff.copy()
pfw["materials"] = [soft["name"], stiff["name"]]
pfw["materialPropertyString"] = soft["materialString"] + stiff["materialString"]

source_block = geom.box(
    "source_block",
    [-0.25, -0.25, -0.25],
    [0.25, 0.25, 0.25],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)

# A zero-volume box adds the matching CPDI subregion for the destination
# material without creating destination particles at t=0.
empty_destination = geom.box(
    "empty_destination_region",
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    vel=[0.0, 0.0, 0.0],
    mat=1,
    group=0,
)
pfw["objects"] = [source_block, empty_destination]

tracers.set_tracers(
    pfw,
    points=[[0.0, 0.0, 0.0]],
    variables=["particleID", "particleMaterialType", "velocityX", "velocityY", "velocityZ"],
    write_interval=stop_time / 100.0,
    output_prefix="materialSwap_center",
)

pfw["useEvents"] = 1
pfw["mpmEventsString"] = f"""
<MaterialSwap
  startTime="{swap_time}"
  sourceRegion="ParticleRegion1"
  destinationRegion="ParticleRegion2"/>
"""

pfw_expected = {
    "variant_label": variant_label,
    "swap_time": swap_time,
    "stop_time": stop_time,
    "expected_initial_material_type": 0,
    "expected_final_material_type": 1,
    "expected_max_displacement": 0.0,
    "source_material": soft["name"],
    "destination_material": stiff["name"],
}
