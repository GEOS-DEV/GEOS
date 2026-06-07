# -*- coding: utf-8 -*-
"""Von Mises elastic-perfectly-plastic uniaxial loading verification.

A single cube is stretched in y with lateral faces released.  For the small
verification strain path, the expected one-dimensional response is linear elastic
until the Von Mises yield strength is reached, followed by an approximately
constant axial stress and increasing equivalent plastic strain.  The test
isolates the VonMisesJ stress update and plastic-strain output.
"""

import os

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
case_name = os.environ.get("VONMISES_CASE_NAME", "vonMisesUniaxialPlasticity")
variant_label = os.environ.get("VONMISES_VARIANT_LABEL", "uniaxial")
final_strain = float(os.environ.get("VONMISES_FINAL_STRAIN", "0.05"))

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("VONMISES_REFINE", "1"))
cpp = int(os.environ.get("VONMISES_CPP", "8"))
ppc = int(os.environ.get("VONMISES_PPC", "2"))
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
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 5.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.10
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 250.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 250.0
pfw["frictionCoefficient"] = 0.0
pfw["updateMethod"] = "PIC"
pfw["outputType"] = "silo"
pfw["particleFileFields"] = ["Velocity", "MaterialType", "PlasticStrainMagnitude"]
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity"]

pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 0, 0]
pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0, 1.0 + final_strain, 1.0]]

material = material_db.verificationVonMises.copy()
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

block = geom.box(
    "von_mises_block",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)
pfw["objects"] = [block]

pfw_expected = {
    "variant_label": variant_label,
    "final_strain": final_strain,
    "young_modulus": material["defaultYoungModulus"],
    "yield_strength": material["defaultYieldStrength"],
}
