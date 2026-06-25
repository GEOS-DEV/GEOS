#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


# -*- coding: utf-8 -*-
"""Oscillating droplet with surface tension

A single drop of liquid subject to surface tension
"""

import os

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
case_name = os.environ.get("PERIODIC_CASE_NAME", "2D")
variant_label = os.environ.get("PERIODIC_VARIANT_LABEL", "2D drop")

pfw["caseName"] = case_name
pfw["runDebug"] = True
pfw["mBatch"] = True
pfw["mWallTime"] = "00:30:00"
pfw["mSubmitJobs"] = True

domainWidth = 1.0
domainHeight = 1.0

refine = 3
cpp = 12
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = int(pfw["ypar"] * cpp/2)
pfw["nK"] = 3
pfw["ppc"] = 2

pfw["xmin"] = -0.5*domainWidth
pfw["xmax"] =  0.5*domainWidth
pfw["ymin"] =  0.0
pfw["ymax"] =  domainHeight

domainLength = domainHeight / (pfw["nJ"]-2)
pfw["zmin"] = -0.5*domainLength
pfw["zmax"] =  0.5*domainLength

stop_time = 100.0
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 200
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-12
pfw["reactionHistory"] = 0
pfw["boxAverageHistory"] = 0
pfw["updateMethod"] = "PIC"

pfw["cpdiDomainScaling"]=0
pfw["damageFieldPartitioning"]=0

pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridSurfaceTensionForce"]
pfw["particleFileFields"] = ["Velocity", "ContactGroup", "MaterialType", "RVectors", "SurfacePosition", "SurfaceNormal"]

pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 1, 1]

# -1 is void
pfw["enableSurfaceTension"]=1
pfw["surfaceTensionPairs"]=[[0, -1, 7.28e-5], # Liquid - Gas
                            [1, -1, 1.9e-6], # Solid - Gas For Si 1.14 (111) - 1.9 (110)
                            [0,  1, 8e-5]] # Liquid - solid (6.5e-8)

rho = 1.0
K = 2.2
mu = 0 # Assume inviscid

pfw["materials"] = ["liquid", "surface"]
pfw["materialPropertyString"] = f"""
<Liquid
    name="liquid"
    defaultDensity="{rho}" 
    defaultBulkModulus="{K}"
    defaultViscosity="{mu}"/>
<ElasticIsotropic
    name="surface"
    defaultDensity="2.329"
    defaultYoungModulus="130"
    defaultPoissonRatio="0.33"/>
"""

thickness= 0.2
radius = 0.3
drop = geom.cylinder(
    "drop",
    [0.0, thickness + radius, 10*pfw["zmin"]],
    [0.0, thickness + radius, 10*pfw["zmax"]],
    r=radius,
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)
substrate=geom.box("substrate",
                   [pfw["xmin"], 0.0],
                   [pfw["xmax"], thickness],
                   vel=[0.0, 0.0, 0.0],
                   mat=1,
                   group=1,
                   dim=2,
                   flaggedSurfaces = [ False, False, False, True]
)
pfw["objects"] = [drop, substrate]

pfw_expected = {
    "variant_label": variant_label,
    "domain": [pfw["xmin"], pfw["xmax"], pfw["ymin"], pfw["ymax"]],
    "stop_time": stop_time,
}
