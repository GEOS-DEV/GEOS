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
pfw["mWallTime"] = "00:10:00"
pfw["mSubmitJobs"] = True

refine = 3
cpp = 15
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppc"] = 2

domainWidth = 1.0
domainHeight = 1.0


pfw["xmin"] = -0.5*domainWidth
pfw["xmax"] =  0.5*domainWidth
pfw["ymin"] = -0.5*domainHeight
pfw["ymax"] =  0.5*domainHeight

domainLength = domainHeight / (pfw["nJ"]-2)
pfw["zmin"] = -0.5*domainLength
pfw["zmax"] =  0.5*domainLength

stop_time = 2.0
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
pfw["particleFileFields"] = ["Velocity", "RVectors", "SurfacePosition", "SurfaceNormal"]

pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 0, 0, 1, 1]

# -1 is void
pfw["enableSurfaceTension"]=1
pfw["surfaceTensionPairs"]=[[0,-1,0.1]] #7.28e-5

rho = 1.0
K = 2.2
mu = 0 # Assume inviscid

pfw["materials"] = ["water"]
pfw["materialPropertyString"] = f"""
<Liquid
    name="water"
    defaultDensity="{rho}" 
    defaultBulkModulus="{K}"
    defaultViscosity="{mu}"/>
"""

drop = geom.cylinder(
    "drop",
    [0.0, 0.0, 10*pfw["zmin"]],
    [0.0, 0.0, 10*pfw["zmax"]],
    r=0.15,
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)
pfw["objects"] = [drop]

pfw_expected = {
    "variant_label": variant_label,
    "domain": [pfw["xmin"], pfw["xmax"], pfw["ymin"], pfw["ymax"]],
    "stop_time": stop_time,
}
