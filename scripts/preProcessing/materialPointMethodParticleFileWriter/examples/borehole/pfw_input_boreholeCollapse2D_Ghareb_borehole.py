# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Plane-strain borehole collapse using the Ghareb Geomechanics material dictionary entry.
# Solver/PFW features: InitializeStress, ConfiningPressure, and BoreholePressure MPM events drive the loading.
# Workflow note: keep this file as a copyable problem definition.  Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

import importlib

import numpy as np

import pfw_geometryObjects as geom

# Material database.  The Ghareb material is intentionally a dictionary entry,
# consistent with the other material/model definitions in pfw_materials.py.
# [pfw_dependency] pfw_materials.py
matdb = importlib.import_module("pfw_materials")

pfw = {}

# Units are mm, mg, microsecond, and GPa.
stopTime = 100.0
rampTime = 0.25 * stopTime

# Annular specimen geometry.
domainSize = 1000.0
boreholeDiameter = 250.0
boreholePressure = 0.0
confiningPressure = 0.008

# Mesh and partitioning.  cppx/cppy are kept low so the example suite runs
# quickly; z has one partition and a small number of cells for plane strain.
refine = 3
cppx = 24
cppy = 24
cppz = 3

pfw["planeStrain"] = 1
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cppx
pfw["nJ"] = pfw["ypar"] * cppy
pfw["nK"] = pfw["zpar"] * cppz
pfw["xmin"] = -0.5 * domainSize
pfw["xmax"] = 0.5 * domainSize
pfw["ymin"] = -0.5 * domainSize
pfw["ymax"] = 0.5 * domainSize

# Plane-strain depth is one grid-cell width.
DZ = (pfw["xmax"] - pfw["xmin"]) / (pfw["nI"] - 2)
pfw["zmin"] = -0.5 * DZ
pfw["zmax"] = 0.5 * DZ

# Job/output controls.  The bank and GEOS executable are supplied by userDefs;
# this input only specifies problem-specific non-default run behavior.
pfw["outputType"] = "silo"
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = True
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime
pfw["restartInterval"] = 2.0 * stopTime

# MPM solver options.
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-16
pfw["cpdiDomainScaling"] = 1
pfw["damageFieldPartitioning"] = 1
pfw["maxSingleFieldStateFractionForSeparability"] = 0.999
pfw["frictionCoefficient"] = 0.20
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = "2"
pfw["useEvents"] = 1

pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "StrengthScale",
    "SurfaceFlag",
    "RVector",
]

def make_objects():
    """Create the annular rock region and apply a Weibull strength field."""
    annulus = geom.cylinder(
        "ring",
        x1 = [0.0, 0.0, pfw["zmin"] - 10.0],
        x2 = [0.0, 0.0, pfw["zmax"] + 10.0],
        ri = 0.5 * boreholeDiameter,
        r = 0.95 * domainSize,
        vel = np.array([0.0, 0.0, 0.0]),
        mat = 0,
        group = 0,
    )

    flawSize = 5.0 * DZ
    return [
        geom.voronoiWeibullBoxWrapper(
            "weibullring",
            subObject = annulus,
            x0 = [pfw["xmin"], pfw["ymin"], pfw["zmin"] - flawSize],
            x1 = [pfw["xmax"], pfw["ymax"], pfw["zmax"] + flawSize],
            flawSize = flawSize,
            weibullVolume = 0.25 * np.pi * 25.12 * 25.12 * 56.6,
            weibullModulus = 7.3,
            weibullSeed = 1,
            vMin = DZ**3,
        )
    ]

pfw["objects"] = make_objects()
pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]

# Contact-style boundary handling keeps the annulus in the computational domain
# while the MPM pressure events apply confining and borehole pressures.
pfw["boundaryConditionTypes"] = [1, 1, 1, 1, 1, 1]

# MPMEvents children only.  particleFileWriter.py writes the surrounding
# <MPMEvents> block in the generated GEOS XML.
pfw["mpmEventsString"] = f"""
    <InitializeStress
        startTime="0.0"
        endTime="0.1"
        targetRegion="all"
        pressure="{confiningPressure}"
        />
    <ConfiningPressure
        startTime="0.0"
        endTime="{rampTime}"
        confiningPressureBoxMin="{{{-0.4 * domainSize},{-0.4 * domainSize},{-2.0 * DZ}}}"
        confiningPressureBoxMax="{{{0.4 * domainSize},{0.4 * domainSize},{2.0 * DZ}}}"
        startPressure="{confiningPressure}"
        endPressure="{confiningPressure}"
        interpType="1"
        />
    <BoreholePressure
        startTime="0.0"
        endTime="{rampTime}"
        boreholeRadius="{0.65 * boreholeDiameter}"
        startPressure="{confiningPressure}"
        endPressure="{boreholePressure}"
        interpType="1"
        />
"""
