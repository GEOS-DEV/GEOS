# -*- coding: utf-8 -*-
"""
Brazilian disk XPIC example for GEOS MPM.

This file is intentionally only a particle-file-writer input. It does not know
about the example-suite driver, run directories, Slurm post-processing, VisIt, or
report generation. A user should be able to copy this file, edit dimensions,
materials, loading, and update method, and run it with particleFileWriter.py.

Problem summary
---------------
A circular disk is embedded in a rectangular plane-strain domain. The x faces are
periodic so the disk can translate horizontally. The y faces are moving
boundaries driven by an F-table; the y scale factor decreases from 1.0 to 0.8,
which dynamically compresses the disk. A Voronoi-Weibull wrapper assigns a
spatially variable strength scale to the quartz material so that damage localizes
and fragments the disk.

Useful source references while editing this file
------------------------------------------------
* particleFileWriter.py defines the accepted pfw dictionary keys and writes the
  GEOS XML, particle file, batch script, and optional submission command.
* pfw_geometryObjects.py defines geometry objects such as cylinder and wrappers
  such as voronoiWeibullBoxWrapper.
* The GEOS MPM solver consumes the generated XML attributes. For example,
  boundaryConditionTypes uses x-, x+, y-, y+, z-, z+ order and the integer
  options Outflow=0, Symmetry=1, Moving=2, Contact=3.
"""

import importlib

import numpy as np

import pfw_geometryObjects as geom


# The particle file writer imports this module and reads the pfw dictionary.
# Keys with XML-facing behavior are mirrored into the generated GEOS input file;
# keys that control discretization/batch behavior are consumed by PFW itself.
pfw = {}

# runDebug=True selects a debug-oriented partition/default behavior in
# particleFileWriter.py. Set to False for production-scale runs where the user
# wants to choose a normal batch queue/partition.
pfw["runDebug"] = True

# End time for this compact example. The wrapper workflow renders the first and
# last plot states, so plotInterval is set equal to stopTime below. Restart
# output is disabled by setting restartInterval greater than the end time.
stopTime = 1.0


# -----------------------------------------------------------------------------
# Material model
# -----------------------------------------------------------------------------
# The [pfw_dependency] tag tells particleFileWriter.py to copy this dependency
# into the run directory before executing the generated problem. The quartz
# dictionary supplies both a material name and a GEOS XML material-property
# string for the ceramic damage model.
# [pfw_dependency] /pfw_materials.py
matFile = importlib.import_module("pfw_materials")
quartz = matFile.quartz

# PFW writes the list of material names and the raw material XML fragment into
# the generated GEOS XML. Geometry objects refer to materials by zero-based
# integer index; mat=0 below selects quartz.
pfw["materials"] = [quartz["name"]]
pfw["materialPropertyString"] = quartz["materialString"]


# -----------------------------------------------------------------------------
# Computational domain and particle/grid resolution
# -----------------------------------------------------------------------------
# refine is the number of MPI partitions in x and y. With zpar=1 this example
# runs on refine*refine MPI ranks. cpp is cells per partition. Decreasing cpp
# makes the example cheaper; increasing cpp gives finer crack/damage resolution.
refine = 3
cpp = 24

pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1

# The domain is 1.5 times wider than tall, so use 1.5*cpp cells per partition in
# x and cpp cells per partition in y. nK must be 3 for planeStrain=1; PFW checks
# this and aborts if plane-strain inputs use a different nK.
pfw["nI"] = round(pfw["xpar"] * cpp * 1.5)
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3

# ppc is particles per cell per active direction. In plane strain PFW creates a
# single layer of particles through thickness, but ppc still controls in-plane
# particle density.
pfw["ppc"] = 2

domainHeight = 1.0
domainWidth = 1.5 * domainHeight

# Choose the z thickness so the generated background cells are approximately
# cubic before the plane-strain reduction is applied.
domainLength = domainHeight * pfw["nK"] / pfw["nJ"]

# DX is used below as the characteristic length for the Weibull flaw spacing.
DX = domainWidth / pfw["nI"]

# Physical domain bounds, before PFW adds ghost cells for non-periodic faces.
pfw["xmin"] = -0.5 * domainWidth
pfw["xmax"] = 0.5 * domainWidth
pfw["ymin"] = 0.0
pfw["ymax"] = domainHeight
pfw["zmin"] = -0.5 * domainLength
pfw["zmax"] = 0.5 * domainLength


# -----------------------------------------------------------------------------
# Batch settings used by particleFileWriter.py
# -----------------------------------------------------------------------------
# mBatch=True asks PFW to write a Slurm batch script. mSubmitJobs=True submits it
# immediately. The bank/account, GEOS executable, default run directory, and
# Python command come from userDefs_$USER.py, not from this input file.
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mCores"] = pfw["xpar"] * pfw["ypar"] * pfw["zpar"]
pfw["mSubmitJobs"] = True

# This example is short enough that the wrapper can simply report failure rather
# than asking PFW to resubmit from restart files.
pfw["autoRestart"] = False


# -----------------------------------------------------------------------------
# GEOS MPM solver controls
# -----------------------------------------------------------------------------
# The generated GEOS event schedule runs until endTime. With plotInterval equal
# to endTime, GEOS writes the initial and final plot states used by VisIt. The
# restart interval is intentionally beyond the end time so no restart files are
# needed for this demonstration problem.
pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime
pfw["restartInterval"] = 2.0 * stopTime

# Silo output is the format used by the VisIt renderer. Use "vtk" instead when
# preparing an example for ParaView.
pfw["outputType"] = "silo"

# ExplicitDynamic is the MPM time integrator used for this dynamic compression
# and fragmentation calculation.
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1e-16

# cpdiDomainScaling enables CPDI domain scaling in GEOS MPM. The damage-field
# partitioning option supports separate damage evolution/averaging in the MPM
# fields used by the ceramic damage model.
pfw["cpdiDomainScaling"] = 1
pfw["damageFieldPartitioning"] = 1

# Plane strain uses one active particle layer through thickness. PFW requires
# nK=3 and zpar=1 for this mode.
pfw["planeStrain"] = 1

# Update method comparison knob. The four Brazilian disk examples differ only in
# this method name. updateOrder=2 is used for FMPM/XPIC-style second-order
# transfers and is harmless to leave in the compact comparison examples.
pfw["updateMethod"] = "XPIC"
pfw["updateOrder"] = 2

# Enable GEOS profiling output and the neighbor-list infrastructure used by the
# contact/damage feature set in this example.
pfw["solverProfiling"] = 1
pfw["needsNeighborList"] = 1

# Write history files for the post-processing script. reactionHistory.csv is used
# for y-face reaction plots; boxAverageHistory.csv is copied as a diagnostic.
pfw["reactionHistory"] = 1
pfw["boxAverageHistory"] = 1

# Friction coefficient for contact-style interactions in the MPM solver.
pfw["frictionCoefficient"] = 0.25

# Particle fields written by PFW into the initial particle file. Velocity gives
# the disk an initial horizontal translation. MaterialType maps particles to
# quartz. ContactGroup separates fields for FMPM/contact logic. StrengthScale
# carries the Weibull variability assigned by the wrapper. SurfaceFlag marks
# boundary/surface particles. RVector stores the reference domain vector used by
# the CPDI particle-domain representation.
pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "StrengthScale",
    "SurfaceFlag",
    "RVector",
]


# -----------------------------------------------------------------------------
# Boundary conditions and loading
# -----------------------------------------------------------------------------
# Do not use a time-dependent boundary-condition-type table; the face types are
# constant throughout the run.
pfw["prescribedBcTable"] = 0

# Boundary-condition order is x-, x+, y-, y+, z-, z+. Integer options in the MPM
# solver are Outflow=0, Symmetry=1, Moving=2, Contact=3. The x faces are also
# periodic, so the x entries are effectively bypassed. The y faces move according
# to the F-table below. The z faces are symmetry planes for plane strain.
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 1, 1]
pfw["periodic"] = [True, False, False]

# prescribedBoundaryFTable enables boundary-driven deformation for faces with
# Moving boundary conditions. The F-table rows are [time, x_scale, y_scale,
# z_scale]. A cosine interpolation avoids a discontinuous acceleration at the
# beginning and end of the loading ramp.
pfw["fTableInterpType"] = "Cosine"
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [
    [0.0,      1.00, 1.00, 1.00],
    [stopTime, 1.00, 0.80, 1.00],
]


# -----------------------------------------------------------------------------
# Geometry and particle attributes
# -----------------------------------------------------------------------------
# A cylinder with its axis through z is a 2D disk in plane strain. x1/x2 are the
# lower/upper axis endpoints, r is the disk radius, vel is the initial particle
# velocity, mat=0 selects quartz from pfw["materials"], and group=0 assigns all
# disk particles to the same contact/material-point group.
disk1 = geom.cylinder(
    "disk1",
    [0.0, domainHeight / 2.0, pfw["zmin"]],
    [0.0, domainHeight / 2.0, pfw["zmax"]],
    domainHeight / 2.0,
    vel=[1.0, 0.0, 0.0],
    mat=0,
    group=0,
)

# Wrap the disk with a Voronoi-Weibull strength-scale field. The wrapper samples
# Voronoi cells over a padded box, clips/assigns the cell information to the disk,
# and computes a strength scale using the material Weibull reference volume and
# modulus. flawSize controls the typical Voronoi spacing; vMin prevents very small
# cells from producing extreme strength scales. dim=3 is used because the
# cylinder and particle positions are three-component quantities even though the
# mechanical response is plane strain.
weibullFlawSize = 6.0 * DX

weibullSample = geom.voronoiWeibullBoxWrapper(
    "weibullSubstrate",
    subObject=disk1,
    x0=np.array([
        pfw["xmin"] - weibullFlawSize,
        pfw["ymin"] - weibullFlawSize,
        pfw["zmin"] - weibullFlawSize,
    ]),
    x1=np.array([
        pfw["xmax"] + weibullFlawSize,
        pfw["ymax"] + weibullFlawSize,
        pfw["zmax"] + weibullFlawSize,
    ]),
    flawSize=weibullFlawSize,
    weibullVolume=quartz["weibullReferenceVolume"],
    weibullModulus=quartz["weibullModulus"],
    weibullSeed=1,
    vMin=DX**3,
    vpts=None,
    dim=3,
    randomMatDir=False,
)

# PFW loops over pfw["objects"] to decide which candidate particles are inside
# the physical material and to query object-provided properties such as velocity,
# material type, contact group, surface flag, and strength scale.
pfw["objects"] = [weibullSample]
