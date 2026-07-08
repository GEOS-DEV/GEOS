# -*- coding: utf-8 -*-
"""PFW input for the PBC compaction momentum verification test.

The folder runTest script dispatches this one human-readable input twice, once
with an elastic material and once with an elastic-plastic material.  Environment
variables are used only for the small set of values that define the subcase; the
comments below show the concrete replacements so a user can copy this file and
make a standalone input without using the run script.
"""

import os

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
pfw["runDebug"] = True

# --------------------------------------------------------------------------------------
# Subcase template values set by runTest.
# --------------------------------------------------------------------------------------
# PBC_COMPACTION_VARIANT_LABEL replacement examples:
#   elastic case:  "Elastic"
#   plastic case:  "Elastic-plastic"
variant_label = os.environ.get("PBC_COMPACTION_VARIANT_LABEL", "Elastic")

# PBC_COMPACTION_CASE_NAME replacement examples:
#   elastic case:  "pbcCompactionMomentum_elastic"
#   plastic case:  "pbcCompactionMomentum_plastic"
case_name = os.environ.get("PBC_COMPACTION_CASE_NAME", "pbcCompactionMomentum_elastic")

# PBC_COMPACTION_MATERIAL replacement examples:
#   elastic case:  "pbcCompactionElastic"
#   plastic case:  "pbcCompactionPlastic"
material_name = os.environ.get("PBC_COMPACTION_MATERIAL", "pbcCompactionElastic")

# PBC_COMPACTION_STOP_TIME replacement examples:
#   suite default: "6.0"
#   longer diagnostic run: "30.0"
stop_time = float(os.environ.get("PBC_COMPACTION_STOP_TIME", "6.0"))

# PBC_COMPACTION_GRID_REFINE, PBC_COMPACTION_X_PARTITIONS,
# PBC_COMPACTION_Y_PARTITIONS, and PBC_COMPACTION_CELLS_PER_PARTITION replacement examples:
#   suite default: refine="2", x_partitions="3", y_partitions="2", cells_per_partition="12"
#                  -> 36 x 24 cells on 3 x 2 MPI partitions
#   diagnostic run: refine="5", x_partitions="5", y_partitions="5", cells_per_partition="20"
#                  -> 100 x 100 cells on 5 x 5 MPI partitions
#
# The x-direction is periodic for this test and requires at least three MPI
# partitions in that periodic direction.  The run script may reduce refine for
# fast runs, but xpar must not drop below 3.
grid_refine = int(os.environ.get("PBC_COMPACTION_GRID_REFINE", "2"))
cells_per_partition = int(os.environ.get("PBC_COMPACTION_CELLS_PER_PARTITION", "12"))
periodic_x_partitions = int(os.environ.get("PBC_COMPACTION_X_PARTITIONS", str(max(3, grid_refine))))
y_partitions = int(os.environ.get("PBC_COMPACTION_Y_PARTITIONS", str(max(1, grid_refine))))
if periodic_x_partitions < 3:
    raise ValueError("PBC compaction requires at least 3 MPI partitions in the periodic x direction")
if y_partitions < 1:
    raise ValueError("PBC compaction requires at least 1 MPI partition in the y direction")

# PBC_COMPACTION_FINAL_YSTRETCH replacement examples:
#   suite default: "0.4"  -> 60 percent imposed height reduction at final time
final_y_stretch = float(os.environ.get("PBC_COMPACTION_FINAL_YSTRETCH", "0.4"))

# PBC_COMPACTION_PARTICLE_TYPE replacement examples:
#   suite default: "2" -> CPDI particles
#   optional:       "1" -> SinglePointBspline particles
particle_type = int(os.environ.get("PBC_COMPACTION_PARTICLE_TYPE", "2"))

seed = int(os.environ.get("PBC_COMPACTION_SEED", "12345"))

# --------------------------------------------------------------------------------------
# Domain and particle resolution.
# --------------------------------------------------------------------------------------
domain_height = 1.0
domain_width = 1.0

pfw["caseName"] = case_name
pfw["xmin"] = -0.5 * domain_width
pfw["xmax"] =  0.5 * domain_width
pfw["ymin"] = -0.5 * domain_height
pfw["ymax"] =  0.5 * domain_height
pfw["planeStrain"] = 1
pfw["periodic"] = [True, False, False]

pfw["xpar"] = periodic_x_partitions
pfw["ypar"] = y_partitions
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cells_per_partition
pfw["nJ"] = pfw["ypar"] * cells_per_partition
pfw["nK"] = 3
pfw["ppcx"] = 2
pfw["ppcy"] = 5

cell_size = domain_height / max(pfw["nJ"], 1)
pfw["zmin"] = -0.5 * cell_size
pfw["zmax"] =  0.5 * cell_size

# --------------------------------------------------------------------------------------
# Batch defaults.  The verification harness may override these with command-line options.
# --------------------------------------------------------------------------------------
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = True
pfw["autoRestart"] = False
pfw["lastRestartBufferInSeconds"] = 300

# --------------------------------------------------------------------------------------
# MPM solver options.
# --------------------------------------------------------------------------------------
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 30.0
pfw["restartInterval"] = 2.0 * stop_time

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-16
pfw["cpdiDomainScaling"] = 1
pfw["damageFieldPartitioning"] = 1
# Default-on solver guard for this severe deformation case: if a particle stencil
# maps outside the local grid, flag that particle for deletion and compact the
# active-ordinal mapping rows rather than letting a stale/bad row reach DFG/P2G/G2P.
pfw["flagParticlesWithBadMappingArraysForDeletion"] = 1

pfw["solverProfiling"] = 0
pfw["needsNeighborList"] = 0
pfw["reactionHistory"] = 0
pfw["boxAverageHistory"] = 0
pfw["useEvents"] = 0
pfw["frictionCoefficient"] = 0.0

pfw["plotUnscaledParticles"] = 1
pfw["maxParticleVelocity"] = 10.0
pfw["minParticleJacobian"] = 0.01
pfw["maxParticleJacobian"] = 10.0

pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 2
pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"] = 1
pfw["explicitSurfaceNormalInfluence"] = 53.2702603781035
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"] = 1

# Current GEOS expects the writeStatistics enum string, not a legacy 0/1 flag.
pfw["writeStatistics"] = "all"

# Momentum verification diagnostics.  The solver writes mpm_momentumHistory.csv.
# Use logMomentum=2 here because this test is intended to diagnose where momentum
# first drifts between P2G, grid update/contact, BCs, FMPM, and G2P.
pfw["logMomentum"] = 2
pfw["logStartCycle"] = 0

pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = [
    "gridMass",
    "gridVelocity",
    "gridMomentum",
    "gridInternalForce",
    "gridContactForce",
]

pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "Damage",
    "StrengthScale",
    "SurfaceFlag",
    "RVector",
    "SurfaceNormal",
    "SurfacePosition",
]

# --------------------------------------------------------------------------------------
# Prescribed compaction.  The x boundaries are periodic; the top and bottom y boundaries
# move symmetrically so the prescribed motion does not introduce net x-momentum.
# --------------------------------------------------------------------------------------
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [0, 0, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [
    [0.0,       1.0, 1.0,             1.0],
    [stop_time, 1.0, final_y_stretch, 1.0],
]

# --------------------------------------------------------------------------------------
# Material properties from pfw_materials.py.
# --------------------------------------------------------------------------------------
material = getattr(material_db, material_name).copy()
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

# --------------------------------------------------------------------------------------
# Geometry.
# --------------------------------------------------------------------------------------
def make_objects():
    """Create a periodic array of circular inclusions using deterministic Poisson sampling."""
    disk_radius = 0.1
    disk_spacing = 1.05 * disk_radius
    x0 = [pfw["xmin"], pfw["ymin"]]
    dx = [pfw["xmax"] - pfw["xmin"], pfw["ymax"] - pfw["ymin"]]

    circles = geom.poisson(2.0 * disk_spacing, x0=x0, dx=dx, seed=seed,
                           numTrials=100, dim=2, periodic=[True, False])
    circles = geom.add_pore_images(circles, x0=x0, dx=dx, dim=2,
                                   periodic=[True, False])[0]

    objects = []
    for n, center in enumerate(circles):
        objects.append(
            geom.cylinder(
                "circle_{}".format(n),
                [center[0], center[1], 10.0 * pfw["zmin"]],
                [center[0], center[1], 10.0 * pfw["zmax"]],
                disk_radius,
                vel=[0.0, 0.0, 0.0],
                mat=0,
                group=0,
                particleType=particle_type,
            )
        )
    return objects


pfw["objects"] = make_objects()
print("Configured {} PBC compaction momentum subcase using material {}".format(variant_label, material_name))
