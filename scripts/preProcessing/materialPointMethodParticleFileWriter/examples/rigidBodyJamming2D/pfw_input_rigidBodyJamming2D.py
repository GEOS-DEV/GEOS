# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: MPI/periodic rigid-body MPM jamming of color-partitioned 2D copper objects followed by continuum compression.
# Solver/PFW features: ParticleColor, shared catch-all ContactGroup, two rigid fields, SinglePointBSpline particles,
#                      periodic x, moving y walls, exact analytic surface normals/positions, event-relative continuum loading.
# Workflow note: keep this file as a copyable problem definition. Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

"""Rigid-body MPM jamming example for copy-and-adapt PFW workflows.

The example contains six initially non-overlapping plane-strain rigid bodies:
three circular disks and three polygonal blocks. All bodies share one catch-all
``ContactGroup=0`` while each body has a distinct ``ParticleColor``. The rigid
event uses only two nodal fields: the deterministic first color at a node and a
weighted overflow field. Thus the grid-field allocation is independent of the
number of global colors and matches a later two-field DFG continuum solve.

The x direction is periodic and decomposed over two MPI partitions; y is split
over three partitions and compacted by opposing moving contact walls. One disk
straddles the x-periodic seam, represented by two PFW geometry images carrying
the same color, so the example exercises the unwrapped rank-zero body solve.

The initial cells satisfy ``DY = 1.5*DX``. PFW's directional particle counts are
``ppcx=2`` and ``ppcy=3``, giving equal x/y particle spacing and square particles
inside the rectangular cells. SinglePointBSpline particles avoid CPDI scaling
while retaining exact imported surface normals and surface-position vectors.

This is a normal ``pfw_input`` file and can be copied into the usual
``runClean_<user>.sh rigidBodyJamming2D`` workflow.
"""

#[pfw_dependency] pfw:pfw_materials.py

import math
import numpy as np

import pfw_geometryObjects as geom
import pfw_materials as matdb

pfw = {}
pfw["runDebug"] = True

# =============================================================================
# User-facing problem controls
# =============================================================================

rigid_compaction_time = 5.0e-5   # us
rigid_event_timeout = 1.5 * rigid_compaction_time
continuum_ramp_time = 2.5e-5     # us, measured from the actual end of rigid packing
stopTime = rigid_event_timeout + continuum_ramp_time
compressionDirection = "y"

# Open-loop y compaction is deliberately aggressive; maxJammingStress normally
# terminates the rigid event before the table reaches this final stretch.
rigid_y_stretch = 0.20

# Copper properties in the PFW material database use mm-us-mg-GPa units.
yield_stress = float(matdb.copper["defaultYieldStrength"])
continuum_target_stress = -2.0 * yield_stress

# Rigid-event safeguards. maxJammingStress has the same stress units as the
# constitutive model. maxContactPenetration is a dimensionless fraction of the
# current plane-strain grid-cell diagonal.
rigid_max_jamming_stress = 0.05 * yield_stress
rigid_max_contact_penetration = 0.05
rigid_contact_cfl = 0.05
rigid_max_time_step = 1.0e-7
rigid_penetration_penalty_beta = 0.02

# =============================================================================
# Domain, MPI partitioning, rectangular cells, and square particles
# =============================================================================

pfw["planeStrain"] = 1
pfw["xmin"] = -0.5
pfw["xmax"] =  0.5
pfw["ymin"] = -0.35
pfw["ymax"] =  0.35
pfw["periodic"] = [True, False, False]

# Exactly 2 x 3 x 1 MPI partitions.
pfw["xpar"] = 2
pfw["ypar"] = 3
pfw["zpar"] = 1

# nI/nJ count total cells, including the two y buffer cells because y is not
# periodic. There are 60 physical x cells and 28 physical y cells:
#   DX = 1.0/60
#   DY = 0.70/28 = 1.5*DX.
pfw["nI"] = 60
pfw["nJ"] = 30
pfw["nK"] = 3

# PFW calls these ppcx/ppcy (particles per cell per direction). With
# ppcx=2, ppcy=3 and DY=1.5*DX, the particle spacing is identical in x and y.
pfw["ppc"] = 2
pfw["ppcx"] = 2
pfw["ppcy"] = 3
pfw["ppcz"] = 1

physical_nx = pfw["nI"]
physical_ny = pfw["nJ"] - 2
dx_grid = (pfw["xmax"] - pfw["xmin"]) / physical_nx
dy_grid = (pfw["ymax"] - pfw["ymin"]) / physical_ny
particle_dx = dx_grid / pfw["ppcx"]
particle_dy = dy_grid / pfw["ppcy"]
assert math.isclose(dy_grid, 1.5 * dx_grid, rel_tol=0.0, abs_tol=1.0e-14)
assert math.isclose(particle_dx, particle_dy, rel_tol=0.0, abs_tol=1.0e-14)

# One physical plane-strain z cell, with a thickness equal to DX.
thickness = dx_grid
pfw["zmin"] = -0.5 * thickness
pfw["zmax"] =  0.5 * thickness
pfw["outputType"] = "silo"

# =============================================================================
# Batch settings used by particleFileWriter.py
# =============================================================================

pfw["mWallTime"] = "00:12:00"
pfw["mCores"] = pfw["xpar"] * pfw["ypar"] * pfw["zpar"]
pfw["mSubmitJobs"] = True

# =============================================================================
# GEOS MPM solver controls
# =============================================================================

pfw["endTime"] = stopTime
pfw["plotInterval"] = stopTime / 12.0
pfw["restartInterval"] = 2.0 * stopTime

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.35
pfw["initialDt"] = 1.0e-12
pfw["updateMethod"] = "PIC"
pfw["useEvents"] = 1
pfw["eventReporting"] = 1

# One shared ContactGroup and continuum DFG with two flags both need two fields;
# hundreds of ParticleColor values do not increase this allocation. During the
# rigid event the two fields are color slots; after the event they revert to the
# ordinary surface-flag/damage DFG partition used by continuum auto-contact.
pfw["rigidBodyMaxGridFields"] = 2
pfw["damageFieldPartitioning"] = 1

# Rigid and continuum contact both use the exact imported surface fields.
pfw["frictionCoefficient"] = 0.05
pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"] = 1
pfw["explicitSurfaceNormalInfluence"] = 25.0
pfw["computeParticleSurfaceNormalsAndPositions"] = 0
pfw["overwriteExistingNormalsAndPositions"] = 0
pfw["cpdiDomainScaling"] = 0
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"] = 0

pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime / 250.0
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stopTime / 500.0
pfw["maxParticleVelocity"] = 10.0
pfw["minParticleJacobian"] = 0.01
pfw["maxParticleJacobian"] = 10.0

pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "ParticleColor",
    "SurfaceFlag",
    "SurfaceNormal",
    "SurfacePosition",
    "RVector",
]

# =============================================================================
# Material
# =============================================================================

pfw["materials"] = [matdb.copper["name"]]
pfw["materialPropertyString"] = matdb.copper["materialString"]
COPPER = 0
SINGLE_POINT_BSPLINE = 1
CONTACT_GROUP = 0

# =============================================================================
# Geometry: six bodies (one periodic disk is represented by two image pieces)
# =============================================================================

def color_object(obj, color):
    """Attach the PFW ParticleColor written by particleFileWriter.py."""
    obj.particleColor = int(color)
    return obj


def use_plane_strain_disk_normals(obj, center, radius):
    """Override through-thickness cylinder normals with exact 2D disk normals."""
    cx, cy = center
    r0 = float(radius)

    def _normal(pt):
        vec = np.array([float(pt[0]) - cx, float(pt[1]) - cy, 0.0])
        mag = np.linalg.norm(vec[:2])
        if mag <= 0.0:
            return np.array([1.0, 0.0, 0.0])
        return vec / mag

    def _surface_position(pt):
        n = _normal(pt)
        radial_distance = np.linalg.norm([float(pt[0]) - cx, float(pt[1]) - cy])
        return (r0 - radial_distance) * n

    obj.getSurfaceNormal = _normal
    obj.getSurfacePosition = _surface_position
    return obj


def rectangle_vertices(center, width, height, angle_degrees=0.0):
    cx, cy = center
    hw = 0.5 * width
    hh = 0.5 * height
    verts = np.array([[-hw, -hh], [hw, -hh], [hw, hh], [-hw, hh]], dtype=float)
    theta = math.radians(angle_degrees)
    rot = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    return (verts @ rot.T) + np.array([cx, cy])


objects = []

objects.append(color_object(use_plane_strain_disk_normals(geom.cylinder(
    "disk_color_0",
    [-0.28, -0.20, pfw["zmin"]],
    [-0.28, -0.20, pfw["zmax"]],
    0.105,
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP,
    particleType=SINGLE_POINT_BSPLINE,
    flaggedSurfaces=[False, False, True, False]), (-0.28, -0.20), 0.105), 0))

objects.append(color_object(geom.polygon(
    "block_color_1",
    rectangle_vertices(center=(0.02, -0.19), width=0.18, height=0.115, angle_degrees=8.0),
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP,
    particleType=SINGLE_POINT_BSPLINE), 1))

# Color 2 straddles x=+/-0.5. The two image cylinders are one rigid body because
# they share ParticleColor=2; their exact normals use their respective image centers.
for image_name, image_center_x in (("right", 0.47), ("left_image", -0.53)):
    objects.append(color_object(use_plane_strain_disk_normals(geom.cylinder(
        f"disk_color_2_{image_name}",
        [image_center_x, -0.18, pfw["zmin"]],
        [image_center_x, -0.18, pfw["zmax"]],
        0.075,
        vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP,
        particleType=SINGLE_POINT_BSPLINE,
        flaggedSurfaces=[False, False, True, False]),
        (image_center_x, -0.18), 0.075), 2))

objects.append(color_object(geom.polygon(
    "block_color_3",
    rectangle_vertices(center=(-0.28, 0.14), width=0.16, height=0.155, angle_degrees=-12.0),
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP,
    particleType=SINGLE_POINT_BSPLINE), 3))

objects.append(color_object(use_plane_strain_disk_normals(geom.cylinder(
    "disk_color_4",
    [0.03, 0.12, pfw["zmin"]],
    [0.03, 0.12, pfw["zmax"]],
    0.115,
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP,
    particleType=SINGLE_POINT_BSPLINE,
    flaggedSurfaces=[False, False, True, False]), (0.03, 0.12), 0.115), 4))

objects.append(color_object(geom.polygon(
    "block_color_5",
    rectangle_vertices(center=(0.34, 0.14), width=0.145, height=0.135, angle_degrees=15.0),
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP,
    particleType=SINGLE_POINT_BSPLINE), 5))

pfw["objects"] = objects

# =============================================================================
# Boundary conditions and staged loading
# =============================================================================

# Boundary order is x-, x+, y-, y+, z-, z+. The x pair is periodic, so no
# physical x-face BC is applied. Opposing y Contact walls move with the F-table;
# z is symmetry for plane strain.
pfw["boundaryConditionTypes"] = [0, 0, 3, 3, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["prescribedFTable"] = 0
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [
    [0.0,                    1.0, 1.0,             1.0],
    [rigid_compaction_time,  1.0, rigid_y_stretch, 1.0],
]

pfw["mpmEventsString"] = f"""
<DeformationUpdate
  name="startRigidCompression"
  startTime="0.0"
  prescribedBoundaryFTable="1"
  prescribedFTable="0"
  stressControl="{{0, 0, 0}}"
  fTableInterpType="Cosine"
  fTable="{{{{0.0, 1.0, 1.0, 1.0}}, {{{rigid_compaction_time}, 1.0, {rigid_y_stretch}, 1.0}}}}"/>

<RigidBodyMPM
  name="rigidJamming"
  startTime="0.0"
  endTime="{rigid_event_timeout}"
  minActiveTime="{0.25 * rigid_compaction_time}"
  maxGridFields="2"
  linearDamping="2.0"
  angularDamping="2.0"
  contactCFL="{rigid_contact_cfl}"
  maxTimeStep="{rigid_max_time_step}"
  rigidBodyPenetrationPenaltyBeta="{rigid_penetration_penalty_beta}"
  maxJammingStress="{rigid_max_jamming_stress}"
  maxContactPenetration="{rigid_max_contact_penetration}"
  writeHistory="1"
  historyWriteInterval="{stopTime / 500.0}"/>

<DeformationUpdate
  name="continuumCompression"
  dependencies="{{ rigidJamming }}"
  delay="0.0"
  duration="{continuum_ramp_time}"
  relativeDeformation="1"
  prescribedBoundaryFTable="0"
  prescribedFTable="0"
  stressControl="{{0, 1, 0}}"
  stressTableInterpType="Cosine"
  stressTable="{{{{0.0, 0.0, 0.0, 0.0}}, {{{continuum_ramp_time}, 0.0, {continuum_target_stress}, 0.0}}}}"/>
"""
