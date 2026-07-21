# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Rigid-body MPM jamming of color-partitioned 2D copper objects followed by continuum compression.
# Solver/PFW features: ParticleColor, analytic surface normals/positions, RigidBodyMPM event, F-table compaction, stress-control continuation, reaction history.
# Workflow note: keep this file as a copyable problem definition.  Run directories,
# Slurm submission, rerun cleanup, and suite reporting belong in runProblem or
# examples/mpm_example_runner.py, not in the pfw dictionary below.
# ---- end example input metadata ----

"""Rigid-body MPM jamming example for copy-and-adapt PFW workflows.

The example contains six initially non-overlapping plane-strain objects: three
circular disks and three polygonal blocks.  All objects use the same copper
material and contact group, but each object sets a distinct ``ParticleColor``.
During the ``RigidBodyMPM`` event the solver forms one rigid body per color and
uses the color-partitioned nodal fields.  The final/overflow field uses the
weighted-lumping policy when more colors meet at a node than ``maxGridFields``.

The geometric objects provide exact analytic surface normals and surface-position
vectors.  PFW writes those fields to the particle file and the solver is told not
to overwrite them.

This file is intentionally a normal ``pfw_input``: it can be copied into any PFW
run location and launched with the usual ``runClean_<user>.sh rigidBodyJamming2D``
workflow.
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

refine = 1
cells_per_partition = 14
particles_per_cell = 2

rigid_compaction_time = 5.0e-5   # us
continuum_ramp_time = 2.5e-5     # us
stopTime = rigid_compaction_time + continuum_ramp_time
compressionDirection = "y"

# F-table compaction in the rigid stage.  The x faces remain fixed contact
# walls; the y contact faces move inward to jam the bodies.  The final stretch
# is intentionally small enough for the moving walls to engage the object pack.
rigid_y_stretch = 0.20

# Copper properties in the PFW material database use mm-us-mg-GPa units.
yield_stress = float(matdb.copper["defaultYieldStrength"])
continuum_target_stress = -2.0 * yield_stress

# =============================================================================
# Domain and grid
# =============================================================================

pfw["planeStrain"] = 1
pfw["xmin"] = -0.5
pfw["xmax"] =  0.5
pfw["ymin"] = -0.35
pfw["ymax"] =  0.35
pfw["periodic"] = [False, False, False]

pfw["xpar"] = 4 * refine
pfw["ypar"] = 4 * refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cells_per_partition
pfw["nJ"] = pfw["ypar"] * cells_per_partition
pfw["nK"] = 3
pfw["ppc"] = particles_per_cell

# Choose the plane-strain thickness so the single physical z cell is comparable
# with the in-plane grid spacing.
dx_grid = (pfw["xmax"] - pfw["xmin"]) / max(pfw["nI"] - 2, 1)
dy_grid = (pfw["ymax"] - pfw["ymin"]) / max(pfw["nJ"] - 2, 1)
thickness = min(dx_grid, dy_grid)
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
# The runProblem wrapper preserves this interval for this case so the VisIt
# renderer has enough states to show color at four times.
pfw["plotInterval"] = stopTime / 12.0
pfw["restartInterval"] = 2.0 * stopTime

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.35
pfw["initialDt"] = 1.0e-12
pfw["updateMethod"] = "PIC"
pfw["useEvents"] = 1
pfw["eventReporting"] = 1

# The RigidBodyMPM event also reserves maxGridFields during solver
# initialization.  The solver-level attribute is retained here as documentation
# and for compatibility with older development branches.
pfw["rigidBodyMaxGridFields"] = 4

# Rigid and continuum contact both use the explicit analytic surface fields from
# the particle file.
pfw["frictionCoefficient"] = 0.05
pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"] = 1
pfw["explicitSurfaceNormalInfluence"] = 25.0
pfw["computeParticleSurfaceNormalsAndPositions"] = 0
pfw["overwriteExistingNormalsAndPositions"] = 0

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

# matdb.copper is a constant-yield J2 copper parameterization.  In this example
# it is used as the elastic-perfectly-plastic copper response for the continuum
# compression after the rigid-body jamming stage.
pfw["materials"] = [matdb.copper["name"]]
pfw["materialPropertyString"] = matdb.copper["materialString"]
COPPER = 0
CPDI = 2
CONTACT_GROUP = 0

# =============================================================================
# Geometry: three disks and three blocks, initially non-overlapping
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

# Disks are represented as through-thickness cylinders.  Only the cylindrical
# wall is flagged as a contact surface in plane strain.
objects.append(color_object(use_plane_strain_disk_normals(geom.cylinder(
    "disk_color_0",
    [-0.34, -0.20, pfw["zmin"]],
    [-0.34, -0.20, pfw["zmax"]],
    0.105,
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP, particleType=CPDI,
    flaggedSurfaces=[False, False, True, False]), (-0.34, -0.20), 0.105), 0))

objects.append(color_object(geom.polygon(
    "block_color_1",
    rectangle_vertices(center=(-0.06, -0.19), width=0.18, height=0.115, angle_degrees=8.0),
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP, particleType=CPDI), 1))

objects.append(color_object(use_plane_strain_disk_normals(geom.cylinder(
    "disk_color_2",
    [0.24, -0.18, pfw["zmin"]],
    [0.24, -0.18, pfw["zmax"]],
    0.090,
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP, particleType=CPDI,
    flaggedSurfaces=[False, False, True, False]), (0.24, -0.18), 0.090), 2))

objects.append(color_object(geom.polygon(
    "block_color_3",
    rectangle_vertices(center=(-0.26, 0.14), width=0.16, height=0.155, angle_degrees=-12.0),
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP, particleType=CPDI), 3))

objects.append(color_object(use_plane_strain_disk_normals(geom.cylinder(
    "disk_color_4",
    [0.05, 0.12, pfw["zmin"]],
    [0.05, 0.12, pfw["zmax"]],
    0.115,
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP, particleType=CPDI,
    flaggedSurfaces=[False, False, True, False]), (0.05, 0.12), 0.115), 4))

objects.append(color_object(geom.polygon(
    "block_color_5",
    rectangle_vertices(center=(0.34, 0.14), width=0.145, height=0.135, angle_degrees=15.0),
    vel=[0.0, 0.0, 0.0], mat=COPPER, group=CONTACT_GROUP, particleType=CPDI), 5))

pfw["objects"] = objects

# =============================================================================
# Boundary conditions and staged loading
# =============================================================================

# Boundary-condition order is x-, x+, y-, y+, z-, z+.
# x faces are fixed contact container walls, y faces are contact walls driven
# by the F-table and then by y-stress control, and z faces are plane-strain
# symmetry faces.  Use Contact rather than Moving for the in-plane walls so the
# walls exchange impulses with free rigid bodies instead of merely prescribing
# grid velocities at already occupied boundary nodes.
pfw["boundaryConditionTypes"] = [3, 3, 3, 3, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["prescribedFTable"] = 0
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [
    [0.0,                    1.0,           1.0,            1.0],
    [rigid_compaction_time,  1.0,           rigid_y_stretch, 1.0],
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
  endTime="{rigid_compaction_time}"
  maxGridFields="2"
  linearDamping="2.0"
  angularDamping="2.0"
  maxForce="5.0e30"
  writeHistory="1"
  historyWriteInterval="{stopTime / 500.0}"/>

<DeformationUpdate
  name="continuumCompression"
  dependencies="{{ rigidJamming }}"
  delay="0.0"
  duration="0.0"
  relativeDeformation="1"
  prescribedBoundaryFTable="0"
  prescribedFTable="0"
  stressControl="{{0, 1, 0}}"
  stressTableInterpType="Cosine"
  stressTable="{{{{0.0, 0.0, 0.0, 0.0}}, {{{continuum_ramp_time}, 0.0, {continuum_target_stress}, 0.0}}}}"/>
"""
