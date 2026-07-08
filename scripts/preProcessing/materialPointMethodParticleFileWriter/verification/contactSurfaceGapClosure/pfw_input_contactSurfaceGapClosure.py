# -*- coding: utf-8 -*-
"""PFW input for the 2D contact surface/gap-closure verification test.

The folder runTest script dispatches this one human-readable input for an array
of contact variants. Environment variables are used only for the small set of
values that define each subcase.  Each templated value is documented with the
concrete replacements used by the run script so a user can copy this file and
replace the defaults to make a standalone input.
"""

import math
import os

import numpy as np

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}
pfw["runDebug"] = True

# --------------------------------------------------------------------------------------
# Variant template values set by runTest.
# --------------------------------------------------------------------------------------
# CONTACT_GAP_CASE_NAME replacement examples:
#   Multi-field, explicit surfaces, softened gap correction:
#       "contactGapClosure_multiField_explicit_Softened"
#   DFG, implicit surfaces, implicit gap correction:
#       "contactGapClosure_dfg_implicit_Implicit"
case_name = os.environ.get("CONTACT_GAP_CASE_NAME", "contactGapClosure_multiField_explicit_Softened")

# CONTACT_GAP_VARIANT_LABEL replacement examples:
#   "Multi-field / explicit surfaces / Softened gap"
#   "DFG / implicit surfaces / Implicit gap"
variant_label = os.environ.get("CONTACT_GAP_VARIANT_LABEL", "Multi-field / explicit surfaces / Softened gap")

# CONTACT_GAP_FIELD_MODE replacement examples:
#   "multiField" -> separate contact groups define the two velocity fields
#   "dfg"        -> damage-field partitioning defines separated fields at the interface
field_mode = os.environ.get("CONTACT_GAP_FIELD_MODE", "multiField")

# CONTACT_GAP_SURFACE_MODE replacement examples:
#   "implicit" -> let the solver infer surface normals/positions from the particle fields
#   "explicit" -> write particle SurfaceNormal and SurfacePosition from the geometry below
surface_mode = os.environ.get("CONTACT_GAP_SURFACE_MODE", "explicit")

# CONTACT_GAP_CORRECTION replacement examples:
#   "Simple", "Implicit", or "Softened"
contact_gap_correction = os.environ.get("CONTACT_GAP_CORRECTION", "Softened")

# CONTACT_GAP_PARTICLE_TYPE replacement examples:
#   suite default: "2" -> CPDI particles
#   optional:       "1" -> SinglePointBspline particles
particle_type = int(os.environ.get("CONTACT_GAP_PARTICLE_TYPE", "2"))

# CONTACT_GAP_MATERIAL replacement examples:
#   suite default: "contactGapClosureHyperelastic"
material_name = os.environ.get("CONTACT_GAP_MATERIAL", "contactGapClosureHyperelastic")

# --------------------------------------------------------------------------------------
# Geometry and loading parameters.
# --------------------------------------------------------------------------------------
domain_width = float(os.environ.get("CONTACT_GAP_DOMAIN_WIDTH", "1.0"))
domain_height = float(os.environ.get("CONTACT_GAP_DOMAIN_HEIGHT", "1.0"))
initial_gap = float(os.environ.get("CONTACT_GAP_INITIAL_GAP", "0.10"))
extra_compression = float(os.environ.get("CONTACT_GAP_EXTRA_COMPRESSION", "0.12"))
notch_radius = float(os.environ.get("CONTACT_GAP_RADIUS", "0.18"))

if initial_gap <= 0.0:
    raise ValueError("CONTACT_GAP_INITIAL_GAP must be positive")
if extra_compression <= 0.0:
    raise ValueError("CONTACT_GAP_EXTRA_COMPRESSION must be positive")
if 2.0 * notch_radius >= domain_width:
    raise ValueError("CONTACT_GAP_RADIUS must be less than half the domain width")
if initial_gap + extra_compression >= 0.5 * domain_height:
    raise ValueError("CONTACT_GAP_INITIAL_GAP + CONTACT_GAP_EXTRA_COMPRESSION is too large for this compact test")

material = getattr(material_db, material_name).copy()
wave_speed = float(material.get("waveSpeed", 1.0))
# Slow-loading quasistatic target: 100 elastic-wave transit times across the solid height.
transit_time = domain_height / max(wave_speed, 1.0e-30)
stop_time = float(os.environ.get("CONTACT_GAP_STOP_TIME", str(100.0 * transit_time)))

# CONTACT_GAP_REFINE / CONTACT_GAP_CPP / CONTACT_GAP_PPC replacement examples:
#   suite default: refine=1, cpp=32, ppc=2 -> 32 x 32 cells on one rank
#   diagnostic:    refine=1, cpp=64, ppc=2 -> higher-resolution one-rank run
# CONTACT_GAP_GRID_CELLS_PER_PARTITION is retained as a backward-compatible alias for cpp.
refine = int(os.environ.get("CONTACT_GAP_REFINE", "1"))
cpp = int(os.environ.get("CONTACT_GAP_CPP", os.environ.get("CONTACT_GAP_GRID_CELLS_PER_PARTITION", "32")))
ppc = int(os.environ.get("CONTACT_GAP_PPC", "2"))
x_partitions = int(os.environ.get("CONTACT_GAP_X_PARTITIONS", str(refine)))
y_partitions = int(os.environ.get("CONTACT_GAP_Y_PARTITIONS", str(refine)))
if refine < 1:
    raise ValueError("CONTACT_GAP_REFINE must be positive")
if cpp < 8:
    raise ValueError("CONTACT_GAP_CPP must be at least 8")
if ppc < 1:
    raise ValueError("CONTACT_GAP_PPC must be positive")

# --------------------------------------------------------------------------------------
# Domain and resolution.
# --------------------------------------------------------------------------------------
pfw["caseName"] = case_name
pfw["xmin"] = -0.5 * domain_width
pfw["xmax"] =  0.5 * domain_width
pfw["ymin"] = -0.5 * domain_height
pfw["ymax"] =  0.5 * domain_height
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]

pfw["xpar"] = x_partitions
pfw["ypar"] = y_partitions
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = 3
pfw["ppcx"] = ppc
pfw["ppcy"] = ppc

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
pfw["plotInterval"] = stop_time / 8.0
pfw["restartInterval"] = 2.0 * stop_time

pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.5
pfw["initialDt"] = 1.0e-16
pfw["cpdiDomainScaling"] = 1
pfw["damageFieldPartitioning"] = 1 if field_mode == "dfg" else 0

pfw["solverProfiling"] = 0
pfw["needsNeighborList"] = 0
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 250.0
pfw["boxAverageHistory"] = 0
pfw["useEvents"] = 0
pfw["frictionCoefficient"] = 0.0
pfw["maxParticleVelocity"] = 10.0
pfw["minParticleJacobian"] = 0.01
pfw["maxParticleJacobian"] = 10.0
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 2
pfw["contactGapCorrection"] = contact_gap_correction
pfw["useSurfacePositionForContact"] = 1 if surface_mode == "explicit" else 0
pfw["explicitSurfaceNormalInfluence"] = 47.14045207910317 if surface_mode == "explicit" else 0
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"] = 0 if surface_mode == "explicit" else 1
pfw["useInternalForceAsFaceReaction"] = 1
pfw["writeStatistics"] = "all"

pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = [
    "gridMass",
    "gridVelocity",
    "gridInternalForce",
    "gridContactForce",
]

pfw["plottableFields"] = [
    "particleID",
    "particleMass",
    "particleVolume",
    "particleDensity",
    "particleMaterialType",
    "particleGroup",
    "particleSurfaceFlag",
    "particleDamage",
    "particleCenter",
    "particleReferencePosition",
    "particleVelocity",
    "particleStress",
    "particleSurfaceNormal",
    "particleSurfacePosition",
    "particleRVectors",
    "particleReferenceRVectors",
    "gridMass",
    "gridVelocity",
]

pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "Damage",
    "SurfaceFlag",
    "RVector",
]
if surface_mode == "explicit":
    pfw["particleFileFields"].extend(["SurfaceNormal", "SurfacePosition"])

# --------------------------------------------------------------------------------------
# Prescribed quasistatic gap closure and compression.
# --------------------------------------------------------------------------------------
final_y_stretch = 1.0 - (initial_gap + extra_compression) / domain_height
closure_y_stretch = 1.0 - initial_gap / domain_height
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [
    [0.0,       1.0, 1.0,             1.0],
    [stop_time, 1.0, final_y_stretch, 1.0],
]

# --------------------------------------------------------------------------------------
# Material properties from pfw_materials.py.
# --------------------------------------------------------------------------------------
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

# --------------------------------------------------------------------------------------
# Geometry.
# --------------------------------------------------------------------------------------
class CurvedGapBlock(geom.Geometry):
    """Two-dimensional block with one sinus-free circular contact surface.

    The upper block has a downward convex semicircular lobe.  The lower block
    has the matching upward-facing concave recess.  The two contact profiles are
    vertical translations of each other by ``initial_gap``, so idealized rigid
    motion closes the entire curved gap at one displacement.
    """

    def __init__(self, name, which, group, mat=0, particleType=2):
        super().__init__(name, vel=[0.0, 0.0, 0.0], mat=mat, group=group, particleType=particleType, damage=0.0)
        if which not in ("lower", "upper"):
            raise ValueError("which must be 'lower' or 'upper'")
        self.which = which
        self.xmin = pfw["xmin"]
        self.xmax = pfw["xmax"]
        self.ymin = pfw["ymin"]
        self.ymax = pfw["ymax"]
        self.radius = notch_radius
        self.gap = initial_gap

    def _flat_y(self):
        return -0.5 * self.gap if self.which == "lower" else 0.5 * self.gap

    def _circle_center(self):
        return np.array([0.0, self._flat_y(), 0.0], dtype=float)

    def _profile_y(self, x):
        ax = abs(float(x))
        if ax >= self.radius:
            sag = 0.0
        else:
            sag = math.sqrt(max(self.radius * self.radius - ax * ax, 0.0))
        return self._flat_y() - sag

    def _profile_dydx(self, x):
        ax = abs(float(x))
        if ax <= 1.0e-12 or ax >= 0.999 * self.radius:
            return 0.0
        return float(x) / math.sqrt(max(self.radius * self.radius - x * x, 1.0e-30))

    def _profile_normal(self, x):
        # Outward normals point from the solid into the open gap.  The lower
        # block recess uses the opposite of the circle radial direction, while
        # the upper lobe uses the radial direction.
        dydx = self._profile_dydx(x)
        if self.which == "lower":
            n = np.array([-dydx, 1.0, 0.0], dtype=float)
        else:
            n = np.array([dydx, -1.0, 0.0], dtype=float)
        norm = np.linalg.norm(n)
        return n / norm if norm > 0.0 else np.array([0.0, 1.0 if self.which == "lower" else -1.0, 0.0])

    def _arc_surface_candidate(self, pt):
        """Closest vector and normal for the circular mating surface.

        The explicit surface-position field stores a vector from the particle to
        its nearest surface.  A vertical projection to y=profile(x) is visibly
        wrong on the round interface, so project radially to the circular arc.
        """
        p = np.array([float(pt[0]), float(pt[1]), 0.0], dtype=float)
        center = self._circle_center()
        radial = p - center
        radial[2] = 0.0
        rnorm = np.linalg.norm(radial[:2])
        if rnorm <= 1.0e-14:
            unit = np.array([0.0, -1.0, 0.0], dtype=float)
        else:
            unit = radial / rnorm
            # The mating interface is the lower semicircle.  If a far interior
            # particle would project to the upper semicircle, clamp to the
            # closest arc endpoint.  Particles close enough to be flagged as
            # surface particles project onto the true curved face.
            if unit[1] > 0.0:
                unit = np.array([1.0 if unit[0] >= 0.0 else -1.0, 0.0, 0.0], dtype=float)

        q = center + self.radius * unit
        delta = q - p
        arc_radial = q - center
        arc_norm = np.linalg.norm(arc_radial[:2])
        if arc_norm > 1.0e-14:
            arc_radial = arc_radial / arc_norm
        else:
            arc_radial = np.array([0.0, -1.0, 0.0], dtype=float)
        normal = -arc_radial if self.which == "lower" else arc_radial
        return (np.linalg.norm(delta[:2]), delta, normal)

    def _inside_xy(self, x, y):
        if x < self.xmin or x >= self.xmax:
            return False
        profile = self._profile_y(x)
        if self.which == "lower":
            return self.ymin <= y < profile
        return profile <= y < self.ymax

    def isInterior(self, pt, skinDepth):
        x = float(pt[0])
        y = float(pt[1])
        if not self._inside_xy(x, y):
            return -1
        if self._distance_to_nearest_surface(pt)[0] <= skinDepth:
            return 2
        return 0

    def _surface_candidates(self, pt):
        x = float(pt[0])
        y = float(pt[1])
        candidates = []
        candidates.append((abs(x - self.xmin), np.array([self.xmin - x, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])))
        candidates.append((abs(self.xmax - x), np.array([self.xmax - x, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])))
        if self.which == "lower":
            candidates.append((abs(y - self.ymin), np.array([0.0, self.ymin - y, 0.0]), np.array([0.0, -1.0, 0.0])))
            if abs(x) < self.radius:
                candidates.append(self._arc_surface_candidate(pt))
            else:
                py = self._profile_y(x)
                candidates.append((abs(py - y), np.array([0.0, py - y, 0.0]), self._profile_normal(x)))
        else:
            candidates.append((abs(self.ymax - y), np.array([0.0, self.ymax - y, 0.0]), np.array([0.0, 1.0, 0.0])))
            if abs(x) < self.radius:
                candidates.append(self._arc_surface_candidate(pt))
            else:
                py = self._profile_y(x)
                candidates.append((abs(y - py), np.array([0.0, py - y, 0.0]), self._profile_normal(x)))
        return candidates

    def _distance_to_nearest_surface(self, pt):
        return min(self._surface_candidates(pt), key=lambda item: item[0])

    def getSurfaceNormal(self, pt):
        return self._distance_to_nearest_surface(pt)[2]

    def getSurfacePosition(self, pt):
        return self._distance_to_nearest_surface(pt)[1]

    def xMin(self):
        return self.xmin

    def xMax(self):
        return self.xmax

    def getSubregions(self):
        return [(self.mat, self.particleType)]


if field_mode == "multiField":
    lower_group = 0
    upper_group = 1
elif field_mode == "dfg":
    lower_group = 0
    upper_group = 0
else:
    raise ValueError("CONTACT_GAP_FIELD_MODE must be 'multiField' or 'dfg'")

if surface_mode not in ("implicit", "explicit"):
    raise ValueError("CONTACT_GAP_SURFACE_MODE must be 'implicit' or 'explicit'")

pfw["objects"] = [
    CurvedGapBlock("lowerBlock", "lower", group=lower_group, mat=0, particleType=particle_type),
    CurvedGapBlock("upperBlock", "upper", group=upper_group, mat=0, particleType=particle_type),
]

print(
    "Configured {label}: field_mode={field}, surface_mode={surface}, gap={gap}, "
    "closure_Fyy={closure:.6g}, final_Fyy={final:.6g}, stop_time={stop:.6g}".format(
        label=variant_label,
        field=field_mode,
        surface=surface_mode,
        gap=contact_gap_correction,
        closure=closure_y_stretch,
        final=final_y_stretch,
        stop=stop_time,
    )
)
