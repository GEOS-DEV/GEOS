# -*- coding: utf-8 -*-
"""MMS vortex verification input template.

The folder-level runTest script runs this single input twice: once as a CPDI-FLIP
case and once as an FMPM/singlePointBspline case.  The run script selects the
variant by setting the environment variables documented below before invoking
particleFileWriter.py through mpm_vv_case_runner.py.

The comments beside each templated value show the concrete values used by
runTest.  A user can copy this file, replace the few os.environ.get defaults
with one of the commented concrete values, and obtain an ordinary single-case
PFW input file.
"""

import math
import os

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db
import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

pfw = {}

# Variant template values -------------------------------------------------------
# MMS_VARIANT_LABEL replacement examples:
#   CPDI-FLIP case:                "CPDI-FLIP"
#   FMPM/singlePointBspline case:  "FMPM/singlePointBspline"
variant_label = os.environ.get("MMS_VARIANT_LABEL", "CPDI-FLIP")

# MMS_CASE_NAME replacement examples:
#   CPDI-FLIP case:                "mmsVortex_cpdiFlip"
#   FMPM/singlePointBspline case:  "mmsVortex_fmpmSinglePointBspline"
case_name = os.environ.get("MMS_CASE_NAME", "mmsVortex_cpdiFlip")

# MMS_UPDATE_METHOD replacement examples:
#   CPDI-FLIP case:                "FLIP"
#   FMPM/singlePointBspline case:  "FMPM"
update_method = os.environ.get("MMS_UPDATE_METHOD", "FLIP")

# MMS_UPDATE_ORDER replacement examples:
#   CPDI-FLIP case:                ""   -> no pfw["updateOrder"] entry
#   FMPM/singlePointBspline case:  "2"  -> pfw["updateOrder"] = 2
update_order_text = os.environ.get("MMS_UPDATE_ORDER", "")

# MMS_PARTICLE_TYPE replacement examples:
#   CPDI-FLIP case:                "2"  -> CPDI particles
#   FMPM/singlePointBspline case:  "1"  -> SinglePointBspline particles
particle_type = int(os.environ.get("MMS_PARTICLE_TYPE", "2"))

# Run identity and batch policy -------------------------------------------------
pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

# Domain and discretization -----------------------------------------------------
# The MMS annulus is 0.75 < R < 1.25 inside a square block.  The verification
# resolution is large enough for the rendered displacement fields to be
# interpretable, while keeping each variant close to the target suite budget.
stopTime = 1.0
pfw["planeStrain"] = 1
pfw["xpar"] = 2
pfw["ypar"] = 2
pfw["zpar"] = 1
pfw["nI"] = 48
pfw["nJ"] = 48
pfw["nK"] = 3
pfw["ppc"] = 2

pfw["xmin"] = -2.0
pfw["xmax"] = 2.0
pfw["ymin"] = -2.0
pfw["ymax"] = 2.0

dx = (pfw["xmax"] - pfw["xmin"]) / float(pfw["nI"] - 2)
dy = (pfw["ymax"] - pfw["ymin"]) / float(pfw["nJ"] - 2)
thickness = 3.0 * min(dx, dy)
pfw["zmin"] = -0.5 * thickness
pfw["zmax"] = 0.5 * thickness

# Solver controls ---------------------------------------------------------------
pfw["endTime"] = stopTime
pfw["plotInterval"] = 0.5 * stopTime
pfw["restartInterval"] = 2.0 * stopTime
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.05
pfw["initialDt"] = 1.0e-12
pfw["maxParticleVelocity"] = 50.0

pfw["updateMethod"] = update_method
if update_order_text.strip():
    pfw["updateOrder"] = int(update_order_text)
pfw["cpdiDomainScaling"] = 0
pfw["bodyForce"] = [0.0, 0.0, 0.0]
pfw["generalizedVortexMMS"] = 1
pfw["reactionHistory"] = 0
pfw["boxAverageHistory"] = 0
# Current GEOS expects the writeStatistics enum string, not a legacy 0/1 flag.
# Allowed values are "none", "iteration", "convergence", and "all".
# The suite report gets scheduler wall time from job metadata; solver statistics
# are also enabled here so detailed GEOS-side timing/iteration diagnostics are
# available in run directories when a case needs deeper inspection.
pfw["writeStatistics"] = "all"

# The output is Silo because this verification post-processor renders VisIt
# images.  Background-grid fields are requested explicitly with current PFW
# syntax; no legacy siloGridFields aliases are used.
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity", "gridDisplacement"]

# Free boundaries: the manufactured body force drives the motion.
pfw["boundaryConditionTypes"] = [0, 0, 0, 0, 0, 0]
pfw["prescribedBcTable"] = 0
pfw["prescribedBoundaryFTable"] = 0

# Material ----------------------------------------------------------------------
# The generalized-vortex material card is centralized in pfw_materials.py.  This
# keeps this verification input focused on the test setup and makes it easy for
# users to reuse or edit the material parameterization in one shared location.
# The current pfw_materials.hyperelasticMMS entry corresponds to the Kamojjala
# et al. generalized-vortex MMS parameters: density=1000, lambda=577, and
# shear modulus=384.6153846153846 in the unit system used by this input.
mms_material = material_db.hyperelasticMMS.copy()
pfw["materials"] = [mms_material["name"]]
pfw["materialPropertyString"] = mms_material["materialString"]

# Particles ---------------------------------------------------------------------
# particle_type is templated above.  The runTest replacements are:
#   CPDI-FLIP case:                particleType=2
#   FMPM/singlePointBspline case:  particleType=1
block = geom.box(
    "mms_vortex_block",
    [pfw["xmin"] + 2.0 * dx, pfw["ymin"] + 2.0 * dy, pfw["zmin"]],
    [pfw["xmax"] - 2.0 * dx, pfw["ymax"] - 2.0 * dy, pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
    particleType=particle_type,
)
pfw["objects"] = [block]

# Lagrangian tracer metric ------------------------------------------------------
# Tracers are seeded in the annulus where the manufactured rotation is nonzero.
# The post-processor evaluates the average norm of their displacement history.
def annular_tracer_points():
    points = []
    for ring in range(5):
        radius = 0.80 + 0.10 * ring
        for point in range(32):
            theta = 2.0 * math.pi * point / 32.0
            points.append((radius * math.cos(theta), radius * math.sin(theta), 0.0))
    return points

tracers.set_tracers(
    pfw,
    annular_tracer_points(),
    variables=["particleID", "speed", "velocityX", "velocityY"],
    write_interval=stopTime / 40.0,
    output_prefix="mmsVortexTracer",
)

# This variable is intentionally unused by PFW.  It keeps the chosen variant
