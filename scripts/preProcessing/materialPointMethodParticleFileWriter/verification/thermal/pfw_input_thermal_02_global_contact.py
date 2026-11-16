"""Global resistance across stationary contact groups.

Two touching halves share one material and use different contact groups.
Change CONTACT_RESISTANCE to 0 for perfect contact or increase it to slow exchange.

Requires the cumulative GEOS MPM thermal extension, including its patched PFW.
Self-contained: no sibling input, mesh, particle file, or material database.
Units: normalized length/mass/time, with absolute temperatures in kelvin.
Properties are demonstration values, not a physical material calibration.
"""
import math
import pfw_geometryObjects as geom

CELLS_X = 16
H = 1.0 / CELLS_X
HEIGHT = 0.25
CELLS_Y = round(HEIGHT / H)
STOP_TIME = 1.0

pfw = {}
pfw["xmin"], pfw["xmax"] = 0.0, 1.0
pfw["ymin"], pfw["ymax"] = 0.0, HEIGHT
pfw["zmin"], pfw["zmax"] = -0.5 * H, 0.5 * H
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]
pfw["xpar"], pfw["ypar"], pfw["zpar"] = 1, 1, 1
# PFW counts one buffer cell on each nonperiodic side.
pfw["nI"], pfw["nJ"], pfw["nK"] = CELLS_X + 2, CELLS_Y + 2, 3
pfw["ppc"] = 2
pfw["mBatch"] = True
pfw["mSubmitJobs"] = True
pfw["autoRestart"] = False
pfw["runDebug"] = True
pfw["mWallTime"] = "00:20:00"
pfw["endTime"] = STOP_TIME
pfw["plotInterval"] = STOP_TIME / 40.0
pfw["restartInterval"] = 10.0 * STOP_TIME
pfw["outputType"] = "silo"  # Change to "silo" if preferred.
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["updateMethod"] = "FLIP"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 0.001
pfw["solverProfiling"] = 0
pfw["subdivideParticles"] = 0
pfw["damageFieldPartitioning"] = 0
pfw["needsNeighborList"] = 0
# Stationary heat transfer is independent of mechanical contact impulses.
pfw["enableContact"] = 0
pfw["frictionCoefficient"] = 0.0
pfw["boundaryConditionTypes"] = [1, 1, 1, 1, 1, 1]  # Symmetry walls.
pfw["particleFileFields"] = [
    "Velocity", "MaterialType", "ContactGroup", "ParticleColor", "SurfaceFlag",
    "Damage", "Temperature", "RVector", "SurfaceNormal", "SurfacePosition",
]
pfw["tracerHistory"] = 1
pfw["tracerCycleInterval"] = 1
pfw["tracerOutputPrefix"] = "thermal_02_global_contact"
pfw["tracerVariables"] = ["particleID", "mass", "temperature"]
pfw["tracerCoordinates"] = [
    [H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 + H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [1.0 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
]

pfw["materials"] = ["solid"]
pfw["materialPropertyString"] = """
<ThermalElasticIsotropic name="solid"
    defaultDensity="1" defaultBulkModulus="1" defaultShearModulus="1"
    thermalSpecificHeat="1" thermalConductivity="0.05"
    thermalReferenceTemperature="300" solverThermalExpansion="0"/>
"""
# Only the internal x faces are flagged. Their normals oppose one another;
# SurfacePosition is a particle-to-interface offset, not an absolute point.
left = geom.box("hot half", [0.0, 0.0], [0.5, HEIGHT], mat=0, group=0,
                dim=2, particleType=2, flaggedSurfaces=[False, False, True, False])
right = geom.box("cold half", [0.5, 0.0], [1.0, HEIGHT], mat=0, group=1,
                 dim=2, particleType=2, flaggedSurfaces=[True, False, False, False])
left = geom.temperatureWrapper("hot", left, 400.0)
right = geom.temperatureWrapper("cold", right, 300.0)
pfw["objects"] = [left, right]
pfw["explicitSurfaceNormalInfluence"] = 1.0

CONTACT_RESISTANCE = 0.5  # Area-specific resistance: length^2 * K / power.
pfw["thermal"] = {
    "enabled": 1,
    "defaultContactMode": "Resistance",
    "defaultContactResistance": CONTACT_RESISTANCE,
    "areaQuadrature": 40,
    "diagnostics": 1,
}
