"""Particle-local heaters and coolers without diffusion.

A per-particle power, a negative specific power, and a color-selected
volumetric source have disjoint supports. Conductivity is zero.
At t=1: heater=310 K, cooler=292.5 K, color 7=305 K, others=300 K.

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
pfw["tracerOutputPrefix"] = "thermal_06_particle_sources"
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
    thermalSpecificHeat="1" thermalConductivity="0"
    thermalReferenceTemperature="300" solverThermalExpansion="0"/>
"""

body = geom.box("source specimen", [0.0, 0.0], [1.0, HEIGHT], dim=2,
                particleType=2, flaggedSurfaces=[False, False, False, False])
body = geom.temperatureWrapper("uniform initial temperature", body, 300.0)
# Attach the color to the outermost wrapper so PFW writes it to each particle.
body.particleColor = lambda pt: 7 if 0.75 <= pt[0] < 0.875 else 0
pfw["objects"] = [body]

# PFW numbers this full rectangle by x column, then y row, starting at 1.
# This ID formula requires serial PFW generation, ppc=2, and no missing particles.
# Both selected particles lie near the center, outside the color-7 stripe.
ROWS = 2 * CELLS_Y
HEATER_ID = (CELLS_X - 1) * ROWS + CELLS_Y
COOLER_ID = (CELLS_X + 1) * ROWS + CELLS_Y
PARTICLE_MASS = (H / 2.0) ** 2 * H  # rho=1; represented plane-strain thickness=H.
pfw["tracerCoordinates"] = [
    [0.5 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 + 3.0 * H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.8, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.25, HEIGHT / 2.0 - H / 4.0, 0.0],
]
pfw["thermal"] = {
    "defaultContactMode": "Insulated",
    "diagnostics": 1,
    "heatSources": [
        {"name": "singleParticleHeater", "particleIDs": [HEATER_ID],
         "power": 20.0 * PARTICLE_MASS, "startTime": 0.0, "endTime": 0.5},
        {"name": "singleParticleCooler", "particleIDs": [COOLER_ID],
         "specificPower": -10.0, "startTime": 0.25, "endTime": 1.0},
        # Selector categories are intersected; entries inside one list are ORed.
        {"name": "volumetricStripe", "targetRegions": ["ParticleRegion1"],
         "contactGroups": [0], "particleColors": [7], "volumetricPower": 5.0},
    ],
}
