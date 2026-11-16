"""Thermal resistance across a closed DFG crack.

One material, one contact group, and an inserted planar crack at x=0.5.
DFG creates the two mechanical fields; a same-material/group rule controls heat flow.
The fully damaged state is prescribed for this stationary transport demonstration.

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
pfw["tracerOutputPrefix"] = "thermal_05_dfg_crack"
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

body = geom.box("single body", [0.0, 0.0], [1.0, HEIGHT], mat=0, group=0,
                dim=2, particleType=2, flaggedSurfaces=[False, False, False, False])
body = geom.interfaceWrapper("closed crack", body, point=[0.5, 0.0],
                             normal=[1.0, 0.0], dim=2,
                             surfaceFlag=geom.SurfaceFlag.Surface)
body = geom.damageWrapper("prescribed damaged state", body, 1.0)
body = geom.temperatureWrapper("hot and cold crack faces", body,
                               lambda obj, pt: 400.0 if pt[0] < 0.5 else 300.0)
pfw["objects"] = [body]
pfw["damageFieldPartitioning"] = 1
pfw["needsNeighborList"] = 1
pfw["neighborRadius"] = 2.0 * H
pfw["explicitSurfaceNormalInfluence"] = 1.0
pfw["surfaceQualityThreshold"] = 0.0
# Deliberately permit separation in this fully damaged, thin test specimen.
pfw["maxSingleFieldStateFractionForSeparability"] = -1.0
pfw["thinFeatureDFGThreshold"] = 10.0
pfw["thermal"] = {
    "defaultContactMode": "Insulated",
    "areaQuadrature": 40,
    "diagnostics": 1,
    "contactResistance": [
        {"name": "sameGroupCrack", "contactGroups": [0, 0],
         "materials": ["solid", "solid"], "mode": "Resistance", "resistance": 0.5},
    ],
}
