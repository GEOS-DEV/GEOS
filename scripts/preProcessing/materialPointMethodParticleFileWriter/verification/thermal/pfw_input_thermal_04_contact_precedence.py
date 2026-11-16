"""Contact modes and material/group rule precedence.

Five separated bars compare, from bottom to top: perfect, material-only,
group override, combined material/group override, and insulated contact.
The same two material instances are reused in every bar.

Requires the cumulative GEOS MPM thermal extension, including its patched PFW.
Self-contained: no sibling input, mesh, particle file, or material database.
Units: normalized length/mass/time, with absolute temperatures in kelvin.
Properties are demonstration values, not a physical material calibration.
"""
import math
import pfw_geometryObjects as geom

CELLS_X = 16
H = 1.0 / CELLS_X
HEIGHT = 2.0
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
pfw["tracerOutputPrefix"] = "thermal_04_contact_precedence"
pfw["tracerVariables"] = ["particleID", "mass", "temperature"]
pfw["tracerCoordinates"] = [
    [H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 + H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [1.0 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
]

pfw["materials"] = ["solidA", "solidB"]
pfw["materialPropertyString"] = """
<ThermalElasticIsotropic name="solidA"
    defaultDensity="1" defaultBulkModulus="1" defaultShearModulus="1"
    thermalSpecificHeat="1" thermalConductivity="0.05" thermalReferenceTemperature="300"/>
<ThermalElasticIsotropic name="solidB"
    defaultDensity="1" defaultBulkModulus="1" defaultShearModulus="1"
    thermalSpecificHeat="1" thermalConductivity="0.05" thermalReferenceTemperature="300"/>
"""
pfw["objects"] = []
pfw["tracerCoordinates"] = []
for lane in range(5):
    y0 = lane * 7.0 * H  # Four cells of solid, then three cells of empty gap.
    y1 = y0 + 4.0 * H
    left = geom.box("hot lane %d" % lane, [0.0, y0], [0.5, y1],
                    mat=0, group=2 * lane, dim=2, particleType=2,
                    flaggedSurfaces=[False, False, True, False])
    right = geom.box("cold lane %d" % lane, [0.5, y0], [1.0, y1],
                     mat=1, group=2 * lane + 1, dim=2, particleType=2,
                     flaggedSurfaces=[True, False, False, False])
    pfw["objects"].extend([geom.temperatureWrapper("hot %d" % lane, left, 400.0),
                           geom.temperatureWrapper("cold %d" % lane, right, 300.0)])
    for x in (0.5 - H / 4.0, 0.5 + H / 4.0):
        pfw["tracerCoordinates"].append([x, y0 + 2.0 * H - H / 4.0, 0.0])
pfw["explicitSurfaceNormalInfluence"] = 1.0
pfw["thermal"] = {
    "defaultContactMode": "Resistance",
    "defaultContactResistance": 2.0,
    "areaQuadrature": 40,
    "diagnostics": 1,
    "contactResistance": [
        {"name": "allAB", "materials": ["solidA", "solidB"], "resistance": 0.8},
        {"name": "perfectBottom", "contactGroups": [0, 1], "mode": "Perfect"},
        # Groups 2/3 have no override: they use the material-pair resistance 0.8.
        {"name": "groupOverride", "contactGroups": [4, 5], "resistance": 0.3},
        {"name": "groupFallback", "contactGroups": [6, 7], "resistance": 0.6},
        {"name": "combinedOverride", "contactGroups": [6, 7],
         "materials": ["solidA", "solidB"], "resistance": 0.1},
        {"name": "insulatedTop", "contactGroups": [8, 9], "mode": "Insulated"},
    ],
}
