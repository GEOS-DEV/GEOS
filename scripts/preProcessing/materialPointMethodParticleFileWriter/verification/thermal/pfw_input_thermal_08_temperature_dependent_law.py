"""Constant and temperature-dependent thermal policies.

Two isolated bars use the same elastic base with different thermal policies.
A uniform pulse (0<=t<0.25) tests the caloric law; a later central heater
creates gradients that exercise temperature-dependent conductivity.

Requires the cumulative GEOS MPM thermal extension, including its patched PFW.
Self-contained: no sibling input, mesh, particle file, or material database.
Units: normalized length/mass/time, with absolute temperatures in kelvin.
Properties are demonstration values, not a physical material calibration.
"""
import math
import pfw_geometryObjects as geom

CELLS_X = 16
H = 1.0 / CELLS_X
HEIGHT = 0.6875
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
pfw["tracerOutputPrefix"] = "thermal_08_temperature_dependent_law"
pfw["tracerVariables"] = ["particleID", "mass", "temperature"]
pfw["tracerCoordinates"] = [
    [H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [0.5 + H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
    [1.0 - H / 4.0, HEIGHT / 2.0 - H / 4.0, 0.0],
]

pfw["materials"] = ["constantLaw", "linearLaw"]
pfw["materialPropertyString"] = """
<ThermalElasticIsotropic name="constantLaw"
    defaultDensity="1" defaultBulkModulus="1" defaultShearModulus="1"
    thermalSpecificHeat="1" thermalConductivity="0.05" thermalReferenceTemperature="300"/>
<LinearThermalElasticIsotropic name="linearLaw"
    defaultDensity="1" defaultBulkModulus="1" defaultShearModulus="1"
    thermalSpecificHeat="1" thermalHeatCapacitySlope="0.02"
    thermalConductivity="0.05" thermalConductivitySlope="0.002"
    thermalReferenceTemperature="300"/>
"""
pfw["objects"] = []
pfw["tracerCoordinates"] = []
for lane in range(2):
    y0 = lane * 7.0 * H
    body = geom.box("law comparison %d" % lane, [0.0, y0], [1.0, y0 + 4.0 * H],
                    mat=lane, group=lane, dim=2, particleType=2,
                    flaggedSurfaces=[False, False, False, False])
    body = geom.temperatureWrapper("initial %d" % lane, body, 300.0)
    body.particleColor = lambda pt: 7 if 0.4375 <= pt[0] < 0.5625 else 0
    pfw["objects"].append(body)
    for x in (0.25 - H / 4.0, 0.5 - H / 4.0):
        pfw["tracerCoordinates"].append([x, y0 + 2.0 * H - H / 4.0, 0.0])
pfw["thermal"] = {
    "defaultContactMode": "Insulated",
    "diagnostics": 1,
    "heatSources": [
        {"name": "uniformCaloricPulse", "specificPower": 20.0,
         "startTime": 0.0, "endTime": 0.25},
        {"name": "localizedTransportPulse", "particleColors": [7],
         "specificPower": 100.0, "startTime": 0.25, "endTime": 1.0},
    ],
}
# During the first pulse, u=20*t in both bars. For the linear law:
# dT=(-1+sqrt(1+2*0.02*u))/0.02; c=1+0.02*dT; k=0.05+0.002*dT.
