# -*- coding: utf-8 -*-
"""TemperatureProfile event and temperature-table interpolation verification.

A stationary elastic block is assigned a global temperature table and a
TemperatureProfile event active over the full run.  The expected solution is a
spatially uniform particle temperature equal to the linearly interpolated table
value, with a piecewise-constant temperature rate.  The test isolates the event
path that copies solver-level temperature histories to particle temperature and
temperature-rate fields.
"""

import os

import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}

case_name = os.environ.get("TEMPERATURE_PROFILE_CASE_NAME", "temperatureProfileEvent")
variant_label = os.environ.get("TEMPERATURE_PROFILE_VARIANT_LABEL", "linear table")
stop_time = float(os.environ.get("TEMPERATURE_PROFILE_STOP_TIME", "1.00"))
mid_time = float(os.environ.get("TEMPERATURE_PROFILE_MID_TIME", "0.40"))
t0 = float(os.environ.get("TEMPERATURE_PROFILE_T0", "300.0"))
t_mid = float(os.environ.get("TEMPERATURE_PROFILE_TMID", "260.0"))
t_final = float(os.environ.get("TEMPERATURE_PROFILE_TFINAL", "380.0"))

pfw["caseName"] = case_name
pfw["runDebug"] = False
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = False

refine = int(os.environ.get("TEMPERATURE_PROFILE_REFINE", "1"))
cpp = int(os.environ.get("TEMPERATURE_PROFILE_CPP", "6"))
ppc = int(os.environ.get("TEMPERATURE_PROFILE_PPC", "2"))
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = refine
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * cpp
pfw["nK"] = pfw["zpar"] * cpp
pfw["ppc"] = ppc

pfw["xmin"] = -0.5
pfw["xmax"] = 0.5
pfw["ymin"] = -0.5
pfw["ymax"] = 0.5
pfw["zmin"] = -0.5
pfw["zmax"] = 0.5

pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 4.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.25
pfw["initialDt"] = 1.0e-12
pfw["writeStatistics"] = "all"
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 100.0
pfw["outputType"] = "silo"
pfw["particleFileFields"] = ["Velocity", "MaterialType", "ContactGroup", "SurfaceFlag", "Temperature", "RVector"]
pfw["plottableFields"] = "{ particleTemperature, particleTemperatureRate, particleVelocity }"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity"]

material = material_db.verificationElastic.copy()
pfw["materials"] = [material["name"]]
pfw["materialPropertyString"] = material["materialString"]

block = geom.box(
    "temperature_block",
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
    vel=[0.0, 0.0, 0.0],
    mat=0,
    group=0,
)
pfw["objects"] = [block]

temperature_table = [
    [0.0, t0],
    [mid_time, t_mid],
    [stop_time, t_final],
]
temperature_table_interp_type = "Linear"

def _temperature_table_xml(table):
    return "{ " + ", ".join("{ " + str(row[0]) + ", " + str(row[1]) + " }" for row in table) + " }"

tracers.set_tracers(
    pfw,
    points=[[0.0, 0.0, 0.0]],
    variables=["particleID", "particleTemperature", "particleTemperatureRate", "velocityX", "velocityY", "velocityZ"],
    write_interval=stop_time / 100.0,
    output_prefix="temperatureProfile_center",
)

pfw["useEvents"] = 1
pfw["mpmEventsString"] = f"""
<TemperatureProfile
  startTime="0.0"
  endTime="{stop_time}"
  temperatureTable="{_temperature_table_xml(temperature_table)}"
  interpolationType="{temperature_table_interp_type}"/>
"""

pfw_expected = {
    "variant_label": variant_label,
    "temperature_table": temperature_table,
    "temperature_table_interp_type": temperature_table_interp_type,
    "expected_uniform_temperature": True,
    "stop_time": stop_time,
}
