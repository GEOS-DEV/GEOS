# -*- coding: utf-8 -*-
"""Surface-informed polymer layer verification input.

This single input is dispatched by runTest in four variants:
continuum tension, continuum compression, cohesive-zone tension, and cohesive-zone compression.
The geometry is a thin polymer layer between two stiff elastic blocks.  The continuum variant
represents the layer with particles; the cohesive-zone variant replaces the layer by a finite-thickness
cohesive interface with the same film thickness.
"""

import os

import pfw_geometryObjects as geom

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

pfw = {}

case_name = os.environ.get("SURFACE_POLYMER_CASE_NAME", "surfaceInformedPolymerLayer_continuum_tension")
variant_label = os.environ.get("SURFACE_POLYMER_VARIANT_LABEL", "Continuum tension")
model_kind = os.environ.get("SURFACE_POLYMER_MODEL_KIND", "continuum")
loading = os.environ.get("SURFACE_POLYMER_LOADING", "tension")
final_global_strain = float(os.environ.get("SURFACE_POLYMER_FINAL_GLOBAL_STRAIN", "0.10"))

if loading == "compression" and final_global_strain > 0.0:
    final_global_strain *= -1.0
if loading == "tension" and final_global_strain < 0.0:
    final_global_strain *= -1.0

pfw["caseName"] = case_name
pfw["runDebug"] = True
pfw["mBatch"] = True
pfw["mWallTime"] = "00:05:00"
pfw["mSubmitJobs"] = True
pfw["autoRestart"] = False

# Geometry and resolution.  Plane strain is used to keep the verification inexpensive.  The y-span is
# the gage length used by the analytical thin-film reconstruction.
domain_width = 0.50
domain_height = 1.00
film_thickness = 0.10
half_film = 0.5 * film_thickness

pfw["xmin"] = -0.5 * domain_width
pfw["xmax"] =  0.5 * domain_width
pfw["ymin"] = -0.5 * domain_height
pfw["ymax"] =  0.5 * domain_height
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]

refine = int(os.environ.get("SURFACE_POLYMER_REFINE", "1"))
cpp = int(os.environ.get("SURFACE_POLYMER_CPP", "16"))
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1
pfw["nI"] = pfw["xpar"] * cpp
pfw["nJ"] = pfw["ypar"] * 2 * cpp
pfw["nK"] = 3
pfw["ppc"] = 2
cell_size = domain_height / max(pfw["nJ"], 1)
pfw["zmin"] = -0.5 * cell_size
pfw["zmax"] =  0.5 * cell_size

stop_time = 1.0
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 4.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.20
pfw["initialDt"] = 1.0e-16
pfw["cpdiDomainScaling"] = 1
pfw["damageFieldPartitioning"] = 1
pfw["solverProfiling"] = 0
pfw["needsNeighborList"] = 0
pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 250.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 250.0
# Average a narrow material layer for the continuum variant.  The post-processor uses
# boundary reactions as the primary stress metric, so this box is mainly a diagnostic
# and a convenient source of plastic-strain history.
pfw["boxAverageMin"] = [pfw["xmin"], -half_film, pfw["zmin"]]
pfw["boxAverageMax"] = [pfw["xmax"],  half_film, pfw["zmax"]]
pfw["boxAverageResizeWithDomain"] = 1
pfw["writeStatistics"] = "all"
pfw["frictionCoefficient"] = 0.0
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 2
pfw["outputType"] = "silo"
pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity", "gridInternalForce"]
pfw["particleFileFields"] = [
    "Velocity",
    "MaterialType",
    "Stress",
    "Damage",
    "equivalentPlasticStrain",
    "PlasticStrainMagnitude",
    "SurfaceFlag",
    "RVector",
    "SurfaceNormal",
    "SurfacePosition",
]

# Uniaxial-strain loading.  Lateral faces follow F_xx=1 while the top and bottom faces impose the
# y-stretch.  The analytical comparison reconstructs the film strain from the total gage strain.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0, 1.0 + final_global_strain, 1.0]]

# Materials.  The stiff blocks are intentionally much stiffer than the polymer film so the imposed
# gage displacement is concentrated in the layer for the analytical comparison.
block_density = 2.70
block_K = 50.0
block_G = 20.0
block_xml = f"""
<ElasticIsotropic
    name="surfacePolymerElasticBlock"
    defaultDensity="{block_density}"
    defaultBulkModulus="{block_K}"
    defaultShearModulus="{block_G}"/>
"""
polymer = material_db.vitonFKM75SurfacePolymer.copy()
cohesive = material_db.vitonFKM75SurfacePolymerCohesiveZone.copy()
cohesive["thickness"] = film_thickness

if model_kind == "continuum":
    pfw["materials"] = ["surfacePolymerElasticBlock", polymer["name"]]
    pfw["materialPropertyString"] = block_xml + polymer["materialString"]
    bottom = geom.box(
        "bottom_block",
        [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
        [pfw["xmax"], -half_film, pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        dim=2,
    )
    layer = geom.box(
        "polymer_layer",
        [pfw["xmin"], -half_film, pfw["zmin"]],
        [pfw["xmax"],  half_film, pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=1,
        group=1,
        dim=2,
    )
    top = geom.box(
        "top_block",
        [pfw["xmin"], half_film, pfw["zmin"]],
        [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        dim=2,
    )
    pfw["objects"] = [bottom, layer, top]
else:
    pfw["useEvents"] = 1
    pfw["materials"] = ["surfacePolymerElasticBlock"]
    pfw["materialPropertyString"] = block_xml + cohesive["materialString"]
    bottom = geom.box(
        "bottom_block",
        [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
        [pfw["xmax"], 0.0, pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        dim=2,
        flaggedSurfaces=[False, False, False, True],
    )
    top = geom.box(
        "top_block",
        [pfw["xmin"], 0.0, pfw["zmin"]],
        [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        dim=2,
        flaggedSurfaces=[False, True, False, False],
    )
    pfw["objects"] = [geom.surfaceFlagWrapper("bottom_interface", bottom, 3), geom.surfaceFlagWrapper("top_interface", top, 3)]
    pfw["cohesiveZoneRegions"] = """
<CohesiveZoneRegion
    name="surfacePolymerCZRegion"
    constitutiveModel="vitonFKM75SurfacePolymerCohesiveZone"
    tag="0"/>"""
    pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="surfacePolymerCZEvent"
    startTime="0.0"
    regionNames="{surfacePolymerCZRegion}"
    czVolumeNormalization="1"/>"""

pfw_expected = {
    "variant_label": variant_label,
    "case_name": case_name,
    "model_kind": model_kind,
    "loading": loading,
    "final_global_strain": final_global_strain,
    "gage_length": domain_height,
    "film_thickness": film_thickness,
    "stop_time": stop_time,
    "bulk_modulus": polymer["defaultBulkModulus"],
    "shear_modulus": polymer["defaultShearModulus"],
    "default_yield_strength": polymer["defaultYieldStrength"],
    "shear_softening_magnitude": polymer["shearSofteningMagnitude"],
    "shear_softening_shape_parameter1": polymer["shearSofteningShapeParameter1"],
    "shear_softening_shape_parameter2": polymer["shearSofteningShapeParameter2"],
    "strain_hardening_slope": polymer["strainHardeningSlope"],
    "hardening_scale_exponent": polymer["hardeningScaleExponent"],
    "maximum_stretch": polymer["maximumStretch"],
    "pressure_asymmetry_amplitude": polymer["pressureAsymmetryAmplitude"],
}
