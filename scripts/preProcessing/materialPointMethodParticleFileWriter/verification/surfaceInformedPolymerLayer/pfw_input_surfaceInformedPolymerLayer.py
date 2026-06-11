# -*- coding: utf-8 -*-
"""Surface-informed polymer layer verification input.

This single input is dispatched by runTest in four variants:
continuum tension, continuum compression, cohesive-zone tension, and cohesive-zone compression.
The continuum variant places a continuum polymer layer between elastic loading bars and bonds the
bar/polymer interfaces with stiff elastic cohesive zones.  The cohesive-zone variant retains two
elastic blocks separated by a finite-thickness polymer cohesive interface with the same film
thickness.
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
pfw["mWallTime"] = "00:10:00"
pfw["mSubmitJobs"] = True
pfw["autoRestart"] = False

# Geometry and resolution.  Plane strain is used to keep the verification inexpensive.
#
#
# Both branches use the same polymer film thickness.  In the continuum branch this is the
# thickness of the continuum polymer layer; in the CZ branch it is the finite-thickness
# polymer cohesive law thickness.
film_thickness = 0.10
half_film = 0.5 * film_thickness

domain_width = film_thickness
domain_height = 10.*film_thickness

pfw["xmin"] = -0.5 * domain_width
pfw["xmax"] =  0.5 * domain_width
pfw["ymin"] = -0.5 * domain_height
pfw["ymax"] =  0.5 * domain_height
pfw["planeStrain"] = 1
pfw["periodic"] = [False, False, False]

pfw["xpar"] = 1
pfw["ypar"] = 1
pfw["zpar"] = 1

pfw["nI"] = 10
pfw["nJ"] = 100
pfw["nK"] = 3
pfw["ppcx"] = 2
pfw["ppcy"] = 3
pfw["ppcz"] = 2

cell_size = domain_height / max(pfw["nJ"], 1)
pfw["zmin"] = -0.5 * cell_size
pfw["zmax"] =  0.5 * cell_size

# The first version of this verification loaded the layer over one microsecond,
# which is shorter than a polymer-layer acoustic transit time for the validation
# material and therefore produced a dynamic reaction curve rather than a
# constitutive thin-film curve.  Use a deliberately slow cosine ramp so the
# boundary reactions approach the quasi-static analytical update.
stop_time = 20.0
pfw["endTime"] = stop_time
pfw["plotInterval"] = stop_time / 4.0
pfw["restartInterval"] = 2.0 * stop_time
pfw["timeIntegrationOption"] = "ExplicitDynamic"
pfw["cflFactor"] = 0.20
pfw["initialDt"] = 1.0e-16

pfw["damageFieldPartitioning"] = 1

pfw["reactionHistory"] = 1
pfw["reactionWriteInterval"] = stop_time / 250.0
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stop_time / 250.0

# Average the central layer using a fixed window sized to the estimated final compressed polymer
# thickness: h = h0 * (1 - 0.8) = 0.02.  Tension may fracture and empty this window; GEOS guards
# against empty box averages separately so diagnostics do not terminate the run.

compressed_film_half_thickness = 0.5 * film_thickness * (1.0 - 0.80)
pfw["boxAverageMin"] = [pfw["xmin"], -compressed_film_half_thickness, pfw["zmin"]]
pfw["boxAverageMax"] = [pfw["xmax"],  compressed_film_half_thickness, pfw["zmax"]]
pfw["boxAverageResizeWithDomain"] = 0

pfw["writeStatistics"] = "all"
pfw["frictionCoefficient"] = 0.0
pfw["updateMethod"] = "FMPM"
pfw["updateOrder"] = 4
pfw["outputType"] = "silo"

pfw["plotGridFields"] = 1
pfw["gridFieldNames"] = ["gridMass", "gridVelocity", "gridInternalForce"]

# GEOS/Silo output fields. These control what VisIt can plot.
# Keep these separate from particleFileFields, which only controls fields
# written to the initial particle file by PFW.
pfw["plottableFields"] = [
    "particleStress",
    "particleDamage",
    "particlePlasticStrain"
]

# PFW particle-file initialization fields.  These are not VisIt/Silo output controls.
# The continuum patch has no initial surface flags, so it only needs the fields required
# to initialize particle motion, material type, and CPDI particle domains.  The CZ branch
# retains surface fields because the cohesive-zone event builds an interface from flagged
# surfaces, surface normals, and surface positions.
pfw["particleFileFields"] = [
        'Velocity', 'MaterialType', 'ContactGroup', 'SurfaceFlag', 'CZTag', 'RVector',
        "SurfaceNormal",
        "SurfacePosition"
]

# Uniaxial-strain loading.  Lateral faces follow F_xx=1 while the top and bottom faces impose the
# y-stretch.  Both branches use the full bar/block gage length, so the post-processor converts the
# imposed global strain to film strain using the common film thickness.
pfw["prescribedBcTable"] = 0
pfw["boundaryConditionTypes"] = [2, 2, 2, 2, 1, 1]
pfw["prescribedBoundaryFTable"] = 1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"] = [[0.0, 1.0, 1.0, 1.0], [stop_time, 1.0, 1.0 + final_global_strain, 1.0]]

# Materials.  The block is a non-polymer elastic loading material.  The continuum branch
# uses a stiff non-polymer elastic cohesive tie at the bar/polymer interfaces so the polymer
# layer can be tested as a continuum material without mixed-cell bar/polymer bonding.
# The CZ branch uses the polymer cohesive law as the material under test.
block = material_db.surfacePolymerVerificationElasticBlock.copy()
polymer = material_db.vitonFKM75SurfacePolymer.copy()
elastic_tie_cz = material_db.surfacePolymerVerificationElasticTieCohesiveZone.copy()
polymer_cz = material_db.vitonFKM75SurfacePolymerCohesiveZone.copy()

polymer_cz["thickness"] = film_thickness

def continuum_interface_tag(pt):
    # Lower bar/polymer interface gets tag 0 and upper bar/polymer interface gets tag 1.
    # Keeping the two tags separate prevents the reference-CZ builder from pairing surfaces
    # across the polymer layer when both interfaces are present in the same run.
    return 0 if pt[1] < 0.0 else 1

if model_kind == "continuum":
    pfw["useEvents"] = 1
    pfw["materials"] = [polymer["name"], block["name"]]
    pfw["materialPropertyString"] = (
        polymer["materialString"]
        + block["materialString"]
        + elastic_tie_cz["materialString"]
    )

    upperbar = geom.box(
        "upperbar",
        [pfw["xmin"], half_film, pfw["zmin"]],
        [pfw["xmax"], pfw["ymax"], pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=1,
        group=0,
        dim=2,
        # y-min face is bonded to the polymer layer by the stiff elastic CZ.
        flaggedSurfaces=[False, True, False, False],
    )
    layer = geom.box(
        "polymer_layer",
        [pfw["xmin"], -half_film, pfw["zmin"]],
        [pfw["xmax"], half_film, pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=0,
        group=0,
        dim=2,
        # y-min and y-max faces are bonded to the lower and upper bars.
        flaggedSurfaces=[False, True, False, True],
    )
    lowerbar = geom.box(
        "lowerbar",
        [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
        [pfw["xmax"], -half_film, pfw["zmax"]],
        vel=[0.0, 0.0, 0.0],
        mat=1,
        group=0,
        dim=2,
        # y-max face is bonded to the polymer layer by the stiff elastic CZ.
        flaggedSurfaces=[False, False, False, True],
    )
    pfw["objects"] = [
        geom.czTagWrapper(
            "polymer_layer_interface_tag",
            geom.surfaceFlagWrapper("polymer_layer_interface", layer, 3),
            continuum_interface_tag,
        ),
        geom.czTagWrapper(
            "lowerbar_interface_tag",
            geom.surfaceFlagWrapper("lowerbar_interface", lowerbar, 3),
            continuum_interface_tag,
        ),
        geom.czTagWrapper(
            "upperbar_interface_tag",
            geom.surfaceFlagWrapper("upperbar_interface", upperbar, 3),
            continuum_interface_tag,
        ),
    ]
    pfw["cohesiveZoneRegions"] = f"""
<CohesiveZoneRegion
    name="surfacePolymerLowerElasticTieRegion"
    constitutiveModel="{elastic_tie_cz['name']}"
    tag="0"/>
<CohesiveZoneRegion
    name="surfacePolymerUpperElasticTieRegion"
    constitutiveModel="{elastic_tie_cz['name']}"
    tag="1"/>"""
    pfw["mpmEventsString"] = """
<ReferenceCohesiveZones
    name="surfacePolymerElasticTieEvent"
    startTime="0.0"
    regionNames="{surfacePolymerLowerElasticTieRegion, surfacePolymerUpperElasticTieRegion}"
    czVolumeNormalization="1"/>"""
else:
    pfw["useEvents"] = 1
    pfw["materials"] = [block["name"]]
    pfw["materialPropertyString"] = block["materialString"] + polymer_cz["materialString"]
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
    pfw["cohesiveZoneRegions"] = f"""
<CohesiveZoneRegion
    name="surfacePolymerCZRegion"
    constitutiveModel="{polymer_cz['name']}"
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
    "compressive_pressure_strengthening_cap": polymer["compressivePressureStrengtheningCap"],
    "block_bulk_modulus": block["defaultBulkModulus"],
    "block_shear_modulus": block["defaultShearModulus"],
    "elastic_tie_normal_force_constant": elastic_tie_cz["normalForceConstant"] if model_kind == "continuum" else None,
    "elastic_tie_shear_force_constant": elastic_tie_cz["shearForceConstant"] if model_kind == "continuum" else None,
}
