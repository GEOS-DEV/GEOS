# MATERIAL PROPERTIES for validation suite.
# units are: mm, us, mg, K, so stress is GPa
#
# Common engineering metals and ceramics, written in the same explicit dictionary style
# as the supplied aluminum and quartz examples.
#
# Notes:
# - Values are representative validation-suite parameters, not certified design allowables.
# - Density is mg/mm^3, numerically equal to g/cm^3.
# - Elastic moduli and strengths are GPa.
# - Fracture toughness is GPa*sqrt(mm).  1 MPa*sqrt(m) = 0.03162277660168379 GPa*sqrt(mm).
# - Fracture energy release rate is GPa*mm.
# - Bulk and shear moduli were computed from representative E and nu where source tables
#   provided Youngs modulus and Poissons ratio rather than K and G directly.
# - Ceramic crackSpeed, maximumStrength, and residual friction slope are validation-model
#   parameters when not directly tabulated material properties.
# - Ceramic Weibull reference volumes in this file are in mm^3.  When the strength source
#   does not publish a specimen geometry, weibullReferenceVolume is set to the full volume
#   of a common ASTM C1161 Type B flexure coupon: 3 mm x 4 mm x 45 mm = 540 mm^3.
#   ASTM C1161 lists 3 x 4 x 45 to 50 mm as a standard-size option for advanced-ceramic
#   flexural testing; the lower end is used here as the nominal reference coupon.
# - For Y-TZP/stabilized zirconia, a reported 3 x 4 x 48 mm four-point-bend bar is used,
#   giving 576 mm^3.
# - For tungsten carbide treated as cemented carbide/hardmetal, an ASTM B406-style
#   transverse-rupture coupon of 0.2 in x 0.25 in x 0.75 in is used, giving
#   5.08 mm x 6.35 mm x 19.05 mm = 614.5149 mm^3.
# - These reference volumes are full coupon/gage volumes, not flexural effective stressed
#   volumes.  Replace them with V_eff values if your size-scaling convention uses Weibull
#   effective volume rather than physical specimen volume.
# - Weibull moduli are representative class values selected from reported size-effect or
#   strength-variability data.  Use batch-specific Weibull statistics when available.
#
# Source groups used while selecting representative properties:
# - MechaniCalc Tables of Material Properties: https://mechanicalc.com/reference/material-properties-tables
# - Engineering ToolBox Ceramic Materials Properties: https://www.engineeringtoolbox.com/ceramics-properties-d_1227.html
# - Superior Technical Ceramics Materials Property Chart: https://ceramics.net/wp-content/uploads/stc-material-property-chart-master-ceramic-full-property-chart-NEW-LOGOS.pdf
# - MakeItFrom material property database: https://www.makeitfrom.com/
# - Vendor/handbook-style datasheets were also used for Ti-6Al-4V, Macor-style glass ceramic,
#   sapphire, fused silica, magnesium AZ31B, Zamak 3, and Co-Cr-Mo F75 entries where needed.
#
# Additional sources used for Weibull reference volumes and moduli:
# - ASTM C1161 summary: https://webstore.ansi.org/standards/astm/astmc1161182023
# - ASTM C1273 summary: https://store.astm.org/c1273-18r24.html
# - ISO 20501 summary: https://www.iso.org/standard/69875.html
# - Morgan Advanced Materials technical-ceramic bending-strength guide:
#   https://www.morgantechnicalceramics.com/en-gb/ceramics-101/mechanical-properties-of-ceramics/bending-strength/
# - Tanaka et al., Weibull parameters for static strengths of ceramics:
#   https://www.jstage.jst.go.jp/article/jsms1963/44/498Appendix/44_498Appendix_51/_pdf/-char/ja
# - 3Y-TZP Weibull characterization summary: https://upcommons.upc.edu/entities/publication/57cb1713-6cde-41c6-b586-b0778e4fff5b
# - Sapphire Weibull flexural-strength summary:
#   https://www.spiedigitallibrary.org/conference-proceedings-of-spie/4452/1/Flexural-strength-of-sapphire--measurements-performed-at-SoRI-UMass/10.1117/12.446895.short
# - Cemented carbide transverse-rupture specimen dimensions, ASTM B406 fixture descriptions:
#   https://www.universalgripco.com/astm-b406

import numpy as np

# -------------------------------------------------------------------------------------------------
# MATERIAL STRING GENERATORS
# -------------------------------------------------------------------------------------------------
# Material entries remain dictionary-compatible, but each entry is finalized below as a
# MaterialProperties object.  MaterialProperties is a dict subclass, so existing code that indexes
# material["defaultDensity"] or material["materialString"] still works.  The materialString value
# is generated from the model-specific generator attached to each finalized entry.

class MaterialProperties(dict):
    """Dictionary-compatible material entry with helper functions exposed as attributes.

    The ``materialString`` field is generated from the model-specific generator every time it is
    requested.  This keeps legacy dictionary use working while avoiding stale XML when a user edits
    a parameter after importing this module.
    """

    def __getitem__(self, key):
        if key == 'materialString':
            return self.refreshMaterialString()
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        if key not in ('materialString',) and getattr(self, '_autoRefreshMaterialString', False):
            self.refreshMaterialString()

    def get(self, key, default=None):
        if key == 'materialString' and self.hasMaterialStringGenerator():
            return self.refreshMaterialString()
        return dict.get(self, key, default)

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __setattr__(self, name, value):
        if name.startswith('_') or name == 'generateMaterialString':
            object.__setattr__(self, name, value)
        else:
            self[name] = value

    def hasMaterialStringGenerator(self):
        return hasattr(self, 'generateMaterialString')

    def refreshMaterialString(self):
        generator = getattr(self, 'generateMaterialString', None)
        if generator is None:
            return dict.__getitem__(self, 'materialString')
        material_string = generator(self)
        dict.__setitem__(self, 'materialString', material_string)
        return material_string

    def copy(self):
        result = MaterialProperties(self)
        for attr, value in getattr(self, '__dict__', {}).items():
            object.__setattr__(result, attr, value)
        return result


def _format_material_value(value):
    if isinstance(value, np.generic):
        value = value.item()
    return str(value)


def _required(material, key):
    if key not in material:
        name = material.get('name', '<unnamed>')
        model = material.get('model', '<unknown model>')
        raise KeyError('Material %s (%s) is missing required parameter %s' % (name, model, key))
    return material[key]


def _format_material_xml(material, attribute_names):
    name = _required(material, 'name')
    model = _required(material, 'model')
    version = material.get('version', 'unknown')
    lines = [
        '<!--%s parameterization of %s model, version: %s-->' % (name, model, version),
        '<%s' % model,
        'name="%s"' % _format_material_value(name),
    ]
    emitted = set(['name'])
    for attr in attribute_names:
        if attr in emitted:
            continue
        if attr in material:
            lines.append('%s="%s"' % (attr, _format_material_value(material[attr])))
            emitted.add(attr)
    lines.append('/>')
    return '\n'.join(lines) + '\n'


def _select_isotropic_elastic_pair(material):
    # GEOS ElasticIsotropic-derived models require one valid pair of elastic constants.  Prefer
    # E/nu when present so a user edit such as material["defaultYoungModulus"] = 2.0 is reflected
    # directly in the generated XML.  Polymer entries below store K/G by default and keep the
    # literature E/nu values under referenceYoungModulus/referencePoissonRatio to avoid over-
    # specifying the GEOS input.
    possible_pairs = [
        ('defaultYoungModulus', 'defaultPoissonRatio'),
        ('defaultBulkModulus', 'defaultShearModulus'),
        ('defaultYoungModulus', 'defaultShearModulus'),
        ('defaultBulkModulus', 'defaultPoissonRatio'),
        ('defaultLambda', 'defaultShearModulus'),
    ]
    for pair in possible_pairs:
        if all(key in material for key in pair):
            return list(pair)
    name = material.get('name', '<unnamed>')
    raise KeyError('Material %s does not define a complete ElasticIsotropic constant pair.' % name)


def generateVonMisesJMaterialString(material):
    return _format_material_xml(material, [
        'defaultDensity',
        'defaultBulkModulus',
        'defaultShearModulus',
        'defaultYieldStrength',
    ])


def generateCeramicDamageMaterialString(material):
    return _format_material_xml(material, [
        'defaultDensity',
        'defaultBulkModulus',
        'defaultShearModulus',
        'tensileStrength',
        'compressiveStrength',
        'maximumStrength',
        'crackSpeed',
        'damagedMaterialFrictionSlope',
        'enableEnergyFailureCriterion',
        'fractureEnergyReleaseRate',
        'fractureToughness',
    ])


def generateElasticIsotropicMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material))


def generateGeomechanicsMaterialString(material):
    return _format_material_xml(material, [
        'defaultDensity',
        'b0', 'b1', 'b2', 'b3', 'b4',
        'dstrendh', 'dfslopedh', 'dpeakI1dh', 'dcrdh',
        'g0', 'g1', 'g2', 'g3', 'g4',
        'p0', 'p1', 'p2', 'p3', 'p4',
        'cr',
        'fluidBulkModulus',
        'fluidInitialPressure',
        't1RateDependence',
        't2RateDependence',
        'peakI1',
        'fSlope',
        'fSlopeFailed',
        'stren',
        'ySlope',
        'beta',
        'fractureEnergyReleaseRate',
        'fractureSofteningExponent',
        'fractureStress',
        'initialTemperature',
        'Q',
        'brittleDuctileTransition',
        'damageEvolutionCriterion',
        'enableBuckling',
        'bucklingLength',
        'bucklingAmplitude',
        'enableCreep',
        'creepC0', 'creepC1', 'creepC2',
        'creepA', 'creepB', 'creepC', 'creepD', 'creepE', 'creepF', 'creepG',
        'strainHardeningN',
        'strainHardeningK',
        'plasticStrainTolerance',
        'stressReturnTolerance',
        'maxAllowedSubcycles',
        'failedStepResponse',
    ])



def generateGraphiteMaterialString(material):
    return _format_material_xml(material, [
        'defaultDensity',
        'defaultDrainedLinearTEC',
        'defaultYoungModulusTransverse',
        'defaultYoungModulusAxial',
        'defaultPoissonRatioTransverse',
        'defaultPoissonRatioAxialTransverse',
        'defaultShearModulusAxialTransverse',
        'defaultYoungModulusTransversePressureDerivative',
        'defaultYoungModulusAxialPressureDerivative',
        'defaultShearModulusAxialTransversePressureDerivative',
        'defaultYoungModulusTransversePressureScale',
        'defaultYoungModulusAxialPressureScale',
        'defaultShearModulusAxialTransversePressureScale',
        'failureStrength',
        'maximumPrincipalStressDamage',
        'crackSpeed',
        'scaleFractureEnergyReleaseRate',
        'fractureEnergyStrengthScaleExponent',
        'basalPlaneFractureEnergyReleaseRate',
        'totalFractureEnergyReleaseRate',
        'damageEvolutionExponent',
        'damagedMaterialFrictionalSlope',
        'distortionShearResponseX2',
        'distortionShearResponseY1',
        'distortionShearResponseY2',
        'distortionShearResponseM1',
        'positiveDistortionShearResponseX2',
        'positiveDistortionShearResponseY1',
        'positiveDistortionShearResponseY2',
        'positiveDistortionShearResponseM1',
        'inPlaneShearResponseX2',
        'inPlaneShearResponseY1',
        'inPlaneShearResponseY2',
        'inPlaneShearResponseM1',
        'coupledShearResponseX2',
        'coupledShearResponseY1',
        'coupledShearResponseY2',
        'coupledShearResponseM1',
        'distortionStrainHardeningC0',
        'inPlaneStrainHardeningC0',
        'coupledStrainHardeningC0',
        'maximumPlasticStrain',
        'alphaL',
        'alphaT',
        'enableCrackTipStressConcentration',
    ])


def generateStrainHardeningPolymerMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material) + [
        'defaultDrainedLinearTEC',
        'defaultYieldStrength',
        'yieldStrengthA',
        'yieldStrengthB',
        'yieldStrengthT0',
        'strainHardeningSlope',
        'strainHardeningSlopeA',
        'strainHardeningSlopeB',
        'strainHardeningSlopeT0',
        'shearSofteningMagnitude',
        'shearSofteningMagnitudeA',
        'shearSofteningMagnitudeB',
        'shearSofteningMagnitudeT0',
        'shearSofteningShapeParameter1',
        'shearSofteningShapeParameter2',
        'maximumStretch',
        'maximumStretchA',
        'maximumStretchB',
        'maximumStretchT0',
        'bulkModulusA',
        'bulkModulusB',
        'bulkModulusT0',
        'shearModulusA',
        'shearModulusB',
        'shearModulusT0',
    ])



def generateSurfaceInformedPolymerMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material) + [
        'defaultDrainedLinearTEC',
        'defaultYieldStrength',
        'shearSofteningMagnitude',
        'shearSofteningShapeParameter1',
        'shearSofteningShapeParameter2',
        'strainHardeningSlope',
        'hardeningScaleExponent',
        'maximumStretch',
        'glassTransitionTemperature',
        'temperatureColdSlope',
        'temperatureHotSlope',
        'temperatureTransitionMagnitude',
        'temperatureTransitionWidth',
        'crystallinity',
        'referenceCrystallinity',
        'crystallinityTransitionWidth',
        'elasticCrystallinityCoeff',
        'yieldStrengthCrystallinityCoeff',
        'pressureAsymmetryAmplitude',
        'pressureAsymmetryWidth',
    ])


def generateSurfaceInformedPolymerCohesiveZoneMaterialString(material):
    return _format_material_xml(material, [
        'thickness',
        'bulkModulus',
        'shearModulus',
        'defaultYieldStrength',
        'shearSofteningMagnitude',
        'shearSofteningShapeParameter1',
        'shearSofteningShapeParameter2',
        'strainHardeningSlope',
        'hardeningScaleExponent',
        'maximumStretch',
        'glassTransitionTemperature',
        'temperatureColdSlope',
        'temperatureHotSlope',
        'temperatureTransitionMagnitude',
        'temperatureTransitionWidth',
        'crystallinity',
        'referenceCrystallinity',
        'crystallinityTransitionWidth',
        'elasticCrystallinityCoeff',
        'yieldStrengthCrystallinityCoeff',
        'pressureAsymmetryAmplitude',
        'pressureAsymmetryWidth',
    ])


def generateHyperelasticMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material) + [
        'defaultDrainedLinearTEC',
    ])


def generateHyperelasticMMSMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material) + [
        'defaultDrainedLinearTEC',
    ])


def generateChiumentiMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material) + [
        'defaultDrainedLinearTEC',
        'criticalLength',
        'failureStrength',
        'energyReleaseRate',
    ])


def generatePerfectlyPlasticMaterialString(material):
    return _format_material_xml(material, ['defaultDensity'] + _select_isotropic_elastic_pair(material) + [
        'defaultDrainedLinearTEC',
        'defaultYieldStress',
    ])


def _select_transverse_isotropic_elastic_set(material):
    engineering_constants = [
        'defaultYoungModulusTransverse',
        'defaultYoungModulusAxial',
        'defaultPoissonRatioTransverse',
        'defaultPoissonRatioAxialTransverse',
        'defaultShearModulusAxialTransverse',
    ]
    stiffness_constants = [
        'defaultC11',
        'defaultC13',
        'defaultC33',
        'defaultC44',
        'defaultC66',
    ]
    if all(key in material for key in engineering_constants):
        return engineering_constants
    if all(key in material for key in stiffness_constants):
        return stiffness_constants
    name = material.get('name', '<unnamed>')
    raise KeyError(
        'Material %s does not define a complete ElasticTransverseIsotropic '
        'engineering-constant or stiffness-constant set.' % name
    )


def generateElasticTransverseIsotropicMaterialString(material):
    return _format_material_xml(material, [
        'defaultDensity',
        'defaultDrainedLinearTEC',
    ] + _select_transverse_isotropic_elastic_set(material))


MATERIAL_STRING_GENERATORS = {
    'VonMisesJ': generateVonMisesJMaterialString,
    'CeramicDamage': generateCeramicDamageMaterialString,
    'ElasticIsotropic': generateElasticIsotropicMaterialString,
    'Geomechanics': generateGeomechanicsMaterialString,
    'Graphite': generateGraphiteMaterialString,
    'StrainHardeningPolymer': generateStrainHardeningPolymerMaterialString,
    'SurfaceInformedPolymer': generateSurfaceInformedPolymerMaterialString,
    'SurfaceInformedPolymerCohesiveZone': generateSurfaceInformedPolymerCohesiveZoneMaterialString,
    'Hyperelastic': generateHyperelasticMaterialString,
    'HyperelasticMMS': generateHyperelasticMMSMaterialString,
    'Chiumenti': generateChiumentiMaterialString,
    'PerfectlyPlastic': generatePerfectlyPlasticMaterialString,
    'ElasticTransverseIsotropic': generateElasticTransverseIsotropicMaterialString,
}


def getMaterialStringGenerator(material_or_model):
    if isinstance(material_or_model, str):
        model = material_or_model
    else:
        model = _required(material_or_model, 'model')
    try:
        return MATERIAL_STRING_GENERATORS[model]
    except KeyError as exc:
        raise KeyError('No material string generator registered for model %s' % model) from exc


def generateMaterialString(material):
    """Generate a GEOS XML material string for a material dictionary."""
    return getMaterialStringGenerator(material)(material)



# -------------------------------------------------------------------------------------------------
# METALS AND ENGINEERING ALLOYS
# Pattern follows the aluminum VonMisesJ example.
# -------------------------------------------------------------------------------------------------

###################################################################################################
# ALUMINUM:
# Generic aluminum: representative density, Young's modulus, Poisson's ratio, and low yield strength.
#
aluminum = {}
aluminum["name"] = "aluminum"
aluminum["version"] = 2605160001
aluminum["model"] = "VonMisesJ"
aluminum["defaultDensity"]=2.7
aluminum["defaultBulkModulus"]=67.647059
aluminum["defaultShearModulus"]=25.93985
aluminum["defaultYieldStrength"]=0.095

aluminum["weibullReferenceVolume"] = 1.0
aluminum["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
aluminum["waveSpeed"] = float(np.sqrt( ( aluminum["defaultBulkModulus"] + 4./3.*aluminum["defaultShearModulus"] ) / aluminum["defaultDensity"] ))

aluminum["materialString"] = generateMaterialString(aluminum)
# #################################################################################################

###################################################################################################
# AL6061 T6:
# Aluminum 6061-T6: representative handbook/table values; yield converted from about 35 ksi.
#
al6061T6 = {}
al6061T6["name"] = "al6061T6"
al6061T6["version"] = 2605160002
al6061T6["model"] = "VonMisesJ"
al6061T6["defaultDensity"]=2.7
al6061T6["defaultBulkModulus"]=67.54902
al6061T6["defaultShearModulus"]=25.902256
al6061T6["defaultYieldStrength"]=0.241

al6061T6["weibullReferenceVolume"] = 1.0
al6061T6["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
al6061T6["waveSpeed"] = float(np.sqrt( ( al6061T6["defaultBulkModulus"] + 4./3.*al6061T6["defaultShearModulus"] ) / al6061T6["defaultDensity"] ))

al6061T6["materialString"] = generateMaterialString(al6061T6)
# #################################################################################################

###################################################################################################
# AL7075 T6:
# Aluminum 7075-T6: representative handbook/table values; yield converted from about 68 ksi.
#
al7075T6 = {}
al7075T6["name"] = "al7075T6"
al7075T6["version"] = 2605160003
al7075T6["model"] = "VonMisesJ"
al7075T6["defaultDensity"]=2.8
al7075T6["defaultBulkModulus"]=69.607843
al7075T6["defaultShearModulus"]=26.691729
al7075T6["defaultYieldStrength"]=0.469

al7075T6["weibullReferenceVolume"] = 1.0
al7075T6["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
al7075T6["waveSpeed"] = float(np.sqrt( ( al7075T6["defaultBulkModulus"] + 4./3.*al7075T6["defaultShearModulus"] ) / al7075T6["defaultDensity"] ))

al7075T6["materialString"] = generateMaterialString(al7075T6)
# #################################################################################################

###################################################################################################
# STEEL:
# Generic mild steel: representative structural-steel density, elastic modulus, and low yield
# strength.
#
steel = {}
steel["name"] = "steel"
steel["version"] = 2605160004
steel["model"] = "VonMisesJ"
steel["defaultDensity"]=7.85
steel["defaultBulkModulus"]=175.0
steel["defaultShearModulus"]=80.769231
steel["defaultYieldStrength"]=0.2

steel["weibullReferenceVolume"] = 1.0
steel["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
steel["waveSpeed"] = float(np.sqrt( ( steel["defaultBulkModulus"] + 4./3.*steel["defaultShearModulus"] ) / steel["defaultDensity"] ))

steel["materialString"] = generateMaterialString(steel)
# #################################################################################################

###################################################################################################
# CARBON STEEL A36:
# ASTM A36-style carbon steel: representative density, elastic modulus, and yield strength.
#
carbonSteelA36 = {}
carbonSteelA36["name"] = "carbonSteelA36"
carbonSteelA36["version"] = 2605160005
carbonSteelA36["model"] = "VonMisesJ"
carbonSteelA36["defaultDensity"]=7.85
carbonSteelA36["defaultBulkModulus"]=166.66667
carbonSteelA36["defaultShearModulus"]=76.923077
carbonSteelA36["defaultYieldStrength"]=0.25

carbonSteelA36["weibullReferenceVolume"] = 1.0
carbonSteelA36["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
carbonSteelA36["waveSpeed"] = float(np.sqrt( ( carbonSteelA36["defaultBulkModulus"] + 4./3.*carbonSteelA36["defaultShearModulus"] ) / carbonSteelA36["defaultDensity"] ))

carbonSteelA36["materialString"] = generateMaterialString(carbonSteelA36)
# #################################################################################################

###################################################################################################
# ALLOY STEEL4140:
# AISI 4140 normalized/annealed-style alloy steel: representative properties from tabulated alloy-
# steel ranges.
#
alloySteel4140 = {}
alloySteel4140["name"] = "alloySteel4140"
alloySteel4140["version"] = 2605160006
alloySteel4140["model"] = "VonMisesJ"
alloySteel4140["defaultDensity"]=7.85
alloySteel4140["defaultBulkModulus"]=162.69841
alloySteel4140["defaultShearModulus"]=79.457364
alloySteel4140["defaultYieldStrength"]=0.415

alloySteel4140["weibullReferenceVolume"] = 1.0
alloySteel4140["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
alloySteel4140["waveSpeed"] = float(np.sqrt( ( alloySteel4140["defaultBulkModulus"] + 4./3.*alloySteel4140["defaultShearModulus"] ) / alloySteel4140["defaultDensity"] ))

alloySteel4140["materialString"] = generateMaterialString(alloySteel4140)
# #################################################################################################

###################################################################################################
# STAINLESS STEEL304:
# AISI 304 stainless steel: representative table values; yield converted from about 30 ksi.
#
stainlessSteel304 = {}
stainlessSteel304["name"] = "stainlessSteel304"
stainlessSteel304["version"] = 2605160007
stainlessSteel304["model"] = "VonMisesJ"
stainlessSteel304["defaultDensity"]=8.0
stainlessSteel304["defaultBulkModulus"]=153.1746
stainlessSteel304["defaultShearModulus"]=74.806202
stainlessSteel304["defaultYieldStrength"]=0.207

stainlessSteel304["weibullReferenceVolume"] = 1.0
stainlessSteel304["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
stainlessSteel304["waveSpeed"] = float(np.sqrt( ( stainlessSteel304["defaultBulkModulus"] + 4./3.*stainlessSteel304["defaultShearModulus"] ) / stainlessSteel304["defaultDensity"] ))

stainlessSteel304["materialString"] = generateMaterialString(stainlessSteel304)
# #################################################################################################

###################################################################################################
# STAINLESS STEEL316:
# AISI 316 stainless steel: representative table values; yield converted from about 30 ksi.
#
stainlessSteel316 = {}
stainlessSteel316["name"] = "stainlessSteel316"
stainlessSteel316["version"] = 2605160008
stainlessSteel316["model"] = "VonMisesJ"
stainlessSteel316["defaultDensity"]=8.0
stainlessSteel316["defaultBulkModulus"]=134.02778
stainlessSteel316["defaultShearModulus"]=76.587302
stainlessSteel316["defaultYieldStrength"]=0.207

stainlessSteel316["weibullReferenceVolume"] = 1.0
stainlessSteel316["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
stainlessSteel316["waveSpeed"] = float(np.sqrt( ( stainlessSteel316["defaultBulkModulus"] + 4./3.*stainlessSteel316["defaultShearModulus"] ) / stainlessSteel316["defaultDensity"] ))

stainlessSteel316["materialString"] = generateMaterialString(stainlessSteel316)
# #################################################################################################

###################################################################################################
# TOOL STEEL A2:
# AISI A2 tool steel: representative density and mechanical properties from tabulated material data.
#
toolSteelA2 = {}
toolSteelA2["name"] = "toolSteelA2"
toolSteelA2["version"] = 2605160009
toolSteelA2["model"] = "VonMisesJ"
toolSteelA2["defaultDensity"]=7.86
toolSteelA2["defaultBulkModulus"]=150.79365
toolSteelA2["defaultShearModulus"]=73.643411
toolSteelA2["defaultYieldStrength"]=0.35

toolSteelA2["weibullReferenceVolume"] = 1.0
toolSteelA2["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
toolSteelA2["waveSpeed"] = float(np.sqrt( ( toolSteelA2["defaultBulkModulus"] + 4./3.*toolSteelA2["defaultShearModulus"] ) / toolSteelA2["defaultDensity"] ))

toolSteelA2["materialString"] = generateMaterialString(toolSteelA2)
# #################################################################################################

###################################################################################################
# MARAGING STEEL300:
# 18Ni maraging steel grade 300: representative value selected from published strength ranges.
#
maragingSteel300 = {}
maragingSteel300["name"] = "maragingSteel300"
maragingSteel300["version"] = 2605160010
maragingSteel300["model"] = "VonMisesJ"
maragingSteel300["defaultDensity"]=8.2
maragingSteel300["defaultBulkModulus"]=158.33333
maragingSteel300["defaultShearModulus"]=73.076923
maragingSteel300["defaultYieldStrength"]=1.4

maragingSteel300["weibullReferenceVolume"] = 1.0
maragingSteel300["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
maragingSteel300["waveSpeed"] = float(np.sqrt( ( maragingSteel300["defaultBulkModulus"] + 4./3.*maragingSteel300["defaultShearModulus"] ) / maragingSteel300["defaultDensity"] ))

maragingSteel300["materialString"] = generateMaterialString(maragingSteel300)
# #################################################################################################

###################################################################################################
# GRAY CAST IRON:
# Gray cast iron: VonMisesJ yield strength is a surrogate for tensile strength because gray iron is
# brittle.
#
grayCastIron = {}
grayCastIron["name"] = "grayCastIron"
grayCastIron["version"] = 2605160011
grayCastIron["model"] = "VonMisesJ"
grayCastIron["defaultDensity"]=7.2
grayCastIron["defaultBulkModulus"]=76.388889
grayCastIron["defaultShearModulus"]=43.650794
grayCastIron["defaultYieldStrength"]=0.17

grayCastIron["weibullReferenceVolume"] = 1.0
grayCastIron["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
grayCastIron["waveSpeed"] = float(np.sqrt( ( grayCastIron["defaultBulkModulus"] + 4./3.*grayCastIron["defaultShearModulus"] ) / grayCastIron["defaultDensity"] ))

grayCastIron["materialString"] = generateMaterialString(grayCastIron)
# #################################################################################################

###################################################################################################
# DUCTILE IRON:
# Ductile/nodular cast iron: representative table values; yield converted from about 40 ksi.
#
ductileIron = {}
ductileIron["name"] = "ductileIron"
ductileIron["version"] = 2605160012
ductileIron["model"] = "VonMisesJ"
ductileIron["defaultDensity"]=7.1
ductileIron["defaultBulkModulus"]=134.12698
ductileIron["defaultShearModulus"]=65.503876
ductileIron["defaultYieldStrength"]=0.276

ductileIron["weibullReferenceVolume"] = 1.0
ductileIron["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
ductileIron["waveSpeed"] = float(np.sqrt( ( ductileIron["defaultBulkModulus"] + 4./3.*ductileIron["defaultShearModulus"] ) / ductileIron["defaultDensity"] ))

ductileIron["materialString"] = generateMaterialString(ductileIron)
# #################################################################################################

###################################################################################################
# TITANIUM GRADE2:
# Commercially pure titanium grade 2: representative density, elastic modulus, and minimum yield.
#
titaniumGrade2 = {}
titaniumGrade2["name"] = "titaniumGrade2"
titaniumGrade2["version"] = 2605160013
titaniumGrade2["model"] = "VonMisesJ"
titaniumGrade2["defaultDensity"]=4.51
titaniumGrade2["defaultBulkModulus"]=109.375
titaniumGrade2["defaultShearModulus"]=39.179104
titaniumGrade2["defaultYieldStrength"]=0.275

titaniumGrade2["weibullReferenceVolume"] = 1.0
titaniumGrade2["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
titaniumGrade2["waveSpeed"] = float(np.sqrt( ( titaniumGrade2["defaultBulkModulus"] + 4./3.*titaniumGrade2["defaultShearModulus"] ) / titaniumGrade2["defaultDensity"] ))

titaniumGrade2["materialString"] = generateMaterialString(titaniumGrade2)
# #################################################################################################

###################################################################################################
# TI64:
# Ti-6Al-4V: representative annealed/AMS-style titanium alloy values; yield converted from about 120
# ksi.
#
ti64 = {}
ti64["name"] = "ti64"
ti64["version"] = 2605160014
ti64["model"] = "VonMisesJ"
ti64["defaultDensity"]=4.43
ti64["defaultBulkModulus"]=118.75
ti64["defaultShearModulus"]=42.537313
ti64["defaultYieldStrength"]=0.828

ti64["weibullReferenceVolume"] = 1.0
ti64["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
ti64["waveSpeed"] = float(np.sqrt( ( ti64["defaultBulkModulus"] + 4./3.*ti64["defaultShearModulus"] ) / ti64["defaultDensity"] ))

ti64["materialString"] = generateMaterialString(ti64)
# #################################################################################################

###################################################################################################
# COPPER:
# Generic copper: representative density, elastic modulus, Poisson's ratio, and yield strength.
#
copper = {}
copper["name"] = "copper"
copper["version"] = 2605160015
copper["model"] = "VonMisesJ"
copper["defaultDensity"]=8.96
copper["defaultBulkModulus"]=121.875
copper["defaultShearModulus"]=43.656716
copper["defaultYieldStrength"]=0.07

copper["weibullReferenceVolume"] = 1.0
copper["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
copper["waveSpeed"] = float(np.sqrt( ( copper["defaultBulkModulus"] + 4./3.*copper["defaultShearModulus"] ) / copper["defaultDensity"] ))

copper["materialString"] = generateMaterialString(copper)
# #################################################################################################

###################################################################################################
# BRASS C260:
# Cartridge brass C26000: yield strength selected from published temper-dependent range.
#
brassC260 = {}
brassC260["name"] = "brassC260"
brassC260["version"] = 2605160016
brassC260["model"] = "VonMisesJ"
brassC260["defaultDensity"]=8.53
brassC260["defaultBulkModulus"]=96.491228
brassC260["defaultShearModulus"]=41.984733
brassC260["defaultYieldStrength"]=0.3

brassC260["weibullReferenceVolume"] = 1.0
brassC260["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
brassC260["waveSpeed"] = float(np.sqrt( ( brassC260["defaultBulkModulus"] + 4./3.*brassC260["defaultShearModulus"] ) / brassC260["defaultDensity"] ))

brassC260["materialString"] = generateMaterialString(brassC260)
# #################################################################################################

###################################################################################################
# PHOSPHOR BRONZE C510:
# Phosphor bronze C510: representative spring-temper value selected from published ranges.
#
phosphorBronzeC510 = {}
phosphorBronzeC510["name"] = "phosphorBronzeC510"
phosphorBronzeC510["version"] = 2605160017
phosphorBronzeC510["model"] = "VonMisesJ"
phosphorBronzeC510["defaultDensity"]=8.86
phosphorBronzeC510["defaultBulkModulus"]=114.58333
phosphorBronzeC510["defaultShearModulus"]=41.044776
phosphorBronzeC510["defaultYieldStrength"]=0.52

phosphorBronzeC510["weibullReferenceVolume"] = 1.0
phosphorBronzeC510["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
phosphorBronzeC510["waveSpeed"] = float(np.sqrt( ( phosphorBronzeC510["defaultBulkModulus"] + 4./3.*phosphorBronzeC510["defaultShearModulus"] ) / phosphorBronzeC510["defaultDensity"] ))

phosphorBronzeC510["materialString"] = generateMaterialString(phosphorBronzeC510)
# #################################################################################################

###################################################################################################
# NICKEL:
# Generic nickel: representative density, elastic modulus, Poisson's ratio, and intermediate yield
# strength.
#
nickel = {}
nickel["name"] = "nickel"
nickel["version"] = 2605160018
nickel["model"] = "VonMisesJ"
nickel["defaultDensity"]=8.9
nickel["defaultBulkModulus"]=187.7193
nickel["defaultShearModulus"]=81.679389
nickel["defaultYieldStrength"]=0.3

nickel["weibullReferenceVolume"] = 1.0
nickel["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
nickel["waveSpeed"] = float(np.sqrt( ( nickel["defaultBulkModulus"] + 4./3.*nickel["defaultShearModulus"] ) / nickel["defaultDensity"] ))

nickel["materialString"] = generateMaterialString(nickel)
# #################################################################################################

###################################################################################################
# INCONEL718:
# Inconel 718 nickel superalloy: representative table values; yield converted from about 120 ksi.
#
inconel718 = {}
inconel718["name"] = "inconel718"
inconel718["version"] = 2605160019
inconel718["model"] = "VonMisesJ"
inconel718["defaultDensity"]=8.22
inconel718["defaultBulkModulus"]=161.11111
inconel718["defaultShearModulus"]=78.682171
inconel718["defaultYieldStrength"]=0.827

inconel718["weibullReferenceVolume"] = 1.0
inconel718["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
inconel718["waveSpeed"] = float(np.sqrt( ( inconel718["defaultBulkModulus"] + 4./3.*inconel718["defaultShearModulus"] ) / inconel718["defaultDensity"] ))

inconel718["materialString"] = generateMaterialString(inconel718)
# #################################################################################################

###################################################################################################
# COBALT CHROME F75:
# Cobalt-chromium-molybdenum alloy: representative property selected from published ranges.
#
cobaltChromeF75 = {}
cobaltChromeF75["name"] = "cobaltChromeF75"
cobaltChromeF75["version"] = 2605160020
cobaltChromeF75["model"] = "VonMisesJ"
cobaltChromeF75["defaultDensity"]=8.4
cobaltChromeF75["defaultBulkModulus"]=182.53968
cobaltChromeF75["defaultShearModulus"]=89.147287
cobaltChromeF75["defaultYieldStrength"]=0.66

cobaltChromeF75["weibullReferenceVolume"] = 1.0
cobaltChromeF75["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
cobaltChromeF75["waveSpeed"] = float(np.sqrt( ( cobaltChromeF75["defaultBulkModulus"] + 4./3.*cobaltChromeF75["defaultShearModulus"] ) / cobaltChromeF75["defaultDensity"] ))

cobaltChromeF75["materialString"] = generateMaterialString(cobaltChromeF75)
# #################################################################################################

###################################################################################################
# MAGNESIUM:
# Generic magnesium: representative density, elastic modulus, Poisson's ratio, and low yield strength.
#
magnesium = {}
magnesium["name"] = "magnesium"
magnesium["version"] = 2605160021
magnesium["model"] = "VonMisesJ"
magnesium["defaultDensity"]=1.74
magnesium["defaultBulkModulus"]=48.888889
magnesium["defaultShearModulus"]=16.296296
magnesium["defaultYieldStrength"]=0.1

magnesium["weibullReferenceVolume"] = 1.0
magnesium["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
magnesium["waveSpeed"] = float(np.sqrt( ( magnesium["defaultBulkModulus"] + 4./3.*magnesium["defaultShearModulus"] ) / magnesium["defaultDensity"] ))

magnesium["materialString"] = generateMaterialString(magnesium)
# #################################################################################################

###################################################################################################
# MAGNESIUM AZ31 B:
# AZ31B magnesium alloy: representative table values for sheet/plate strength.
#
magnesiumAZ31B = {}
magnesiumAZ31B["name"] = "magnesiumAZ31B"
magnesiumAZ31B["version"] = 2605160022
magnesiumAZ31B["model"] = "VonMisesJ"
magnesiumAZ31B["defaultDensity"]=1.77
magnesiumAZ31B["defaultBulkModulus"]=49.777778
magnesiumAZ31B["defaultShearModulus"]=16.592593
magnesiumAZ31B["defaultYieldStrength"]=0.2

magnesiumAZ31B["weibullReferenceVolume"] = 1.0
magnesiumAZ31B["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
magnesiumAZ31B["waveSpeed"] = float(np.sqrt( ( magnesiumAZ31B["defaultBulkModulus"] + 4./3.*magnesiumAZ31B["defaultShearModulus"] ) / magnesiumAZ31B["defaultDensity"] ))

magnesiumAZ31B["materialString"] = generateMaterialString(magnesiumAZ31B)
# #################################################################################################

###################################################################################################
# ZINC:
# Generic zinc: representative density, elastic modulus, Poisson's ratio, and yield strength.
#
zinc = {}
zinc["name"] = "zinc"
zinc["version"] = 2605160023
zinc["model"] = "VonMisesJ"
zinc["defaultDensity"]=7.0
zinc["defaultBulkModulus"]=92.222222
zinc["defaultShearModulus"]=30.740741
zinc["defaultYieldStrength"]=0.12

zinc["weibullReferenceVolume"] = 1.0
zinc["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
zinc["waveSpeed"] = float(np.sqrt( ( zinc["defaultBulkModulus"] + 4./3.*zinc["defaultShearModulus"] ) / zinc["defaultDensity"] ))

zinc["materialString"] = generateMaterialString(zinc)
# #################################################################################################

###################################################################################################
# ZAMAK3:
# Zamak 3 zinc die-casting alloy: representative density, elastic modulus, and yield strength.
#
zamak3 = {}
zamak3["name"] = "zamak3"
zamak3["version"] = 2605160024
zamak3["model"] = "VonMisesJ"
zamak3["defaultDensity"]=6.7
zamak3["defaultBulkModulus"]=64.0
zamak3["defaultShearModulus"]=38.4
zamak3["defaultYieldStrength"]=0.208

zamak3["weibullReferenceVolume"] = 1.0
zamak3["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
zamak3["waveSpeed"] = float(np.sqrt( ( zamak3["defaultBulkModulus"] + 4./3.*zamak3["defaultShearModulus"] ) / zamak3["defaultDensity"] ))

zamak3["materialString"] = generateMaterialString(zamak3)
# #################################################################################################

###################################################################################################
# TUNGSTEN:
# Generic tungsten: representative refractory-metal elastic properties; yield is a conservative
# validation surrogate.
#
tungsten = {}
tungsten["name"] = "tungsten"
tungsten["version"] = 2605160025
tungsten["model"] = "VonMisesJ"
tungsten["defaultDensity"]=19.3
tungsten["defaultBulkModulus"]=261.36364
tungsten["defaultShearModulus"]=134.76562
tungsten["defaultYieldStrength"]=1.0

tungsten["weibullReferenceVolume"] = 1.0
tungsten["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
tungsten["waveSpeed"] = float(np.sqrt( ( tungsten["defaultBulkModulus"] + 4./3.*tungsten["defaultShearModulus"] ) / tungsten["defaultDensity"] ))

tungsten["materialString"] = generateMaterialString(tungsten)
# #################################################################################################

###################################################################################################
# MOLYBDENUM:
# Generic molybdenum: representative refractory-metal properties; yield is an intermediate validation
# value.
#
molybdenum = {}
molybdenum["name"] = "molybdenum"
molybdenum["version"] = 2605160026
molybdenum["model"] = "VonMisesJ"
molybdenum["defaultDensity"]=10.22
molybdenum["defaultBulkModulus"]=255.55556
molybdenum["defaultShearModulus"]=104.54545
molybdenum["defaultYieldStrength"]=0.55

molybdenum["weibullReferenceVolume"] = 1.0
molybdenum["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
molybdenum["waveSpeed"] = float(np.sqrt( ( molybdenum["defaultBulkModulus"] + 4./3.*molybdenum["defaultShearModulus"] ) / molybdenum["defaultDensity"] ))

molybdenum["materialString"] = generateMaterialString(molybdenum)
# #################################################################################################

###################################################################################################
# TANTALUM:
# Generic tantalum: representative density and elastic properties; yield is a validation surrogate.
#
tantalum = {}
tantalum["name"] = "tantalum"
tantalum["version"] = 2605160027
tantalum["model"] = "VonMisesJ"
tantalum["defaultDensity"]=16.6
tantalum["defaultBulkModulus"]=206.66667
tantalum["defaultShearModulus"]=68.888889
tantalum["defaultYieldStrength"]=0.34

tantalum["weibullReferenceVolume"] = 1.0
tantalum["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
tantalum["waveSpeed"] = float(np.sqrt( ( tantalum["defaultBulkModulus"] + 4./3.*tantalum["defaultShearModulus"] ) / tantalum["defaultDensity"] ))

tantalum["materialString"] = generateMaterialString(tantalum)
# #################################################################################################

###################################################################################################
# NIOBIUM:
# Generic niobium: representative refractory-metal density and elastic properties; yield is a
# validation surrogate.
#
niobium = {}
niobium["name"] = "niobium"
niobium["version"] = 2605160028
niobium["model"] = "VonMisesJ"
niobium["defaultDensity"]=8.57
niobium["defaultBulkModulus"]=171.66667
niobium["defaultShearModulus"]=36.785714
niobium["defaultYieldStrength"]=0.2

niobium["weibullReferenceVolume"] = 1.0
niobium["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
niobium["waveSpeed"] = float(np.sqrt( ( niobium["defaultBulkModulus"] + 4./3.*niobium["defaultShearModulus"] ) / niobium["defaultDensity"] ))

niobium["materialString"] = generateMaterialString(niobium)
# #################################################################################################

###################################################################################################
# LEAD:
# Generic lead: representative density, elastic modulus, Poisson's ratio, and low yield strength.
#
lead = {}
lead["name"] = "lead"
lead["version"] = 2605160029
lead["model"] = "VonMisesJ"
lead["defaultDensity"]=11.35
lead["defaultBulkModulus"]=32.857143
lead["defaultShearModulus"]=4.8251748
lead["defaultYieldStrength"]=0.015

lead["weibullReferenceVolume"] = 1.0
lead["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
lead["waveSpeed"] = float(np.sqrt( ( lead["defaultBulkModulus"] + 4./3.*lead["defaultShearModulus"] ) / lead["defaultDensity"] ))

lead["materialString"] = generateMaterialString(lead)
# #################################################################################################

###################################################################################################
# TIN:
# Generic tin: representative density, elastic modulus, Poisson's ratio, and low yield strength.
#
tin = {}
tin["name"] = "tin"
tin["version"] = 2605160030
tin["model"] = "VonMisesJ"
tin["defaultDensity"]=7.31
tin["defaultBulkModulus"]=40.196078
tin["defaultShearModulus"]=15.413534
tin["defaultYieldStrength"]=0.012

tin["weibullReferenceVolume"] = 1.0
tin["weibullModulus"] = 0.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
tin["waveSpeed"] = float(np.sqrt( ( tin["defaultBulkModulus"] + 4./3.*tin["defaultShearModulus"] ) / tin["defaultDensity"] ))

tin["materialString"] = generateMaterialString(tin)
# #################################################################################################


# -------------------------------------------------------------------------------------------------
# ENGINEERING CERAMICS
# Pattern follows the quartz CeramicDamage example.
# -------------------------------------------------------------------------------------------------

###################################################################################################
# ALUMINA995:
# 99.5 pct alumina: representative density, elastic modulus, strengths, and fracture toughness from
# ceramics property tables.
#
alumina995 = {}
alumina995["name"] = "alumina995"
alumina995["version"] = 2605160101
alumina995["model"] = "CeramicDamage"

# Constitutive model parameters:
alumina995["defaultDensity"]=3.88
alumina995["defaultBulkModulus"]=233.95062
alumina995["defaultShearModulus"]=154.06504

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
alumina995["tensileStrength"]=0.172
alumina995["compressiveStrength"]=2.137
alumina995["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
alumina995["crackSpeed"]=3.1924266

# Damaged/residual material behavior:
alumina995["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
alumina995["enableEnergyFailureCriterion"]=1
alumina995["fractureEnergyReleaseRate"]=5.0603625e-05
alumina995["fractureToughness"]=0.14230249

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the full ASTM C1161 Type B flexure coupon volume:
# 3 mm x 4 mm x 45 mm = 540 mm^3; the source property table did not publish a specimen geometry.
# m=10 is a representative dense alumina value within typical technical-ceramic variability.
alumina995["weibullReferenceVolume"] = 540.0
alumina995["weibullModulus"] = 10.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
alumina995["waveSpeed"] = float(np.sqrt( ( alumina995["defaultBulkModulus"] + 4./3.*alumina995["defaultShearModulus"] ) / alumina995["defaultDensity"] ))

alumina995["materialString"] = generateMaterialString(alumina995)
# #################################################################################################

###################################################################################################
# ZIRCONIA YTZP:
# Yttria-stabilized zirconia: representative structural ceramic values; tensile strength approximated
# from flexural/tensile ranges.
#
zirconiaYTZP = {}
zirconiaYTZP["name"] = "zirconiaYTZP"
zirconiaYTZP["version"] = 2605160102
zirconiaYTZP["model"] = "CeramicDamage"

# Constitutive model parameters:
zirconiaYTZP["defaultDensity"]=6.07
zirconiaYTZP["defaultBulkModulus"]=175.0
zirconiaYTZP["defaultShearModulus"]=80.769231

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
zirconiaYTZP["tensileStrength"]=0.69
zirconiaYTZP["compressiveStrength"]=2.485
zirconiaYTZP["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
zirconiaYTZP["crackSpeed"]=2.0473114

# Damaged/residual material behavior:
zirconiaYTZP["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
zirconiaYTZP["enableEnergyFailureCriterion"]=1
zirconiaYTZP["fractureEnergyReleaseRate"]=0.00043333334
zirconiaYTZP["fractureToughness"]=0.31622777

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume uses a reported 3Y-TZP four-point-bend bar:
# 3 mm x 4 mm x 48 mm = 576 mm^3.
# m=10.9 is the reported as-received 3Y-TZP flexural Weibull modulus.
zirconiaYTZP["weibullReferenceVolume"] = 576.0
zirconiaYTZP["weibullModulus"] = 10.9

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
zirconiaYTZP["waveSpeed"] = float(np.sqrt( ( zirconiaYTZP["defaultBulkModulus"] + 4./3.*zirconiaYTZP["defaultShearModulus"] ) / zirconiaYTZP["defaultDensity"] ))

zirconiaYTZP["materialString"] = generateMaterialString(zirconiaYTZP)
# #################################################################################################

###################################################################################################
# SILICON CARBIDE:
# Silicon carbide: representative table values for density, elastic modulus, strengths, and fracture
# toughness.
#
siliconCarbide = {}
siliconCarbide["name"] = "siliconCarbide"
siliconCarbide["version"] = 2605160103
siliconCarbide["model"] = "CeramicDamage"

# Constitutive model parameters:
siliconCarbide["defaultDensity"]=3.15
siliconCarbide["defaultBulkModulus"]=226.26263
siliconCarbide["defaultShearModulus"]=191.45299

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
siliconCarbide["tensileStrength"]=0.241
siliconCarbide["compressiveStrength"]=3.306
siliconCarbide["maximumStrength"]=6.612

# Dynamic/time-to-failure style parameter:
siliconCarbide["crackSpeed"]=3.7091905

# Damaged/residual material behavior:
siliconCarbide["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
siliconCarbide["enableEnergyFailureCriterion"]=1
siliconCarbide["fractureEnergyReleaseRate"]=3.4682145e-05
siliconCarbide["fractureToughness"]=0.12649111

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=10.7 is from pooled sintered silicon-carbide static-strength Weibull data.
siliconCarbide["weibullReferenceVolume"] = 540.0
siliconCarbide["weibullModulus"] = 10.7

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
siliconCarbide["waveSpeed"] = float(np.sqrt( ( siliconCarbide["defaultBulkModulus"] + 4./3.*siliconCarbide["defaultShearModulus"] ) / siliconCarbide["defaultDensity"] ))

siliconCarbide["materialString"] = generateMaterialString(siliconCarbide)
# #################################################################################################

###################################################################################################
# SILICON NITRIDE:
# Silicon nitride: representative table values; compressive strength selected as a typical structural-
# ceramic value.
#
siliconNitride = {}
siliconNitride["name"] = "siliconNitride"
siliconNitride["version"] = 2605160104
siliconNitride["model"] = "CeramicDamage"

# Constitutive model parameters:
siliconNitride["defaultDensity"]=3.25
siliconNitride["defaultBulkModulus"]=227.27273
siliconNitride["defaultShearModulus"]=117.1875

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
siliconNitride["tensileStrength"]=0.537
siliconNitride["compressiveStrength"]=2.5
siliconNitride["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
siliconNitride["crackSpeed"]=3.2589307

# Damaged/residual material behavior:
siliconNitride["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
siliconNitride["enableEnergyFailureCriterion"]=1
siliconNitride["fractureEnergyReleaseRate"]=0.000110592
siliconNitride["fractureToughness"]=0.18973666

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=16.5 is from pooled sintered silicon-nitride static-strength Weibull data.
siliconNitride["weibullReferenceVolume"] = 540.0
siliconNitride["weibullModulus"] = 16.5

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
siliconNitride["waveSpeed"] = float(np.sqrt( ( siliconNitride["defaultBulkModulus"] + 4./3.*siliconNitride["defaultShearModulus"] ) / siliconNitride["defaultDensity"] ))

siliconNitride["materialString"] = generateMaterialString(siliconNitride)
# #################################################################################################

###################################################################################################
# BORON CARBIDE:
# Boron carbide: representative values selected from published density, elastic, strength, and K_IC
# ranges.
#
boronCarbide = {}
boronCarbide["name"] = "boronCarbide"
boronCarbide["version"] = 2605160105
boronCarbide["model"] = "CeramicDamage"

# Constitutive model parameters:
boronCarbide["defaultDensity"]=2.52
boronCarbide["defaultBulkModulus"]=208.33333
boronCarbide["defaultShearModulus"]=169.49153

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
boronCarbide["tensileStrength"]=0.35
boronCarbide["compressiveStrength"]=1.95
boronCarbide["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
boronCarbide["crackSpeed"]=3.9384643

# Damaged/residual material behavior:
boronCarbide["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
boronCarbide["enableEnergyFailureCriterion"]=1
boronCarbide["fractureEnergyReleaseRate"]=2.1771e-05
boronCarbide["fractureToughness"]=0.09486833

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=12 is a representative engineering value for dense boron-carbide strength variability.
boronCarbide["weibullReferenceVolume"] = 540.0
boronCarbide["weibullModulus"] = 12.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
boronCarbide["waveSpeed"] = float(np.sqrt( ( boronCarbide["defaultBulkModulus"] + 4./3.*boronCarbide["defaultShearModulus"] ) / boronCarbide["defaultDensity"] ))

boronCarbide["materialString"] = generateMaterialString(boronCarbide)
# #################################################################################################

###################################################################################################
# ALUMINUM NITRIDE:
# Aluminum nitride: representative table values; compressive strength selected from common ceramic-
# property tables.
#
aluminumNitride = {}
aluminumNitride["name"] = "aluminumNitride"
aluminumNitride["version"] = 2605160106
aluminumNitride["model"] = "CeramicDamage"

# Constitutive model parameters:
aluminumNitride["defaultDensity"]=3.2
aluminumNitride["defaultBulkModulus"]=198.71795
aluminumNitride["defaultShearModulus"]=125.0

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
aluminumNitride["tensileStrength"]=0.139
aluminumNitride["compressiveStrength"]=2.1
aluminumNitride["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
aluminumNitride["crackSpeed"]=3.2056891

# Damaged/residual material behavior:
aluminumNitride["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
aluminumNitride["enableEnergyFailureCriterion"]=1
aluminumNitride["fractureEnergyReleaseRate"]=2.736e-05
aluminumNitride["fractureToughness"]=0.09486833

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=10 is a conservative representative value for AlN; optimized AlN can be higher.
aluminumNitride["weibullReferenceVolume"] = 540.0
aluminumNitride["weibullModulus"] = 10.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
aluminumNitride["waveSpeed"] = float(np.sqrt( ( aluminumNitride["defaultBulkModulus"] + 4./3.*aluminumNitride["defaultShearModulus"] ) / aluminumNitride["defaultDensity"] ))

aluminumNitride["materialString"] = generateMaterialString(aluminumNitride)
# #################################################################################################

###################################################################################################
# BORON NITRIDE:
# Hot-pressed/hexagonal boron nitride: representative values selected from broad published ranges;
# K_IC is an estimate.
#
boronNitride = {}
boronNitride["name"] = "boronNitride"
boronNitride["version"] = 2605160107
boronNitride["model"] = "CeramicDamage"

# Constitutive model parameters:
boronNitride["defaultDensity"]=2.0
boronNitride["defaultBulkModulus"]=17.094017
boronNitride["defaultShearModulus"]=18.018018

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
boronNitride["tensileStrength"]=0.06
boronNitride["compressiveStrength"]=0.15
boronNitride["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
boronNitride["crackSpeed"]=1.3602617

# Damaged/residual material behavior:
boronNitride["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
boronNitride["enableEnergyFailureCriterion"]=1
boronNitride["fractureEnergyReleaseRate"]=5.5569375e-05
boronNitride["fractureToughness"]=0.047434165

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=5 reflects higher expected scatter for hot-pressed/hexagonal BN because of anisotropy and porosity.
boronNitride["weibullReferenceVolume"] = 540.0
boronNitride["weibullModulus"] = 5.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
boronNitride["waveSpeed"] = float(np.sqrt( ( boronNitride["defaultBulkModulus"] + 4./3.*boronNitride["defaultShearModulus"] ) / boronNitride["defaultDensity"] ))

boronNitride["materialString"] = generateMaterialString(boronNitride)
# #################################################################################################

###################################################################################################
# TUNGSTEN CARBIDE:
# Tungsten carbide/cemented carbide style entry: representative density, elastic modulus, strengths,
# and K_IC.
#
tungstenCarbide = {}
tungstenCarbide["name"] = "tungstenCarbide"
tungstenCarbide["version"] = 2605160108
tungstenCarbide["model"] = "CeramicDamage"

# Constitutive model parameters:
tungstenCarbide["defaultDensity"]=15.6
tungstenCarbide["defaultBulkModulus"]=403.84615
tungstenCarbide["defaultShearModulus"]=254.03226

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
tungstenCarbide["tensileStrength"]=0.35
tungstenCarbide["compressiveStrength"]=4.78
tungstenCarbide["maximumStrength"]=9.56

# Dynamic/time-to-failure style parameter:
tungstenCarbide["crackSpeed"]=2.0697768

# Damaged/residual material behavior:
tungstenCarbide["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
tungstenCarbide["enableEnergyFailureCriterion"]=1
tungstenCarbide["fractureEnergyReleaseRate"]=0.00021540572
tungstenCarbide["fractureToughness"]=0.37947332

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume uses a cemented-carbide transverse-rupture coupon:
# 0.2 in x 0.25 in x 0.75 in = 5.08 mm x 6.35 mm x 19.05 mm = 614.5149 mm^3.
# m=12 is a representative value for a dense carbide/cemented-carbide engineering entry.
tungstenCarbide["weibullReferenceVolume"] = 614.5149
tungstenCarbide["weibullModulus"] = 12.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
tungstenCarbide["waveSpeed"] = float(np.sqrt( ( tungstenCarbide["defaultBulkModulus"] + 4./3.*tungstenCarbide["defaultShearModulus"] ) / tungstenCarbide["defaultDensity"] ))

tungstenCarbide["materialString"] = generateMaterialString(tungstenCarbide)
# #################################################################################################

###################################################################################################
# TITANIUM CARBIDE:
# Titanium carbide: density, elastic modulus, shear modulus, Poisson's ratio, and tensile strength
# from published table; compression and K_IC are validation estimates.
#
titaniumCarbide = {}
titaniumCarbide["name"] = "titaniumCarbide"
titaniumCarbide["version"] = 2605160109
titaniumCarbide["model"] = "CeramicDamage"

# Constitutive model parameters:
titaniumCarbide["defaultDensity"]=4.93
titaniumCarbide["defaultBulkModulus"]=238.09524
titaniumCarbide["defaultShearModulus"]=189.87342

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
titaniumCarbide["tensileStrength"]=0.258
titaniumCarbide["compressiveStrength"]=2.5
titaniumCarbide["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
titaniumCarbide["crackSpeed"]=2.9947006

# Damaged/residual material behavior:
titaniumCarbide["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
titaniumCarbide["enableEnergyFailureCriterion"]=1
titaniumCarbide["fractureEnergyReleaseRate"]=3.4338669e-05
titaniumCarbide["fractureToughness"]=0.12649111

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=11 is a representative dense carbide value selected from typical technical-ceramic variability.
titaniumCarbide["weibullReferenceVolume"] = 540.0
titaniumCarbide["weibullModulus"] = 11.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
titaniumCarbide["waveSpeed"] = float(np.sqrt( ( titaniumCarbide["defaultBulkModulus"] + 4./3.*titaniumCarbide["defaultShearModulus"] ) / titaniumCarbide["defaultDensity"] ))

titaniumCarbide["materialString"] = generateMaterialString(titaniumCarbide)
# #################################################################################################

###################################################################################################
# TITANIUM NITRIDE:
# Titanium nitride: density and elastic constants from published table; K_IC selected from reported
# coating/film ranges; strengths are validation estimates.
#
titaniumNitride = {}
titaniumNitride["name"] = "titaniumNitride"
titaniumNitride["version"] = 2605160110
titaniumNitride["model"] = "CeramicDamage"

# Constitutive model parameters:
titaniumNitride["defaultDensity"]=5.21
titaniumNitride["defaultBulkModulus"]=366.66667
titaniumNitride["defaultShearModulus"]=220.0

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
titaniumNitride["tensileStrength"]=0.4
titaniumNitride["compressiveStrength"]=4.0
titaniumNitride["maximumStrength"]=8.0

# Dynamic/time-to-failure style parameter:
titaniumNitride["crackSpeed"]=3.3765591

# Damaged/residual material behavior:
titaniumNitride["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
titaniumNitride["enableEnergyFailureCriterion"]=1
titaniumNitride["fractureEnergyReleaseRate"]=1.5340909e-05
titaniumNitride["fractureToughness"]=0.09486833

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=11 is a representative dense nitride/cermet-like value selected from typical technical-ceramic variability.
titaniumNitride["weibullReferenceVolume"] = 540.0
titaniumNitride["weibullModulus"] = 11.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
titaniumNitride["waveSpeed"] = float(np.sqrt( ( titaniumNitride["defaultBulkModulus"] + 4./3.*titaniumNitride["defaultShearModulus"] ) / titaniumNitride["defaultDensity"] ))

titaniumNitride["materialString"] = generateMaterialString(titaniumNitride)
# #################################################################################################

###################################################################################################
# FUSED SILICA:
# Fused silica/quartz glass: density, elastic constants, tensile and compressive strengths from
# manufacturer data; K_IC from published fused-silica tests.
#
fusedSilica = {}
fusedSilica["name"] = "fusedSilica"
fusedSilica["version"] = 2605160111
fusedSilica["model"] = "CeramicDamage"

# Constitutive model parameters:
fusedSilica["defaultDensity"]=2.2
fusedSilica["defaultBulkModulus"]=36.868687
fusedSilica["defaultShearModulus"]=31.196581

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
fusedSilica["tensileStrength"]=0.054
fusedSilica["compressiveStrength"]=1.14
fusedSilica["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
fusedSilica["crackSpeed"]=1.7916183

# Damaged/residual material behavior:
fusedSilica["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
fusedSilica["enableEnergyFailureCriterion"]=1
fusedSilica["fractureEnergyReleaseRate"]=7.4827908e-06
fusedSilica["fractureToughness"]=0.023717082

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=5 reflects surface-finish dominated strength scatter typical of glass/silica strength measurements.
fusedSilica["weibullReferenceVolume"] = 540.0
fusedSilica["weibullModulus"] = 5.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
fusedSilica["waveSpeed"] = float(np.sqrt( ( fusedSilica["defaultBulkModulus"] + 4./3.*fusedSilica["defaultShearModulus"] ) / fusedSilica["defaultDensity"] ))

fusedSilica["materialString"] = generateMaterialString(fusedSilica)
# #################################################################################################

###################################################################################################
# QUARTZ:
# Single-crystal alpha-quartz entry copied from the provided validation-suite example.
#
quartz = {}
quartz["name"] = "quartz"
quartz["version"] = 2605160112
quartz["model"] = "CeramicDamage"

# Constitutive model parameters:
quartz["defaultDensity"]=2.648
quartz["defaultBulkModulus"]=37.8
quartz["defaultShearModulus"]=44.4

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
quartz["tensileStrength"]=0.183
quartz["compressiveStrength"]=1.46
quartz["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
quartz["crackSpeed"]=1.8

# Damaged/residual material behavior:
quartz["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
quartz["enableEnergyFailureCriterion"]=1
quartz["fractureEnergyReleaseRate"]=9.5e-06
quartz["fractureToughness"]=0.03

# Used by pfw to compute strength-scale for the material:
# The original quartz validation entry did not identify a tensile-test specimen geometry.
# Use the ASTM C1161 Type B full coupon volume fallback: 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=5 is retained as a conservative flaw-scatter value from the original example.
quartz["weibullReferenceVolume"] = 540.0
quartz["weibullModulus"] = 5.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
quartz["waveSpeed"] = float(np.sqrt( ( quartz["defaultBulkModulus"] + 4./3.*quartz["defaultShearModulus"] ) / quartz["defaultDensity"] ))

quartz["materialString"] = generateMaterialString(quartz)
# #################################################################################################

###################################################################################################
# SAPPHIRE:
# Sapphire/single-crystal alumina: representative orientation-averaged mechanical data from
# manufacturer table.
#
sapphire = {}
sapphire["name"] = "sapphire"
sapphire["version"] = 2605160113
sapphire["model"] = "CeramicDamage"

# Constitutive model parameters:
sapphire["defaultDensity"]=3.98
sapphire["defaultBulkModulus"]=317.46032
sapphire["defaultShearModulus"]=155.03876

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
sapphire["tensileStrength"]=0.35
sapphire["compressiveStrength"]=2.0
sapphire["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
sapphire["crackSpeed"]=3.4428602

# Damaged/residual material behavior:
sapphire["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
sapphire["enableEnergyFailureCriterion"]=1
sapphire["fractureEnergyReleaseRate"]=9.1589999e-06
sapphire["fractureToughness"]=0.063245553

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=4 is based on reported sapphire flexural Weibull moduli near 3.4 to 4.1, reflecting surface/orientation scatter.
sapphire["weibullReferenceVolume"] = 540.0
sapphire["weibullModulus"] = 4.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
sapphire["waveSpeed"] = float(np.sqrt( ( sapphire["defaultBulkModulus"] + 4./3.*sapphire["defaultShearModulus"] ) / sapphire["defaultDensity"] ))

sapphire["materialString"] = generateMaterialString(sapphire)
# #################################################################################################

###################################################################################################
# MULLITE:
# Mullite: representative values from engineering-ceramic property tables; K_IC is a typical
# validation value.
#
mullite = {}
mullite["name"] = "mullite"
mullite["version"] = 2605160114
mullite["model"] = "CeramicDamage"

# Constitutive model parameters:
mullite["defaultDensity"]=3.0
mullite["defaultBulkModulus"]=114.74359
mullite["defaultShearModulus"]=72.177419

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
mullite["tensileStrength"]=0.138
mullite["compressiveStrength"]=1.034
mullite["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
mullite["crackSpeed"]=2.5158308

# Damaged/residual material behavior:
mullite["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
mullite["enableEnergyFailureCriterion"]=1
mullite["fractureEnergyReleaseRate"]=4.738324e-05
mullite["fractureToughness"]=0.09486833

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=8 reflects moderate scatter typical of refractory/traditional structural ceramics.
mullite["weibullReferenceVolume"] = 540.0
mullite["weibullModulus"] = 8.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
mullite["waveSpeed"] = float(np.sqrt( ( mullite["defaultBulkModulus"] + 4./3.*mullite["defaultShearModulus"] ) / mullite["defaultDensity"] ))

mullite["materialString"] = generateMaterialString(mullite)
# #################################################################################################

###################################################################################################
# CORDIERITE:
# Cordierite: representative values from engineering-ceramic property tables; tensile strength and
# K_IC are validation estimates.
#
cordierite = {}
cordierite["name"] = "cordierite"
cordierite["version"] = 2605160115
cordierite["model"] = "CeramicDamage"

# Constitutive model parameters:
cordierite["defaultDensity"]=2.0
cordierite["defaultBulkModulus"]=66.025641
cordierite["defaultShearModulus"]=41.532258

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
cordierite["tensileStrength"]=0.019
cordierite["compressiveStrength"]=0.165
cordierite["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
cordierite["crackSpeed"]=2.3373253

# Damaged/residual material behavior:
cordierite["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
cordierite["enableEnergyFailureCriterion"]=1
cordierite["fractureEnergyReleaseRate"]=2.0586408e-05
cordierite["fractureToughness"]=0.047434165

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=6 reflects higher scatter expected for low-strength, thermal-shock ceramic bodies.
cordierite["weibullReferenceVolume"] = 540.0
cordierite["weibullModulus"] = 6.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
cordierite["waveSpeed"] = float(np.sqrt( ( cordierite["defaultBulkModulus"] + 4./3.*cordierite["defaultShearModulus"] ) / cordierite["defaultDensity"] ))

cordierite["materialString"] = generateMaterialString(cordierite)
# #################################################################################################

###################################################################################################
# STEATITE:
# Steatite: representative values from engineering-ceramic property tables; K_IC is a validation
# estimate.
#
steatite = {}
steatite["name"] = "steatite"
steatite["version"] = 2605160116
steatite["model"] = "CeramicDamage"

# Constitutive model parameters:
steatite["defaultDensity"]=2.75
steatite["defaultBulkModulus"]=66.025641
steatite["defaultShearModulus"]=41.532258

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
steatite["tensileStrength"]=0.103
steatite["compressiveStrength"]=0.586
steatite["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
steatite["crackSpeed"]=1.9932777

# Damaged/residual material behavior:
steatite["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
steatite["enableEnergyFailureCriterion"]=1
steatite["fractureEnergyReleaseRate"]=2.0586408e-05
steatite["fractureToughness"]=0.047434165

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=8 reflects moderate scatter for electrical steatite ceramics.
steatite["weibullReferenceVolume"] = 540.0
steatite["weibullModulus"] = 8.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
steatite["waveSpeed"] = float(np.sqrt( ( steatite["defaultBulkModulus"] + 4./3.*steatite["defaultShearModulus"] ) / steatite["defaultDensity"] ))

steatite["materialString"] = generateMaterialString(steatite)
# #################################################################################################

###################################################################################################
# PORCELAIN:
# Technical porcelain: representative engineering-ceramic values; K_IC and tensile strength are
# validation estimates.
#
porcelain = {}
porcelain["name"] = "porcelain"
porcelain["version"] = 2605160117
porcelain["model"] = "CeramicDamage"

# Constitutive model parameters:
porcelain["defaultDensity"]=2.4
porcelain["defaultBulkModulus"]=41.666667
porcelain["defaultShearModulus"]=28.688525

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
porcelain["tensileStrength"]=0.025
porcelain["compressiveStrength"]=0.25
porcelain["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
porcelain["crackSpeed"]=1.7311633

# Damaged/residual material behavior:
porcelain["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
porcelain["enableEnergyFailureCriterion"]=1
porcelain["fractureEnergyReleaseRate"]=1.3594286e-05
porcelain["fractureToughness"]=0.031622777

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=10 is selected from reported porcelain/dental-porcelain flexural Weibull values around 9 to 14.
porcelain["weibullReferenceVolume"] = 540.0
porcelain["weibullModulus"] = 10.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
porcelain["waveSpeed"] = float(np.sqrt( ( porcelain["defaultBulkModulus"] + 4./3.*porcelain["defaultShearModulus"] ) / porcelain["defaultDensity"] ))

porcelain["materialString"] = generateMaterialString(porcelain)
# #################################################################################################

###################################################################################################
# MACOR GLASS CERAMIC:
# Machinable glass ceramic/Macor-style entry: density, flexural/compressive strength, modulus, and
# K_IC from vendor data; tensile is estimated below flexural strength.
#
macorGlassCeramic = {}
macorGlassCeramic["name"] = "macorGlassCeramic"
macorGlassCeramic["version"] = 2605160118
macorGlassCeramic["model"] = "CeramicDamage"

# Constitutive model parameters:
macorGlassCeramic["defaultDensity"]=2.52
macorGlassCeramic["defaultBulkModulus"]=52.380952
macorGlassCeramic["defaultShearModulus"]=25.581395

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
macorGlassCeramic["tensileStrength"]=0.05
macorGlassCeramic["compressiveStrength"]=0.345
macorGlassCeramic["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
macorGlassCeramic["crackSpeed"]=1.7575295

# Damaged/residual material behavior:
macorGlassCeramic["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
macorGlassCeramic["enableEnergyFailureCriterion"]=1
macorGlassCeramic["fractureEnergyReleaseRate"]=3.1223864e-05
macorGlassCeramic["fractureToughness"]=0.047434165

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=5 reflects relatively high scatter for machinable glass-ceramic strength governed by machining/surface flaws.
macorGlassCeramic["weibullReferenceVolume"] = 540.0
macorGlassCeramic["weibullModulus"] = 5.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
macorGlassCeramic["waveSpeed"] = float(np.sqrt( ( macorGlassCeramic["defaultBulkModulus"] + 4./3.*macorGlassCeramic["defaultShearModulus"] ) / macorGlassCeramic["defaultDensity"] ))

macorGlassCeramic["materialString"] = generateMaterialString(macorGlassCeramic)
# #################################################################################################

###################################################################################################
# MAGNESIUM OXIDE:
# Magnesium oxide: representative density and elastic properties from ceramic property tables;
# compression and K_IC are validation estimates.
#
magnesiumOxide = {}
magnesiumOxide["name"] = "magnesiumOxide"
magnesiumOxide["version"] = 2605160119
magnesiumOxide["model"] = "CeramicDamage"

# Constitutive model parameters:
magnesiumOxide["defaultDensity"]=3.58
magnesiumOxide["defaultBulkModulus"]=129.6875
magnesiumOxide["defaultShearModulus"]=105.50847

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
magnesiumOxide["tensileStrength"]=0.06
magnesiumOxide["compressiveStrength"]=1.0
magnesiumOxide["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
magnesiumOxide["crackSpeed"]=2.6070859

# Damaged/residual material behavior:
magnesiumOxide["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
magnesiumOxide["enableEnergyFailureCriterion"]=1
magnesiumOxide["fractureEnergyReleaseRate"]=8.7433735e-06
magnesiumOxide["fractureToughness"]=0.047434165

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume is the ASTM C1161 Type B full coupon volume fallback:
# 3 mm x 4 mm x 45 mm = 540 mm^3.
# m=5 reflects conservative scatter for MgO refractory ceramic strength.
magnesiumOxide["weibullReferenceVolume"] = 540.0
magnesiumOxide["weibullModulus"] = 5.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
magnesiumOxide["waveSpeed"] = float(np.sqrt( ( magnesiumOxide["defaultBulkModulus"] + 4./3.*magnesiumOxide["defaultShearModulus"] ) / magnesiumOxide["defaultDensity"] ))

magnesiumOxide["materialString"] = generateMaterialString(magnesiumOxide)
# #################################################################################################

###################################################################################################
# STABILIZED ZIRCONIA:
# Stabilized zirconia: representative structural-ceramic values similar to Y-TZP; tensile strength is
# a validation estimate.
#
stabilizedZirconia = {}
stabilizedZirconia["name"] = "stabilizedZirconia"
stabilizedZirconia["version"] = 2605160120
stabilizedZirconia["model"] = "CeramicDamage"

# Constitutive model parameters:
stabilizedZirconia["defaultDensity"]=6.0
stabilizedZirconia["defaultBulkModulus"]=175.0
stabilizedZirconia["defaultShearModulus"]=80.769231

# Effective flaw-controlled strengths, not intrinsic perfect-crystal strengths:
stabilizedZirconia["tensileStrength"]=0.55
stabilizedZirconia["compressiveStrength"]=2.485
stabilizedZirconia["maximumStrength"]=5.0

# Dynamic/time-to-failure style parameter:
stabilizedZirconia["crackSpeed"]=2.0592194

# Damaged/residual material behavior:
stabilizedZirconia["damagedMaterialFrictionSlope"]=0.25

# Energy regularization:
stabilizedZirconia["enableEnergyFailureCriterion"]=1
stabilizedZirconia["fractureEnergyReleaseRate"]=0.00043333334
stabilizedZirconia["fractureToughness"]=0.31622777

# Used by pfw to compute strength-scale for the material:
# Weibull reference volume uses a reported zirconia-family four-point-bend bar:
# 3 mm x 4 mm x 48 mm = 576 mm^3.
# m=15 reflects lower scatter/higher reliability typical of toughened/stabilized zirconias.
stabilizedZirconia["weibullReferenceVolume"] = 576.0
stabilizedZirconia["weibullModulus"] = 15.0

# Estimate wave speed for use in timing calculations.  This is not used by GEOS, but may be used
# in pfw, e.g. to make stopTime some multiple of the transit time for an elastic wave in the material
stabilizedZirconia["waveSpeed"] = float(np.sqrt( ( stabilizedZirconia["defaultBulkModulus"] + 4./3.*stabilizedZirconia["defaultShearModulus"] ) / stabilizedZirconia["defaultDensity"] ))

stabilizedZirconia["materialString"] = generateMaterialString(stabilizedZirconia)
# #################################################################################################

# -------------------------------------------------------------------------------------------------
# GRAPHITE AND ENGINEERING POLYMERS
# -------------------------------------------------------------------------------------------------
# Graphite entries use the transversely isotropic Graphite model parameters.  The elastic constants
# distinguish ideal single-crystal / pyrolytic graphite from a fine-grain near-isotropic engineering
# graphite.  Strength, fracture-energy, and plastic-response parameters are representative validation
# values chosen to be consistent with the model constraints and published engineering data; recalibrate
# them for production analysis against the exact graphite grade and stress-strain/fracture data.


def _set_graphite_auxiliary_properties(material):
    transverse_shear = material['defaultYoungModulusTransverse']/(2.0*(1.0 + material['defaultPoissonRatioTransverse']))
    material['defaultEffectiveBulkModulus'] = float((2.0*material['defaultYoungModulusTransverse'] + material['defaultYoungModulusAxial'])/9.0)
    material['defaultEffectiveShearModulus'] = float((2.0*material['defaultShearModulusAxialTransverse'] + transverse_shear)/3.0)
    material['weibullReferenceVolume'] = 1.0
    material['weibullModulus'] = 0.0
    material['waveSpeed'] = float(np.sqrt((material['defaultEffectiveBulkModulus'] + 4.0/3.0*material['defaultEffectiveShearModulus'])/material['defaultDensity']))


###################################################################################################
# GRAPHITE SINGLE CRYSTAL:
# Room-temperature single-crystal graphite approximation.  The transverse plane corresponds to the
# basal plane; the axial direction corresponds to the c axis.
#
graphiteSingleCrystal = {}
graphiteSingleCrystal['name'] = 'graphiteSingleCrystal'
graphiteSingleCrystal['version'] = 2605190001
graphiteSingleCrystal['model'] = 'Graphite'
graphiteSingleCrystal['defaultDensity'] = 2.267
graphiteSingleCrystal['defaultYoungModulusTransverse'] = 1090.0
graphiteSingleCrystal['defaultYoungModulusAxial'] = 38.7
graphiteSingleCrystal['defaultPoissonRatioTransverse'] = 0.125
graphiteSingleCrystal['defaultPoissonRatioAxialTransverse'] = 0.0
graphiteSingleCrystal['defaultShearModulusAxialTransverse'] = 4.95
graphiteSingleCrystal['defaultYoungModulusTransversePressureDerivative'] = 0.0
graphiteSingleCrystal['defaultYoungModulusAxialPressureDerivative'] = 0.0
graphiteSingleCrystal['defaultShearModulusAxialTransversePressureDerivative'] = 0.0
graphiteSingleCrystal['failureStrength'] = 0.035
graphiteSingleCrystal['maximumPrincipalStressDamage'] = 1
graphiteSingleCrystal['crackSpeed'] = 4.14
graphiteSingleCrystal['scaleFractureEnergyReleaseRate'] = 0
graphiteSingleCrystal['fractureEnergyStrengthScaleExponent'] = 0.0
graphiteSingleCrystal['basalPlaneFractureEnergyReleaseRate'] = 4.0e-7
graphiteSingleCrystal['totalFractureEnergyReleaseRate'] = 1.0e-5
graphiteSingleCrystal['damageEvolutionExponent'] = 32.0
graphiteSingleCrystal['damagedMaterialFrictionalSlope'] = 0.30
graphiteSingleCrystal['distortionShearResponseX2'] = 0.10
graphiteSingleCrystal['distortionShearResponseY1'] = 0.025
graphiteSingleCrystal['distortionShearResponseY2'] = 0.080
graphiteSingleCrystal['distortionShearResponseM1'] = 0.60
graphiteSingleCrystal['inPlaneShearResponseX2'] = 0.20
graphiteSingleCrystal['inPlaneShearResponseY1'] = 0.080
graphiteSingleCrystal['inPlaneShearResponseY2'] = 0.250
graphiteSingleCrystal['inPlaneShearResponseM1'] = 0.90
graphiteSingleCrystal['coupledShearResponseX2'] = 0.05
graphiteSingleCrystal['coupledShearResponseY1'] = 0.005
graphiteSingleCrystal['coupledShearResponseY2'] = 0.025
graphiteSingleCrystal['coupledShearResponseM1'] = 0.45
graphiteSingleCrystal['distortionStrainHardeningC0'] = 1.00
graphiteSingleCrystal['inPlaneStrainHardeningC0'] = 1.10
graphiteSingleCrystal['coupledStrainHardeningC0'] = 1.00
graphiteSingleCrystal['maximumPlasticStrain'] = 0.020
_set_graphite_auxiliary_properties(graphiteSingleCrystal)
# #################################################################################################

###################################################################################################
# PYROLYTIC GRAPHITE:
# Strongly anisotropic graphite approximation with slightly reduced basal-plane stiffness relative
# to the ideal single-crystal entry.
#
graphitePyrolytic = {}
graphitePyrolytic['name'] = 'graphitePyrolytic'
graphitePyrolytic['version'] = 2605190002
graphitePyrolytic['model'] = 'Graphite'
graphitePyrolytic['defaultDensity'] = 2.20
graphitePyrolytic['defaultYoungModulusTransverse'] = 900.0
graphitePyrolytic['defaultYoungModulusAxial'] = 36.5
graphitePyrolytic['defaultPoissonRatioTransverse'] = 0.13
graphitePyrolytic['defaultPoissonRatioAxialTransverse'] = 0.02
graphitePyrolytic['defaultShearModulusAxialTransverse'] = 4.50
graphitePyrolytic['defaultYoungModulusTransversePressureDerivative'] = 0.0
graphitePyrolytic['defaultYoungModulusAxialPressureDerivative'] = 0.0
graphitePyrolytic['defaultShearModulusAxialTransversePressureDerivative'] = 0.0
graphitePyrolytic['failureStrength'] = 0.030
graphitePyrolytic['maximumPrincipalStressDamage'] = 1
graphitePyrolytic['crackSpeed'] = 4.00
graphitePyrolytic['scaleFractureEnergyReleaseRate'] = 0
graphitePyrolytic['fractureEnergyStrengthScaleExponent'] = 0.0
graphitePyrolytic['basalPlaneFractureEnergyReleaseRate'] = 8.0e-7
graphitePyrolytic['totalFractureEnergyReleaseRate'] = 1.5e-5
graphitePyrolytic['damageEvolutionExponent'] = 32.0
graphitePyrolytic['damagedMaterialFrictionalSlope'] = 0.30
graphitePyrolytic['distortionShearResponseX2'] = 0.10
graphitePyrolytic['distortionShearResponseY1'] = 0.020
graphitePyrolytic['distortionShearResponseY2'] = 0.070
graphitePyrolytic['distortionShearResponseM1'] = 0.55
graphitePyrolytic['inPlaneShearResponseX2'] = 0.20
graphitePyrolytic['inPlaneShearResponseY1'] = 0.060
graphitePyrolytic['inPlaneShearResponseY2'] = 0.200
graphitePyrolytic['inPlaneShearResponseM1'] = 0.80
graphitePyrolytic['coupledShearResponseX2'] = 0.05
graphitePyrolytic['coupledShearResponseY1'] = 0.004
graphitePyrolytic['coupledShearResponseY2'] = 0.020
graphitePyrolytic['coupledShearResponseM1'] = 0.40
graphitePyrolytic['distortionStrainHardeningC0'] = 1.00
graphitePyrolytic['inPlaneStrainHardeningC0'] = 1.10
graphitePyrolytic['coupledStrainHardeningC0'] = 1.00
graphitePyrolytic['maximumPlasticStrain'] = 0.020
_set_graphite_auxiliary_properties(graphitePyrolytic)
# #################################################################################################

###################################################################################################
# PRESSURE-SENSITIVE GRAPHITIC CRYSTAL:
# Representative layered graphitic crystal parameterization using the signed-distortion
# strength split and saturating pressure-dependent elastic moduli.  The pressure scales below recover
# nearly linear stiffening when large and saturating stiffening when comparable to the pressure range.
# Fracture and post-failure parameters are placeholder validation values and should be calibrated for
# production simulations involving fracture propagation or fully comminuted flow.
#
graphitePressureSensitiveCrystal = {}
graphitePressureSensitiveCrystal['name'] = 'graphitePressureSensitiveCrystal'
graphitePressureSensitiveCrystal['version'] = 2606090002
graphitePressureSensitiveCrystal['model'] = 'Graphite'
graphitePressureSensitiveCrystal['defaultDensity'] = 1.50
graphitePressureSensitiveCrystal['defaultYoungModulusTransverse'] = 43.69846
graphitePressureSensitiveCrystal['defaultYoungModulusAxial'] = 15.46671
graphitePressureSensitiveCrystal['defaultPoissonRatioTransverse'] = 0.0
graphitePressureSensitiveCrystal['defaultPoissonRatioAxialTransverse'] = 0.2471611
graphitePressureSensitiveCrystal['defaultShearModulusAxialTransverse'] = 0.9568250
graphitePressureSensitiveCrystal['defaultYoungModulusTransversePressureDerivative'] = 18.31736
graphitePressureSensitiveCrystal['defaultYoungModulusAxialPressureDerivative'] = 9.316295
graphitePressureSensitiveCrystal['defaultShearModulusAxialTransversePressureDerivative'] = 0.5230170
graphitePressureSensitiveCrystal['defaultYoungModulusTransversePressureScale'] = 5.952545
graphitePressureSensitiveCrystal['defaultYoungModulusAxialPressureScale'] = 71.25111
graphitePressureSensitiveCrystal['defaultShearModulusAxialTransversePressureScale'] = 33.57501
graphitePressureSensitiveCrystal['failureStrength'] = 0.25
graphitePressureSensitiveCrystal['maximumPrincipalStressDamage'] = 1
graphitePressureSensitiveCrystal['crackSpeed'] = 4.0
graphitePressureSensitiveCrystal['scaleFractureEnergyReleaseRate'] = 0
graphitePressureSensitiveCrystal['fractureEnergyStrengthScaleExponent'] = 0.0
graphitePressureSensitiveCrystal['basalPlaneFractureEnergyReleaseRate'] = 1.0e-4
graphitePressureSensitiveCrystal['totalFractureEnergyReleaseRate'] = 1.0e-4
graphitePressureSensitiveCrystal['damageEvolutionExponent'] = 32.0
graphitePressureSensitiveCrystal['damagedMaterialFrictionalSlope'] = 0.20
graphitePressureSensitiveCrystal['distortionShearResponseX2'] = 31.726211302217223
graphitePressureSensitiveCrystal['distortionShearResponseY1'] = 1.1436120660995532
graphitePressureSensitiveCrystal['distortionShearResponseY2'] = 4.136203368288201
graphitePressureSensitiveCrystal['distortionShearResponseM1'] = 0.2
graphitePressureSensitiveCrystal['positiveDistortionShearResponseX2'] = 31.726211302217223
graphitePressureSensitiveCrystal['positiveDistortionShearResponseY1'] = 0.0
graphitePressureSensitiveCrystal['positiveDistortionShearResponseY2'] = 0.6093923874824587
graphitePressureSensitiveCrystal['positiveDistortionShearResponseM1'] = 0.2
graphitePressureSensitiveCrystal['inPlaneShearResponseX2'] = 31.726211302217223
graphitePressureSensitiveCrystal['inPlaneShearResponseY1'] = 0.1360119464627108
graphitePressureSensitiveCrystal['inPlaneShearResponseY2'] = 4.426039768025678
graphitePressureSensitiveCrystal['inPlaneShearResponseM1'] = 0.2
graphitePressureSensitiveCrystal['coupledShearResponseX2'] = 31.726211302217223
graphitePressureSensitiveCrystal['coupledShearResponseY1'] = 0.3855532087924691
graphitePressureSensitiveCrystal['coupledShearResponseY2'] = 3.7463651992725477
graphitePressureSensitiveCrystal['coupledShearResponseM1'] = 0.2
graphitePressureSensitiveCrystal['distortionStrainHardeningC0'] = 1.00
graphitePressureSensitiveCrystal['inPlaneStrainHardeningC0'] = 1.00
graphitePressureSensitiveCrystal['coupledStrainHardeningC0'] = 1.00
graphitePressureSensitiveCrystal['maximumPlasticStrain'] = 0.020
_set_graphite_auxiliary_properties(graphitePressureSensitiveCrystal)
# #################################################################################################

###################################################################################################
# FINE-GRAIN ISOTROPIC GRAPHITE:
# Engineering-grade near-isotropic graphite.  The Graphite model is still transversely isotropic, so
# the axial and transverse moduli are set close together rather than exactly equal.
#
graphiteFineGrainIso = {}
graphiteFineGrainIso['name'] = 'graphiteFineGrainIso'
graphiteFineGrainIso['version'] = 2605190003
graphiteFineGrainIso['model'] = 'Graphite'
graphiteFineGrainIso['defaultDensity'] = 1.82
graphiteFineGrainIso['defaultYoungModulusTransverse'] = 11.0
graphiteFineGrainIso['defaultYoungModulusAxial'] = 9.0
graphiteFineGrainIso['defaultPoissonRatioTransverse'] = 0.20
graphiteFineGrainIso['defaultPoissonRatioAxialTransverse'] = 0.15
graphiteFineGrainIso['defaultShearModulusAxialTransverse'] = 3.80
graphiteFineGrainIso['defaultYoungModulusTransversePressureDerivative'] = 0.0
graphiteFineGrainIso['defaultYoungModulusAxialPressureDerivative'] = 0.0
graphiteFineGrainIso['defaultShearModulusAxialTransversePressureDerivative'] = 0.0
graphiteFineGrainIso['failureStrength'] = 0.035
graphiteFineGrainIso['maximumPrincipalStressDamage'] = 1
graphiteFineGrainIso['crackSpeed'] = 2.50
graphiteFineGrainIso['scaleFractureEnergyReleaseRate'] = 0
graphiteFineGrainIso['fractureEnergyStrengthScaleExponent'] = 0.0
graphiteFineGrainIso['basalPlaneFractureEnergyReleaseRate'] = 2.0e-5
graphiteFineGrainIso['totalFractureEnergyReleaseRate'] = 3.0e-5
graphiteFineGrainIso['damageEvolutionExponent'] = 32.0
graphiteFineGrainIso['damagedMaterialFrictionalSlope'] = 0.35
graphiteFineGrainIso['distortionShearResponseX2'] = 0.10
graphiteFineGrainIso['distortionShearResponseY1'] = 0.035
graphiteFineGrainIso['distortionShearResponseY2'] = 0.120
graphiteFineGrainIso['distortionShearResponseM1'] = 0.85
graphiteFineGrainIso['inPlaneShearResponseX2'] = 0.10
graphiteFineGrainIso['inPlaneShearResponseY1'] = 0.040
graphiteFineGrainIso['inPlaneShearResponseY2'] = 0.130
graphiteFineGrainIso['inPlaneShearResponseM1'] = 0.90
graphiteFineGrainIso['coupledShearResponseX2'] = 0.08
graphiteFineGrainIso['coupledShearResponseY1'] = 0.025
graphiteFineGrainIso['coupledShearResponseY2'] = 0.080
graphiteFineGrainIso['coupledShearResponseM1'] = 0.70
graphiteFineGrainIso['distortionStrainHardeningC0'] = 1.10
graphiteFineGrainIso['inPlaneStrainHardeningC0'] = 1.10
graphiteFineGrainIso['coupledStrainHardeningC0'] = 1.05
graphiteFineGrainIso['maximumPlasticStrain'] = 0.030
_set_graphite_auxiliary_properties(graphiteFineGrainIso)
# #################################################################################################

# Polymer entries use the StrainHardeningPolymer model.  Published density, modulus, yield strength,
# and elongation data set the elastic constants, reference yield strength, and maximum stretch.  The
# model's zero-plastic-strain flow stress is defaultYieldStrength + shearSofteningMagnitude, so the
# defaultYieldStrength below is chosen to make that initial flow stress match the reference yield.
# The hardening and softening shape parameters are representative validation values and should be
# calibrated against full stress-strain curves for production studies.


def _set_polymer_elastic_constants(material, young_modulus, poisson_ratio):
    material['referenceYoungModulus'] = young_modulus
    material['referencePoissonRatio'] = poisson_ratio
    material['defaultBulkModulus'] = float(young_modulus/(3.0*(1.0 - 2.0*poisson_ratio)))
    material['defaultShearModulus'] = float(young_modulus/(2.0*(1.0 + poisson_ratio)))


def _set_polymer_reference_yield(material, reference_yield_strength):
    material['referenceYieldStrength'] = reference_yield_strength
    material['defaultYieldStrength'] = max(reference_yield_strength - material['shearSofteningMagnitude'], 0.0)


def _set_polymer_temperature_defaults(material):
    for prefix in ['yieldStrength', 'strainHardeningSlope', 'shearSofteningMagnitude', 'maximumStretch', 'bulkModulus', 'shearModulus']:
        material[prefix + 'A'] = 0.0
        material[prefix + 'B'] = 0.0
        material[prefix + 'T0'] = 300.0


def _set_polymer_wave_properties(material):
    material['weibullReferenceVolume'] = 1.0
    material['weibullModulus'] = 0.0
    material['waveSpeed'] = float(np.sqrt((material['defaultBulkModulus'] + 4.0/3.0*material['defaultShearModulus'])/material['defaultDensity']))



def _set_surface_polymer_model_defaults(material):
    # Parameters are in the MPM validation unit system: mm, us, mg, K, so stress is GPa.
    # The thermal scale is normalized to unity at glassTransitionTemperature.  Default values here
    # leave room-temperature response unchanged above the mechanical transition and stiffen only for
    # temperatures below the transition when a cold slope is supplied.
    material['hardeningScaleExponent'] = 1.0
    material['glassTransitionTemperature'] = 258.15
    material['temperatureColdSlope'] = 0.030
    material['temperatureHotSlope'] = 0.0
    material['temperatureTransitionMagnitude'] = 0.0
    material['temperatureTransitionWidth'] = 10.0
    material['crystallinity'] = 0.0
    material['referenceCrystallinity'] = 0.0
    material['crystallinityTransitionWidth'] = 10.0
    material['elasticCrystallinityCoeff'] = 0.0
    material['yieldStrengthCrystallinityCoeff'] = 0.0
    material['pressureAsymmetryAmplitude'] = 0.0
    material['pressureAsymmetryWidth'] = 10.0


###################################################################################################
# VITON/FKM 75 SHORE A FLUOROELASTOMER:
# Representative 75A fluoroelastomer parameterization for the SurfaceInformedPolymer model.  The
# elastic shear modulus is estimated from a published 100% modulus using an incompressible rubber
# relation for nominal tensile stress at lambda=2.  Strength and maximum stretch are selected from
# published tensile strength and elongation-at-break values for 75A FKM/Viton compound data.  The
# plastic softening and hardening coefficients are validation parameters chosen so that the reduced
# one-dimensional response passes through the same order of stress at 100% strain and at failure.
# Use batch-specific test data when available.
#
vitonFKM75SurfacePolymer = {}
vitonFKM75SurfacePolymer['name'] = 'vitonFKM75SurfacePolymer'
vitonFKM75SurfacePolymer['version'] = 2606091329
vitonFKM75SurfacePolymer['model'] = 'SurfaceInformedPolymer'
vitonFKM75SurfacePolymer['defaultDensity'] = 1.85
_set_polymer_elastic_constants(vitonFKM75SurfacePolymer, 0.01577, 0.49)
vitonFKM75SurfacePolymer['defaultDrainedLinearTEC'] = 1.8e-4
vitonFKM75SurfacePolymer['defaultYieldStrength'] = 0.0030
vitonFKM75SurfacePolymer['shearSofteningMagnitude'] = 0.0030
vitonFKM75SurfacePolymer['shearSofteningShapeParameter1'] = 0.30
vitonFKM75SurfacePolymer['shearSofteningShapeParameter2'] = 1.25
vitonFKM75SurfacePolymer['strainHardeningSlope'] = 0.0020
vitonFKM75SurfacePolymer['maximumStretch'] = 2.60
_set_surface_polymer_model_defaults(vitonFKM75SurfacePolymer)
_set_polymer_wave_properties(vitonFKM75SurfacePolymer)
# #################################################################################################

###################################################################################################
# VITON/FKM 75 SHORE A FLUOROELASTOMER COHESIVE ZONE:
# Thin-film cohesive projection of vitonFKM75SurfacePolymer.  The thickness is a validation default
# and should normally be changed to the physical film or bond-line thickness used in a model.
#
vitonFKM75SurfacePolymerCohesiveZone = {}
vitonFKM75SurfacePolymerCohesiveZone['name'] = 'vitonFKM75SurfacePolymerCohesiveZone'
vitonFKM75SurfacePolymerCohesiveZone['version'] = 2606091329
vitonFKM75SurfacePolymerCohesiveZone['model'] = 'SurfaceInformedPolymerCohesiveZone'
vitonFKM75SurfacePolymerCohesiveZone['thickness'] = 0.10
vitonFKM75SurfacePolymerCohesiveZone['bulkModulus'] = vitonFKM75SurfacePolymer['defaultBulkModulus']
vitonFKM75SurfacePolymerCohesiveZone['shearModulus'] = vitonFKM75SurfacePolymer['defaultShearModulus']
vitonFKM75SurfacePolymerCohesiveZone['defaultYieldStrength'] = vitonFKM75SurfacePolymer['defaultYieldStrength']
vitonFKM75SurfacePolymerCohesiveZone['shearSofteningMagnitude'] = vitonFKM75SurfacePolymer['shearSofteningMagnitude']
vitonFKM75SurfacePolymerCohesiveZone['shearSofteningShapeParameter1'] = vitonFKM75SurfacePolymer['shearSofteningShapeParameter1']
vitonFKM75SurfacePolymerCohesiveZone['shearSofteningShapeParameter2'] = vitonFKM75SurfacePolymer['shearSofteningShapeParameter2']
vitonFKM75SurfacePolymerCohesiveZone['strainHardeningSlope'] = vitonFKM75SurfacePolymer['strainHardeningSlope']
vitonFKM75SurfacePolymerCohesiveZone['hardeningScaleExponent'] = vitonFKM75SurfacePolymer['hardeningScaleExponent']
vitonFKM75SurfacePolymerCohesiveZone['maximumStretch'] = vitonFKM75SurfacePolymer['maximumStretch']
vitonFKM75SurfacePolymerCohesiveZone['glassTransitionTemperature'] = vitonFKM75SurfacePolymer['glassTransitionTemperature']
vitonFKM75SurfacePolymerCohesiveZone['temperatureColdSlope'] = vitonFKM75SurfacePolymer['temperatureColdSlope']
vitonFKM75SurfacePolymerCohesiveZone['temperatureHotSlope'] = vitonFKM75SurfacePolymer['temperatureHotSlope']
vitonFKM75SurfacePolymerCohesiveZone['temperatureTransitionMagnitude'] = vitonFKM75SurfacePolymer['temperatureTransitionMagnitude']
vitonFKM75SurfacePolymerCohesiveZone['temperatureTransitionWidth'] = vitonFKM75SurfacePolymer['temperatureTransitionWidth']
vitonFKM75SurfacePolymerCohesiveZone['crystallinity'] = vitonFKM75SurfacePolymer['crystallinity']
vitonFKM75SurfacePolymerCohesiveZone['referenceCrystallinity'] = vitonFKM75SurfacePolymer['referenceCrystallinity']
vitonFKM75SurfacePolymerCohesiveZone['crystallinityTransitionWidth'] = vitonFKM75SurfacePolymer['crystallinityTransitionWidth']
vitonFKM75SurfacePolymerCohesiveZone['elasticCrystallinityCoeff'] = vitonFKM75SurfacePolymer['elasticCrystallinityCoeff']
vitonFKM75SurfacePolymerCohesiveZone['yieldStrengthCrystallinityCoeff'] = vitonFKM75SurfacePolymer['yieldStrengthCrystallinityCoeff']
vitonFKM75SurfacePolymerCohesiveZone['pressureAsymmetryAmplitude'] = vitonFKM75SurfacePolymer['pressureAsymmetryAmplitude']
vitonFKM75SurfacePolymerCohesiveZone['pressureAsymmetryWidth'] = vitonFKM75SurfacePolymer['pressureAsymmetryWidth']
# #################################################################################################


###################################################################################################
# POLYCARBONATE:
# Ductile engineering thermoplastic with high impact resistance.
#
polymerPolycarbonate = {}
polymerPolycarbonate['name'] = 'polymerPolycarbonate'
polymerPolycarbonate['version'] = 2605190101
polymerPolycarbonate['model'] = 'StrainHardeningPolymer'
polymerPolycarbonate['defaultDensity'] = 1.20
_set_polymer_elastic_constants(polymerPolycarbonate, 2.35, 0.37)
polymerPolycarbonate['strainHardeningSlope'] = 0.012
polymerPolycarbonate['shearSofteningMagnitude'] = 0.010
polymerPolycarbonate['shearSofteningShapeParameter1'] = 0.10
polymerPolycarbonate['shearSofteningShapeParameter2'] = 1.10
polymerPolycarbonate['maximumStretch'] = 1.80
_set_polymer_reference_yield(polymerPolycarbonate, 0.060)
_set_polymer_temperature_defaults(polymerPolycarbonate)
_set_polymer_wave_properties(polymerPolycarbonate)
# #################################################################################################

###################################################################################################
# ABS:
# Acrylonitrile butadiene styrene engineering thermoplastic.
#
polymerABS = {}
polymerABS['name'] = 'polymerABS'
polymerABS['version'] = 2605190102
polymerABS['model'] = 'StrainHardeningPolymer'
polymerABS['defaultDensity'] = 1.05
_set_polymer_elastic_constants(polymerABS, 2.30, 0.35)
polymerABS['strainHardeningSlope'] = 0.010
polymerABS['shearSofteningMagnitude'] = 0.0075
polymerABS['shearSofteningShapeParameter1'] = 0.08
polymerABS['shearSofteningShapeParameter2'] = 1.00
polymerABS['maximumStretch'] = 1.24
_set_polymer_reference_yield(polymerABS, 0.045)
_set_polymer_temperature_defaults(polymerABS)
_set_polymer_wave_properties(polymerABS)
# #################################################################################################

###################################################################################################
# NYLON 6/6:
# Polyamide 6/6 engineering thermoplastic.
#
polymerNylon66 = {}
polymerNylon66['name'] = 'polymerNylon66'
polymerNylon66['version'] = 2605190103
polymerNylon66['model'] = 'StrainHardeningPolymer'
polymerNylon66['defaultDensity'] = 1.15
_set_polymer_elastic_constants(polymerNylon66, 2.93, 0.39)
polymerNylon66['strainHardeningSlope'] = 0.018
polymerNylon66['shearSofteningMagnitude'] = 0.010
polymerNylon66['shearSofteningShapeParameter1'] = 0.12
polymerNylon66['shearSofteningShapeParameter2'] = 1.00
polymerNylon66['maximumStretch'] = 1.50
_set_polymer_reference_yield(polymerNylon66, 0.083)
_set_polymer_temperature_defaults(polymerNylon66)
_set_polymer_wave_properties(polymerNylon66)
# #################################################################################################

###################################################################################################
# ACETAL / POM-C:
# Polyoxymethylene acetal copolymer.
#
polymerAcetalPOM = {}
polymerAcetalPOM['name'] = 'polymerAcetalPOM'
polymerAcetalPOM['version'] = 2605190104
polymerAcetalPOM['model'] = 'StrainHardeningPolymer'
polymerAcetalPOM['defaultDensity'] = 1.41
_set_polymer_elastic_constants(polymerAcetalPOM, 2.80, 0.35)
polymerAcetalPOM['strainHardeningSlope'] = 0.012
polymerAcetalPOM['shearSofteningMagnitude'] = 0.008
polymerAcetalPOM['shearSofteningShapeParameter1'] = 0.10
polymerAcetalPOM['shearSofteningShapeParameter2'] = 1.00
polymerAcetalPOM['maximumStretch'] = 1.32
_set_polymer_reference_yield(polymerAcetalPOM, 0.067)
_set_polymer_temperature_defaults(polymerAcetalPOM)
_set_polymer_wave_properties(polymerAcetalPOM)
# #################################################################################################

###################################################################################################
# PEEK:
# Polyether ether ketone high-performance engineering thermoplastic.
#
polymerPEEK = {}
polymerPEEK['name'] = 'polymerPEEK'
polymerPEEK['version'] = 2605190105
polymerPEEK['model'] = 'StrainHardeningPolymer'
polymerPEEK['defaultDensity'] = 1.32
_set_polymer_elastic_constants(polymerPEEK, 3.70, 0.37)
polymerPEEK['strainHardeningSlope'] = 0.025
polymerPEEK['shearSofteningMagnitude'] = 0.012
polymerPEEK['shearSofteningShapeParameter1'] = 0.12
polymerPEEK['shearSofteningShapeParameter2'] = 1.00
polymerPEEK['maximumStretch'] = 1.45
_set_polymer_reference_yield(polymerPEEK, 0.115)
_set_polymer_temperature_defaults(polymerPEEK)
_set_polymer_wave_properties(polymerPEEK)
# #################################################################################################

###################################################################################################
# PMMA:
# Acrylic, comparatively brittle polymer entry.
#
polymerPMMA = {}
polymerPMMA['name'] = 'polymerPMMA'
polymerPMMA['version'] = 2605190106
polymerPMMA['model'] = 'StrainHardeningPolymer'
polymerPMMA['defaultDensity'] = 1.18
_set_polymer_elastic_constants(polymerPMMA, 2.90, 0.35)
polymerPMMA['strainHardeningSlope'] = 0.004
polymerPMMA['shearSofteningMagnitude'] = 0.003
polymerPMMA['shearSofteningShapeParameter1'] = 0.04
polymerPMMA['shearSofteningShapeParameter2'] = 1.00
polymerPMMA['maximumStretch'] = 1.04
_set_polymer_reference_yield(polymerPMMA, 0.065)
_set_polymer_temperature_defaults(polymerPMMA)
_set_polymer_wave_properties(polymerPMMA)
# #################################################################################################

###################################################################################################
# HDPE:
# High-density polyethylene, highly ductile semicrystalline polyolefin.
#
polymerHDPE = {}
polymerHDPE['name'] = 'polymerHDPE'
polymerHDPE['version'] = 2605190107
polymerHDPE['model'] = 'StrainHardeningPolymer'
polymerHDPE['defaultDensity'] = 0.948
_set_polymer_elastic_constants(polymerHDPE, 0.86, 0.42)
polymerHDPE['strainHardeningSlope'] = 0.0010
polymerHDPE['shearSofteningMagnitude'] = 0.003
polymerHDPE['shearSofteningShapeParameter1'] = 0.25
polymerHDPE['shearSofteningShapeParameter2'] = 1.00
polymerHDPE['maximumStretch'] = 6.00
_set_polymer_reference_yield(polymerHDPE, 0.022)
_set_polymer_temperature_defaults(polymerHDPE)
_set_polymer_wave_properties(polymerHDPE)
# #################################################################################################

###################################################################################################
# POLYPROPYLENE:
# Ductile polypropylene homopolymer / general engineering polypropylene approximation.
#
polymerPolypropylene = {}
polymerPolypropylene['name'] = 'polymerPolypropylene'
polymerPolypropylene['version'] = 2605190108
polymerPolypropylene['model'] = 'StrainHardeningPolymer'
polymerPolypropylene['defaultDensity'] = 0.905
_set_polymer_elastic_constants(polymerPolypropylene, 1.50, 0.42)
polymerPolypropylene['strainHardeningSlope'] = 0.0025
polymerPolypropylene['shearSofteningMagnitude'] = 0.004
polymerPolypropylene['shearSofteningShapeParameter1'] = 0.18
polymerPolypropylene['shearSofteningShapeParameter2'] = 1.00
polymerPolypropylene['maximumStretch'] = 3.00
_set_polymer_reference_yield(polymerPolypropylene, 0.033)
_set_polymer_temperature_defaults(polymerPolypropylene)
_set_polymer_wave_properties(polymerPolypropylene)
# #################################################################################################

# -------------------------------------------------------------------------------------------------
# ADDITIONAL EXPLICIT-MPM SOLID MODEL EXAMPLES
# -------------------------------------------------------------------------------------------------
# These entries cover GEOS solid constitutive models that are explicitly dispatched by the MPM
# solver but were not represented in the original PFW material database.  They are intended as
# starting cards for verification/example inputs.  The Chiumenti entries use the same units as the
# rest of this file: density in mg/mm^3, stress in GPa, and fracture energy in GPa*mm
# (equivalently J/mm^2).


def _set_lame_wave_properties(material):
    bulk_modulus = material['defaultLambda'] + 2.0 * material['defaultShearModulus'] / 3.0
    material['referenceBulkModulus'] = bulk_modulus
    material['waveSpeed'] = float(np.sqrt((bulk_modulus + 4.0/3.0*material['defaultShearModulus']) / material['defaultDensity']))


def _set_isotropic_wave_properties(material):
    if 'defaultBulkModulus' in material and 'defaultShearModulus' in material:
        bulk_modulus = material['defaultBulkModulus']
        shear_modulus = material['defaultShearModulus']
    elif 'defaultYoungModulus' in material and 'defaultPoissonRatio' in material:
        bulk_modulus = material['defaultYoungModulus']/(3.0*(1.0 - 2.0*material['defaultPoissonRatio']))
        shear_modulus = material['defaultYoungModulus']/(2.0*(1.0 + material['defaultPoissonRatio']))
        material['referenceBulkModulus'] = bulk_modulus
        material['referenceShearModulus'] = shear_modulus
    else:
        return
    material['waveSpeed'] = float(np.sqrt((bulk_modulus + 4.0/3.0*shear_modulus) / material['defaultDensity']))


###################################################################################################
# HYPERELASTIC MMS - GENERALIZED VORTEX VERIFICATION:
# Manufactured-solution material used by the generalized-vortex MMS of Kamojjala et al.,
# "Verification tests in solid mechanics," Engineering with Computers 31:193-213 (2015),
# doi:10.1007/s00366-013-0342-x.  The paper lists rho=1000, E=1000, and nu=0.3,
# corresponding to lambda=577 and mu=385 after rounding.  The values below preserve the
# historical GEOS/PFW generalized-vortex input values used by this verification suite.
# This is a verification material card rather than a physical material recommendation.
#
hyperelasticMMS = {}
hyperelasticMMS['name'] = 'hyperelasticMMS'
hyperelasticMMS['version'] = 2605200001
hyperelasticMMS['model'] = 'HyperelasticMMS'
hyperelasticMMS['defaultDensity'] = 1000.0
hyperelasticMMS['defaultLambda'] = 577.0
hyperelasticMMS['defaultShearModulus'] = 384.615384615384585
_set_lame_wave_properties(hyperelasticMMS)
# #################################################################################################

###################################################################################################
# CHIUMENTI - HOMEL/HERBOLD KALTHOFF-WINKLER DFG PARAMETERIZATION:
# Homel and Herbold used K = 158.333 GPa, G = 73.0769 GPa, density = 8 g/cm^3,
# tensile strength = 570 MPa, fracture energy = 22.3 mJ/mm^2, and l_ch = 1.414 mm
# for their Kalthoff-Winkler DFG comparison.  The GEOS Chiumenti model multiplies
# criticalLength by the particle lengthScale, so this value should be revisited for
# a production mesh if particle volumes are not near 1 mm^3.
#
chiumentiHomelHerboldKalthoff = {}
chiumentiHomelHerboldKalthoff['name'] = 'chiumentiHomelHerboldKalthoff'
chiumentiHomelHerboldKalthoff['version'] = 2605200011
chiumentiHomelHerboldKalthoff['model'] = 'Chiumenti'
chiumentiHomelHerboldKalthoff['defaultDensity'] = 8.0
chiumentiHomelHerboldKalthoff['defaultBulkModulus'] = 158.333
chiumentiHomelHerboldKalthoff['defaultShearModulus'] = 73.0769
chiumentiHomelHerboldKalthoff['failureStrength'] = 0.570
chiumentiHomelHerboldKalthoff['energyReleaseRate'] = 0.0223
chiumentiHomelHerboldKalthoff['criticalLength'] = 1.414
_set_isotropic_wave_properties(chiumentiHomelHerboldKalthoff)
# #################################################################################################

###################################################################################################
# CHIUMENTI - HOMEL/HERBOLD CHARPY DFG PARAMETERIZATION:
# Same elastic constants and fracture energy as the Homel/Herbold Kalthoff-Winkler case, but with
# the simplified Charpy impact target tensile strength of 1.0 GPa from the DFG paper.
#
chiumentiHomelHerboldCharpy = {}
chiumentiHomelHerboldCharpy['name'] = 'chiumentiHomelHerboldCharpy'
chiumentiHomelHerboldCharpy['version'] = 2605200012
chiumentiHomelHerboldCharpy['model'] = 'Chiumenti'
chiumentiHomelHerboldCharpy['defaultDensity'] = 8.0
chiumentiHomelHerboldCharpy['defaultBulkModulus'] = 158.333
chiumentiHomelHerboldCharpy['defaultShearModulus'] = 73.0769
chiumentiHomelHerboldCharpy['failureStrength'] = 1.000
chiumentiHomelHerboldCharpy['energyReleaseRate'] = 0.0223
chiumentiHomelHerboldCharpy['criticalLength'] = 1.414
_set_isotropic_wave_properties(chiumentiHomelHerboldCharpy)
# #################################################################################################

###################################################################################################
# CHIUMENTI - CERVERA/CHIUMENTI 2006 BENCHMARK PARAMETERIZATION:
# Cervera and Chiumenti's 2006 mesh-objective tensile-cracking examples use E = 30 GPa,
# nu = 0.2, tensile strength = 2 MPa, and mode-I fracture energy = 100 J/m^2.
# Their examples use exponential softening; this GEOS card maps the same reference material
# properties onto the Chiumenti linear Rankine-damage update.  Density was not part of that
# quasi-static benchmark, so a concrete-like 2.4 mg/mm^3 is used as an explicit-MPM placeholder.
#
chiumentiCerveraChiumenti2006 = {}
chiumentiCerveraChiumenti2006['name'] = 'chiumentiCerveraChiumenti2006'
chiumentiCerveraChiumenti2006['version'] = 2605200013
chiumentiCerveraChiumenti2006['model'] = 'Chiumenti'
chiumentiCerveraChiumenti2006['defaultDensity'] = 2.4
chiumentiCerveraChiumenti2006['defaultYoungModulus'] = 30.0
chiumentiCerveraChiumenti2006['defaultPoissonRatio'] = 0.2
chiumentiCerveraChiumenti2006['failureStrength'] = 0.002
chiumentiCerveraChiumenti2006['energyReleaseRate'] = 1.0e-4
chiumentiCerveraChiumenti2006['criticalLength'] = 1.0
_set_isotropic_wave_properties(chiumentiCerveraChiumenti2006)
# #################################################################################################

###################################################################################################
# HYPERELASTIC NATURAL RUBBER:
# Compressible neo-Hookean-style Hyperelastic card for a common elastomer.  The nearly
# incompressible Poisson ratio is kept below GEOS's strict nu < 0.5 check.
#
hyperelasticNaturalRubber = {}
hyperelasticNaturalRubber['name'] = 'hyperelasticNaturalRubber'
hyperelasticNaturalRubber['version'] = 2605200021
hyperelasticNaturalRubber['model'] = 'Hyperelastic'
hyperelasticNaturalRubber['defaultDensity'] = 1.10
hyperelasticNaturalRubber['defaultYoungModulus'] = 0.010
hyperelasticNaturalRubber['defaultPoissonRatio'] = 0.49
_set_isotropic_wave_properties(hyperelasticNaturalRubber)
# #################################################################################################

###################################################################################################
# PERFECTLY PLASTIC MILD STEEL:
# Common isotropic elastic-perfectly-plastic steel approximation for the MPM PerfectlyPlastic model.
#
perfectlyPlasticMildSteel = {}
perfectlyPlasticMildSteel['name'] = 'perfectlyPlasticMildSteel'
perfectlyPlasticMildSteel['version'] = 2605200022
perfectlyPlasticMildSteel['model'] = 'PerfectlyPlastic'
perfectlyPlasticMildSteel['defaultDensity'] = 7.85
perfectlyPlasticMildSteel['defaultYoungModulus'] = 200.0
perfectlyPlasticMildSteel['defaultPoissonRatio'] = 0.26
perfectlyPlasticMildSteel['defaultYieldStress'] = 0.250
_set_isotropic_wave_properties(perfectlyPlasticMildSteel)
# #################################################################################################

###################################################################################################
# ELASTIC TRANSVERSE ISOTROPIC CARBON-FIBER/EPOXY:
# Unidirectional carbon-fiber/epoxy engineering-constant approximation for the MPM
# ElasticTransverseIsotropic model.  The axial direction follows the particle materialDirection.
#
transverseIsotropicCarbonFiberEpoxy = {}
transverseIsotropicCarbonFiberEpoxy['name'] = 'transverseIsotropicCarbonFiberEpoxy'
transverseIsotropicCarbonFiberEpoxy['version'] = 2605200023
transverseIsotropicCarbonFiberEpoxy['model'] = 'ElasticTransverseIsotropic'
transverseIsotropicCarbonFiberEpoxy['defaultDensity'] = 1.60
transverseIsotropicCarbonFiberEpoxy['defaultYoungModulusTransverse'] = 10.0
transverseIsotropicCarbonFiberEpoxy['defaultYoungModulusAxial'] = 135.0
transverseIsotropicCarbonFiberEpoxy['defaultPoissonRatioTransverse'] = 0.35
transverseIsotropicCarbonFiberEpoxy['defaultPoissonRatioAxialTransverse'] = 0.28
transverseIsotropicCarbonFiberEpoxy['defaultShearModulusAxialTransverse'] = 5.0
# Effective wave speed estimate only; GEOS computes effective moduli during model initialization.
transverseIsotropicCarbonFiberEpoxy['waveSpeed'] = float(np.sqrt((135.0 / (3.0*(1.0 - 2.0*0.28))) / transverseIsotropicCarbonFiberEpoxy['defaultDensity']))
# #################################################################################################

mpmExplicitSolidMaterials = [
    hyperelasticMMS,
    chiumentiHomelHerboldKalthoff,
    chiumentiHomelHerboldCharpy,
    chiumentiCerveraChiumenti2006,
    hyperelasticNaturalRubber,
    perfectlyPlasticMildSteel,
    transverseIsotropicCarbonFiberEpoxy,
]

engineeringMetals = [
    aluminum,
    al6061T6,
    al7075T6,
    steel,
    carbonSteelA36,
    alloySteel4140,
    stainlessSteel304,
    stainlessSteel316,
    toolSteelA2,
    maragingSteel300,
    grayCastIron,
    ductileIron,
    titaniumGrade2,
    ti64,
    copper,
    brassC260,
    phosphorBronzeC510,
    nickel,
    inconel718,
    cobaltChromeF75,
    magnesium,
    magnesiumAZ31B,
    zinc,
    zamak3,
    tungsten,
    molybdenum,
    tantalum,
    niobium,
    lead,
    tin,
]

engineeringCeramics = [
    alumina995,
    zirconiaYTZP,
    siliconCarbide,
    siliconNitride,
    boronCarbide,
    aluminumNitride,
    boronNitride,
    tungstenCarbide,
    titaniumCarbide,
    titaniumNitride,
    fusedSilica,
    quartz,
    sapphire,
    mullite,
    cordierite,
    steatite,
    porcelain,
    macorGlassCeramic,
    magnesiumOxide,
    stabilizedZirconia,
]

graphiteMaterials = [
    graphiteSingleCrystal,
    graphitePyrolytic,
    graphitePressureSensitiveCrystal,
    graphiteFineGrainIso,
]

engineeringPolymers = [
    vitonFKM75SurfacePolymer,
    polymerPolycarbonate,
    polymerABS,
    polymerNylon66,
    polymerAcetalPOM,
    polymerPEEK,
    polymerPMMA,
    polymerHDPE,
    polymerPolypropylene,
]

cohesiveZoneMaterials = [
    vitonFKM75SurfacePolymerCohesiveZone,
]

engineeringMaterials = engineeringMetals + engineeringCeramics + graphiteMaterials + engineeringPolymers

# -------------------------------------------------------------------------------------------------
# EXAMPLE-SUITE MATERIAL DICTIONARY ENTRIES
# -------------------------------------------------------------------------------------------------
# These entries are written in the same explicit dictionary style as the other material entries in
# this file.  They are used by the examples suite so pfw_input files can stay concise while still
# referring to named, reusable material parameterizations.

###################################################################################################
# ELASTIC DEMO MATERIAL:
# Small-unit-system ElasticIsotropic material for compact 2D demonstration problems such as
# elasticDisk and collidingDisks.  Units follow the MPM example convention: mm, us, mg, and GPa.
#
###################################################################################################
# ELASTIC DEMO MATERIAL:
# Generic ElasticIsotropic material for simple example problems.  This entry intentionally specifies
# exactly two independent elastic constants in the XML, defaultYoungModulus and defaultPoissonRatio,
# because current GEOS requires one valid elastic-constant pair and rejects over/under-specified
# ElasticIsotropic inputs.
#
elasticDemo = {}
elasticDemo["name"] = "elasticDemo"
elasticDemo["version"] = 2605170001
elasticDemo["model"] = "ElasticIsotropic"
elasticDemo["defaultDensity"] = 1.0
elasticDemo["defaultYoungModulus"] = 1.0
elasticDemo["defaultPoissonRatio"] = 0.25
elasticDemo["defaultBulkModulus"] = elasticDemo["defaultYoungModulus"]/(3.0*(1.0 - 2.0*elasticDemo["defaultPoissonRatio"]))
elasticDemo["defaultShearModulus"] = elasticDemo["defaultYoungModulus"]/(2.0*(1.0 + elasticDemo["defaultPoissonRatio"]))
elasticDemo["waveSpeed"] = float(np.sqrt((elasticDemo["defaultBulkModulus"] + 4.0/3.0*elasticDemo["defaultShearModulus"])/elasticDemo["defaultDensity"]))
elasticDemo["materialString"] = generateMaterialString(elasticDemo)
# #################################################################################################

# #################################################################################################

###################################################################################################
# ELASTIC ALUMINUM, SI-STYLE:
# ElasticIsotropic aluminum-like material for examples that already use SI-style density and moduli
# such as the bar and plate-impact examples.
#
###################################################################################################
# ELASTIC ALUMINUM MATERIAL:
# Reusable ElasticIsotropic aluminum-like material for example problems that only need an elastic
# solid.  Units follow the other MPM material examples: density in mg/mm^3 and stiffness in GPa.
# The XML again specifies exactly two independent constants, E and nu, to satisfy GEOS validation.
#
elasticAluminumSI = {}
elasticAluminumSI["name"] = "elasticAluminumSI"
elasticAluminumSI["version"] = 2605170002
elasticAluminumSI["model"] = "ElasticIsotropic"
elasticAluminumSI["defaultDensity"] = 2.70
elasticAluminumSI["defaultYoungModulus"] = 69.0
elasticAluminumSI["defaultPoissonRatio"] = 0.33
elasticAluminumSI["defaultBulkModulus"] = elasticAluminumSI["defaultYoungModulus"]/(3.0*(1.0 - 2.0*elasticAluminumSI["defaultPoissonRatio"]))
elasticAluminumSI["defaultShearModulus"] = elasticAluminumSI["defaultYoungModulus"]/(2.0*(1.0 + elasticAluminumSI["defaultPoissonRatio"]))
elasticAluminumSI["waveSpeed"] = float(np.sqrt((elasticAluminumSI["defaultBulkModulus"] + 4.0/3.0*elasticAluminumSI["defaultShearModulus"])/elasticAluminumSI["defaultDensity"]))
elasticAluminumSI["materialString"] = generateMaterialString(elasticAluminumSI)
# #################################################################################################

# #################################################################################################

###################################################################################################
# GHAREB GEOMECHANICS MATERIAL:
# Ghareb chalk parameterization for the borehole examples.  Units follow the existing MPM
# borehole examples: mm, microseconds, mg, and GPa.  The entry is intentionally written in
# the same explicit dictionary style as the other entries in this file.
#
# Notes on current GEOS inputs:
# - The current Geomechanics model requires the hardening-direction parameters dstrendh,
#   dfslopedh, dpeakI1dh, and dcrdh.  The historical examples predate those XML attributes.
# - strainHardeningK and strainHardeningN still define the scalar hardening variable.
#   Setting dstrendh = 1.0 preserves the older STREN_i + K*(1-exp(-N*gamma_p)) behavior.
# - The remaining hardening-direction parameters are zero so the historical examples do not
#   harden fSlope, peakI1, or cr unless a user edits this dictionary explicitly.
# - p0 corresponds to the default borehole confining pressure scale used by the example
#   input, max in compression around I1 = -3*p.
#
###################################################################################################
# GHAREB GEOMECHANICS MATERIAL:
# Ghareb chalk parameterization for the borehole examples.  Units follow the existing MPM borehole
# examples: mm, microseconds, mg, and GPa.  This is written as a normal explicit material dictionary,
# matching the style of the other entries in this file.
#
# Current GEOS Geomechanics validation checks require, among other things:
#   b0 > 0, b0+b1 > 0, g0 > 0, p0 < 0, p1 > 0, p3 > 0,
#   0 < cr < 1, fSlopeFailed <= fSlope, beta > 0,
#   fluidBulkModulus = fluidInitialPressure = 0,
#   enableCreep in {0,1}, and if creep is enabled then creepB > 0 and creepG > 0.
# The historical Ghareb example enabled creep but did not provide creepG, so this entry explicitly
# sets creepG = 1.0 to satisfy the current check while preserving the legacy creepB value.
#
ghareb = {}
ghareb["name"] = "ghareb"
ghareb["version"] = 2605170003
ghareb["model"] = "Geomechanics"
ghareb["defaultDensity"] = 1.57

# Nonlinear pressure-dependent bulk modulus parameters.
ghareb["b0"] = 1.67
ghareb["b1"] = 30.0
ghareb["b2"] = 0.3
ghareb["b3"] = 1.42
ghareb["b4"] = 0.015

# Hardening-direction parameters required by current GEOS.
ghareb["dstrendh"] = 1.0
ghareb["dfslopedh"] = 0.0
ghareb["dpeakI1dh"] = 0.0
ghareb["dcrdh"] = 0.0

# Pressure-dependent shear modulus / Poisson-ratio parameters.
ghareb["g1"] = 0.1
ghareb["g2"] = -0.0001
ghareb["g0"] = 1.5*ghareb["b0"]*(1.0 - 2.0*(ghareb["g1"] + ghareb["g2"]))/(1.0 + (ghareb["g1"] + ghareb["g2"]))
ghareb["g3"] = 0.3
ghareb["g4"] = 0.0

# Crush curve / cap parameters.
ghareb["p0"] = -0.030001
ghareb["p1"] = 0.75
ghareb["p2"] = 0.0
ghareb["p3"] = 0.255
ghareb["p4"] = 0.0
ghareb["cr"] = 0.2

# Fluid effects are disabled by current GEOS Geomechanics.
ghareb["fluidBulkModulus"] = 0.0
ghareb["fluidInitialPressure"] = 0.0

# Rate dependence terms are disabled in this example material.
ghareb["t1RateDependence"] = 0.0
ghareb["t2RateDependence"] = 0.0

# Yield-surface parameters.  These satisfy the current nonlinear Drucker-Prager limit-surface check:
# fSlope > ySlope > 0, stren > ySlope*peakI1, and peakI1 >= 0.
ghareb["peakI1"] = 0.0322/2.0
ghareb["fSlope"] = 0.18
ghareb["fSlopeFailed"] = ghareb["fSlope"] - 0.05
ghareb["stren"] = 0.50
ghareb["ySlope"] = 0.002
ghareb["beta"] = 1.4

# Damage parameters.
ghareb["fractureEnergyReleaseRate"] = 1.5e-8
ghareb["fractureSofteningExponent"] = 0.75
ghareb["fractureStress"] = 0.0184
ghareb["damageEvolutionCriterion"] = 1
ghareb["brittleDuctileTransition"] = 0.020

# Thermal/activation-energy parameters.  Q is zero here, but initialTemperature is kept positive.
ghareb["initialTemperature"] = 300.0
ghareb["Q"] = 0.0

# Buckling disabled, with valid placeholder parameters.
ghareb["enableBuckling"] = 0
ghareb["bucklingLength"] = 1.0
ghareb["bucklingAmplitude"] = 0.0

# Creep parameters from the historical Ghareb example plus the required positive creepG.
ghareb["enableCreep"] = 1
ghareb["creepC0"] = 10.0
ghareb["creepC1"] = 1.229
ghareb["creepC2"] = 7.45e-08
ghareb["creepA"] = 0.2205
ghareb["creepB"] = 0.000628
ghareb["creepC"] = 0.308e-10
ghareb["creepD"] = 2.307
ghareb["creepE"] = 0.0
ghareb["creepF"] = 0.564
ghareb["creepG"] = 1.0

# Deviatoric hardening law. dstrendh=1 above means this K value changes STREN directly.
ghareb["strainHardeningN"] = 250.0
ghareb["strainHardeningK"] = 0.017

# Constitutive-update controls.  These are normally optional in GEOS, but are written explicitly so
# the example remains robust if defaults change.
ghareb["plasticStrainTolerance"] = 1.0e-10
ghareb["stressReturnTolerance"] = 1.0e-6
ghareb["maxAllowedSubcycles"] = 256
ghareb["failedStepResponse"] = 2

ghareb["waveSpeed"] = float(np.sqrt((ghareb["b0"] + 4.0/3.0*ghareb["g0"])/ghareb["defaultDensity"]))

ghareb["materialString"] = generateMaterialString(ghareb)
# #################################################################################################

# #################################################################################################


###################################################################################################
# VERIFICATION ELASTIC MATERIAL:
# Simple dimensionless ElasticIsotropic entry for verification smoke/analytic cases.  The dictionary
# owns all parameters and materialString is consumed directly by PFW inputs, matching the style of
# other material database entries.
verificationElastic = {}
verificationElastic["name"] = "verificationElastic"
verificationElastic["version"] = 2605180001
verificationElastic["model"] = "ElasticIsotropic"
verificationElastic["defaultDensity"] = 1.0
verificationElastic["defaultYoungModulus"] = 1.0
verificationElastic["defaultPoissonRatio"] = 0.25
verificationElastic["defaultBulkModulus"] = verificationElastic["defaultYoungModulus"]/(3.0*(1.0 - 2.0*verificationElastic["defaultPoissonRatio"]))
verificationElastic["defaultShearModulus"] = verificationElastic["defaultYoungModulus"]/(2.0*(1.0 + verificationElastic["defaultPoissonRatio"]))
verificationElastic["waveSpeed"] = float(np.sqrt((verificationElastic["defaultBulkModulus"] + 4.0/3.0*verificationElastic["defaultShearModulus"])/verificationElastic["defaultDensity"]))
verificationElastic["materialString"] = generateMaterialString(verificationElastic)
# #################################################################################################

###################################################################################################
# VERIFICATION QUARTZ DAMAGE MATERIAL:
# Quartz is not a separate GEOS constitutive-model catalog entry in this MPM source tree; this is a
# deterministic quartz-like CeramicDamage card for verification cases that should not rely on
# stochastic Weibull fields unless the input explicitly creates them.
verificationQuartz = {}
verificationQuartz["name"] = "verificationQuartz"
verificationQuartz["version"] = 2605180003
verificationQuartz["model"] = "CeramicDamage"
verificationQuartz["defaultDensity"] = 2.65
verificationQuartz["defaultBulkModulus"] = 37.0
verificationQuartz["defaultShearModulus"] = 44.0
verificationQuartz["tensileStrength"] = 0.030
verificationQuartz["compressiveStrength"] = 0.300
verificationQuartz["maximumStrength"] = 5.0
verificationQuartz["crackSpeed"] = 1.0e16
verificationQuartz["damagedMaterialFrictionSlope"] = 0.25
verificationQuartz["enableEnergyFailureCriterion"] = 1
verificationQuartz["fractureEnergyReleaseRate"] = 1.0e-5
verificationQuartz["fractureToughness"] = 0.03
verificationQuartz["weibullReferenceVolume"] = 1.0
verificationQuartz["weibullModulus"] = 0.0
verificationQuartz["waveSpeed"] = float(np.sqrt((verificationQuartz["defaultBulkModulus"] + 4.0/3.0*verificationQuartz["defaultShearModulus"])/verificationQuartz["defaultDensity"]))
verificationQuartz["materialString"] = generateMaterialString(verificationQuartz)
# #################################################################################################


###################################################################################################
# VERIFICATION VON MISES MATERIAL:
# Small, dimensionless VonMisesJ card for the uniaxial plasticity verification.  It shares the
# verificationElastic stiffness so the expected elastic slope and yield point are simple, and uses a
# low yield strength so the fast single-block test reaches plastic flow at small strain.
verificationVonMises = {}
verificationVonMises["name"] = "verificationVonMises"
verificationVonMises["version"] = 2606050001
verificationVonMises["model"] = "VonMisesJ"
verificationVonMises["defaultDensity"] = 1.0
verificationVonMises["defaultYoungModulus"] = 1.0
verificationVonMises["defaultPoissonRatio"] = 0.25
verificationVonMises["defaultBulkModulus"] = verificationVonMises["defaultYoungModulus"]/(3.0*(1.0 - 2.0*verificationVonMises["defaultPoissonRatio"]))
verificationVonMises["defaultShearModulus"] = verificationVonMises["defaultYoungModulus"]/(2.0*(1.0 + verificationVonMises["defaultPoissonRatio"]))
verificationVonMises["defaultYieldStrength"] = 0.02
verificationVonMises["waveSpeed"] = float(np.sqrt((verificationVonMises["defaultBulkModulus"] + 4.0/3.0*verificationVonMises["defaultShearModulus"])/verificationVonMises["defaultDensity"]))
verificationVonMises["materialString"] = generateMaterialString(verificationVonMises)
# #################################################################################################



###################################################################################################
# VERIFICATION MATERIAL-SWAP MATERIALS:
# Dimensionless elastic cards with identical density and Poisson ratio but different stiffness.  The
# material-swap verification uses them to make the event visible in material type and stress/strain
# diagnostics without introducing a density-driven momentum jump.
verificationMaterialSwapSoft = {}
verificationMaterialSwapSoft["name"] = "verificationMaterialSwapSoft"
verificationMaterialSwapSoft["version"] = 2606050002
verificationMaterialSwapSoft["model"] = "ElasticIsotropic"
verificationMaterialSwapSoft["defaultDensity"] = 1.0
verificationMaterialSwapSoft["defaultYoungModulus"] = 1.0
verificationMaterialSwapSoft["defaultPoissonRatio"] = 0.25
verificationMaterialSwapSoft["defaultBulkModulus"] = verificationMaterialSwapSoft["defaultYoungModulus"]/(3.0*(1.0 - 2.0*verificationMaterialSwapSoft["defaultPoissonRatio"]))
verificationMaterialSwapSoft["defaultShearModulus"] = verificationMaterialSwapSoft["defaultYoungModulus"]/(2.0*(1.0 + verificationMaterialSwapSoft["defaultPoissonRatio"]))
verificationMaterialSwapSoft["waveSpeed"] = float(np.sqrt((verificationMaterialSwapSoft["defaultBulkModulus"] + 4.0/3.0*verificationMaterialSwapSoft["defaultShearModulus"])/verificationMaterialSwapSoft["defaultDensity"]))
verificationMaterialSwapSoft["materialString"] = generateMaterialString(verificationMaterialSwapSoft)

verificationMaterialSwapStiff = {}
verificationMaterialSwapStiff["name"] = "verificationMaterialSwapStiff"
verificationMaterialSwapStiff["version"] = 2606050003
verificationMaterialSwapStiff["model"] = "ElasticIsotropic"
verificationMaterialSwapStiff["defaultDensity"] = verificationMaterialSwapSoft["defaultDensity"]
verificationMaterialSwapStiff["defaultYoungModulus"] = 4.0 * verificationMaterialSwapSoft["defaultYoungModulus"]
verificationMaterialSwapStiff["defaultPoissonRatio"] = verificationMaterialSwapSoft["defaultPoissonRatio"]
verificationMaterialSwapStiff["defaultBulkModulus"] = verificationMaterialSwapStiff["defaultYoungModulus"]/(3.0*(1.0 - 2.0*verificationMaterialSwapStiff["defaultPoissonRatio"]))
verificationMaterialSwapStiff["defaultShearModulus"] = verificationMaterialSwapStiff["defaultYoungModulus"]/(2.0*(1.0 + verificationMaterialSwapStiff["defaultPoissonRatio"]))
verificationMaterialSwapStiff["waveSpeed"] = float(np.sqrt((verificationMaterialSwapStiff["defaultBulkModulus"] + 4.0/3.0*verificationMaterialSwapStiff["defaultShearModulus"])/verificationMaterialSwapStiff["defaultDensity"]))
verificationMaterialSwapStiff["materialString"] = generateMaterialString(verificationMaterialSwapStiff)
# #################################################################################################


###################################################################################################
# CONTACT SURFACE/GAP-CLOSURE HYPERELASTIC VERIFICATION MATERIAL:
# Dimensionless Neo-Hookean material used by verification/contactSurfaceGapClosure.  The values are
# intentionally soft so the compact quasistatic contact variants run quickly while retaining a
# clear uniaxial-strain reaction slope after the curved profiles mate.
contactGapClosureHyperelastic = {}
contactGapClosureHyperelastic["name"] = "contactGapClosureHyperelastic"
contactGapClosureHyperelastic["version"] = 2606050004
contactGapClosureHyperelastic["model"] = "HyperelasticMMS"
contactGapClosureHyperelastic["defaultDensity"] = 1.0
contactGapClosureHyperelastic["defaultYoungModulus"] = 10.0
contactGapClosureHyperelastic["defaultPoissonRatio"] = 0.25
contactGapClosureHyperelastic["defaultLambda"] = (contactGapClosureHyperelastic["defaultYoungModulus"] * contactGapClosureHyperelastic["defaultPoissonRatio"] /
                                                  ((1.0 + contactGapClosureHyperelastic["defaultPoissonRatio"]) *
                                                   (1.0 - 2.0 * contactGapClosureHyperelastic["defaultPoissonRatio"])))
contactGapClosureHyperelastic["defaultShearModulus"] = contactGapClosureHyperelastic["defaultYoungModulus"]/(2.0*(1.0 + contactGapClosureHyperelastic["defaultPoissonRatio"]))
contactGapClosureHyperelastic["waveSpeed"] = float(np.sqrt((contactGapClosureHyperelastic["defaultLambda"] + 2.0*contactGapClosureHyperelastic["defaultShearModulus"])/contactGapClosureHyperelastic["defaultDensity"]))
contactGapClosureHyperelastic["materialString"] = generateMaterialString(contactGapClosureHyperelastic)
# #################################################################################################


###################################################################################################
# PBC COMPACTION MOMENTUM VERIFICATION MATERIALS:
# Dimensionless, intentionally soft material cards used by verification/pbcCompactionMomentumTest.
# The elastic and plastic subcases share density, Young's modulus, and Poisson's ratio so that any
# difference in total x-momentum drift can be attributed to the constitutive update/contact state,
# not to a different elastic wave speed.
pbcCompactionElastic = {}
pbcCompactionElastic["name"] = "pbcCompactionElastic"
pbcCompactionElastic["version"] = 2606040001
pbcCompactionElastic["model"] = "ElasticIsotropic"
pbcCompactionElastic["defaultDensity"] = 1.0
pbcCompactionElastic["defaultYoungModulus"] = 10.0
pbcCompactionElastic["defaultPoissonRatio"] = 0.22
pbcCompactionElastic["defaultBulkModulus"] = pbcCompactionElastic["defaultYoungModulus"]/(3.0*(1.0 - 2.0*pbcCompactionElastic["defaultPoissonRatio"]))
pbcCompactionElastic["defaultShearModulus"] = pbcCompactionElastic["defaultYoungModulus"]/(2.0*(1.0 + pbcCompactionElastic["defaultPoissonRatio"]))
pbcCompactionElastic["waveSpeed"] = float(np.sqrt((pbcCompactionElastic["defaultBulkModulus"] + 4.0/3.0*pbcCompactionElastic["defaultShearModulus"])/pbcCompactionElastic["defaultDensity"]))
pbcCompactionElastic["materialString"] = generateMaterialString(pbcCompactionElastic)

pbcCompactionPlastic = {}
pbcCompactionPlastic["name"] = "pbcCompactionPlastic"
pbcCompactionPlastic["version"] = 2606040002
pbcCompactionPlastic["model"] = "VonMisesJ"
pbcCompactionPlastic["defaultDensity"] = pbcCompactionElastic["defaultDensity"]
pbcCompactionPlastic["defaultYoungModulus"] = pbcCompactionElastic["defaultYoungModulus"]
pbcCompactionPlastic["defaultPoissonRatio"] = pbcCompactionElastic["defaultPoissonRatio"]
pbcCompactionPlastic["defaultBulkModulus"] = pbcCompactionElastic["defaultBulkModulus"]
pbcCompactionPlastic["defaultShearModulus"] = pbcCompactionElastic["defaultShearModulus"]
pbcCompactionPlastic["defaultYieldStrength"] = 0.5
pbcCompactionPlastic["waveSpeed"] = pbcCompactionElastic["waveSpeed"]
pbcCompactionPlastic["materialString"] = generateMaterialString(pbcCompactionPlastic)
# #################################################################################################

# -------------------------------------------------------------------------------------------------
# FINALIZE MATERIAL DICTIONARIES
# -------------------------------------------------------------------------------------------------
# Regenerate each materialString from its model-specific generator and attach that generator directly
# to the dictionary-compatible material object.  This intentionally happens at the end of the module
# so every exported material entry uses the canonical generator for its current model.


def _finalize_material(material):
    generator = getMaterialStringGenerator(material)
    finalized = MaterialProperties(material)
    # Discard any pre-finalization XML before attaching the generator.  From this point on,
    # finalized["materialString"], finalized.get("materialString"), and
    # finalized.materialString are all generated from the current dictionary values.
    dict.pop(finalized, 'materialString', None)
    object.__setattr__(finalized, 'generateMaterialString', generator)
    object.__setattr__(finalized, '_autoRefreshMaterialString', True)
    finalized.refreshMaterialString()
    return finalized


def _material_variable_names():
    names = []
    for name, value in list(globals().items()):
        if isinstance(value, dict) and 'name' in value and 'model' in value:
            names.append(name)
    return names


def _refresh_material_list(material_list):
    refreshed = []
    for material in material_list:
        refreshed.append(globals()[material['name']])
    return refreshed


for _material_variable_name in _material_variable_names():
    globals()[_material_variable_name] = _finalize_material(globals()[_material_variable_name])

engineeringMetals = _refresh_material_list(engineeringMetals)
engineeringCeramics = _refresh_material_list(engineeringCeramics)
graphiteMaterials = _refresh_material_list(graphiteMaterials)
engineeringPolymers = _refresh_material_list(engineeringPolymers)
mpmExplicitSolidMaterials = _refresh_material_list(mpmExplicitSolidMaterials)
cohesiveZoneMaterials = _refresh_material_list(cohesiveZoneMaterials)
engineeringMaterials = engineeringMetals + engineeringCeramics + graphiteMaterials + engineeringPolymers

exampleSuiteMaterials = [
    elasticDemo,
    elasticAluminumSI,
    ghareb,
]

verificationMaterials = [
    verificationElastic,
    contactGapClosureHyperelastic,
    verificationMaterialSwapSoft,
    verificationMaterialSwapStiff,
    verificationQuartz,
    verificationVonMises,
    pbcCompactionElastic,
    pbcCompactionPlastic,
]

allMaterials = engineeringMaterials + mpmExplicitSolidMaterials + cohesiveZoneMaterials + exampleSuiteMaterials + verificationMaterials
materialDatabase = {material['name']: material for material in allMaterials}


def finalizeMaterialEntry(material):
    """Return a MaterialProperties entry with its generator attached and materialString refreshed."""
    return _finalize_material(material)


# Convenience camel-case alias for users who refer to the model name rather than the XML instance name.
hyperElasticMMS = hyperelasticMMS

try:
    del _material_variable_name
except NameError:
    pass

