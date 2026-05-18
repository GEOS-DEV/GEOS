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

aluminum["materialString"]="<!--"+aluminum["name"]+" parameterization of "+aluminum["model"]+" model, version: "+str(aluminum["version"])+"""-->
<"""+aluminum["model"]+"""
name="""+'"'+aluminum["name"]+'"'+"""
defaultDensity=""" + '"' + str(aluminum["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(aluminum["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(aluminum["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(aluminum["defaultYieldStrength"]) + '"' + """
/>
"""
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

al6061T6["materialString"]="<!--"+al6061T6["name"]+" parameterization of "+al6061T6["model"]+" model, version: "+str(al6061T6["version"])+"""-->
<"""+al6061T6["model"]+"""
name="""+'"'+al6061T6["name"]+'"'+"""
defaultDensity=""" + '"' + str(al6061T6["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(al6061T6["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(al6061T6["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(al6061T6["defaultYieldStrength"]) + '"' + """
/>
"""
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

al7075T6["materialString"]="<!--"+al7075T6["name"]+" parameterization of "+al7075T6["model"]+" model, version: "+str(al7075T6["version"])+"""-->
<"""+al7075T6["model"]+"""
name="""+'"'+al7075T6["name"]+'"'+"""
defaultDensity=""" + '"' + str(al7075T6["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(al7075T6["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(al7075T6["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(al7075T6["defaultYieldStrength"]) + '"' + """
/>
"""
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

steel["materialString"]="<!--"+steel["name"]+" parameterization of "+steel["model"]+" model, version: "+str(steel["version"])+"""-->
<"""+steel["model"]+"""
name="""+'"'+steel["name"]+'"'+"""
defaultDensity=""" + '"' + str(steel["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(steel["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(steel["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(steel["defaultYieldStrength"]) + '"' + """
/>
"""
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

carbonSteelA36["materialString"]="<!--"+carbonSteelA36["name"]+" parameterization of "+carbonSteelA36["model"]+" model, version: "+str(carbonSteelA36["version"])+"""-->
<"""+carbonSteelA36["model"]+"""
name="""+'"'+carbonSteelA36["name"]+'"'+"""
defaultDensity=""" + '"' + str(carbonSteelA36["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(carbonSteelA36["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(carbonSteelA36["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(carbonSteelA36["defaultYieldStrength"]) + '"' + """
/>
"""
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

alloySteel4140["materialString"]="<!--"+alloySteel4140["name"]+" parameterization of "+alloySteel4140["model"]+" model, version: "+str(alloySteel4140["version"])+"""-->
<"""+alloySteel4140["model"]+"""
name="""+'"'+alloySteel4140["name"]+'"'+"""
defaultDensity=""" + '"' + str(alloySteel4140["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(alloySteel4140["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(alloySteel4140["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(alloySteel4140["defaultYieldStrength"]) + '"' + """
/>
"""
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

stainlessSteel304["materialString"]="<!--"+stainlessSteel304["name"]+" parameterization of "+stainlessSteel304["model"]+" model, version: "+str(stainlessSteel304["version"])+"""-->
<"""+stainlessSteel304["model"]+"""
name="""+'"'+stainlessSteel304["name"]+'"'+"""
defaultDensity=""" + '"' + str(stainlessSteel304["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(stainlessSteel304["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(stainlessSteel304["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(stainlessSteel304["defaultYieldStrength"]) + '"' + """
/>
"""
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

stainlessSteel316["materialString"]="<!--"+stainlessSteel316["name"]+" parameterization of "+stainlessSteel316["model"]+" model, version: "+str(stainlessSteel316["version"])+"""-->
<"""+stainlessSteel316["model"]+"""
name="""+'"'+stainlessSteel316["name"]+'"'+"""
defaultDensity=""" + '"' + str(stainlessSteel316["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(stainlessSteel316["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(stainlessSteel316["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(stainlessSteel316["defaultYieldStrength"]) + '"' + """
/>
"""
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

toolSteelA2["materialString"]="<!--"+toolSteelA2["name"]+" parameterization of "+toolSteelA2["model"]+" model, version: "+str(toolSteelA2["version"])+"""-->
<"""+toolSteelA2["model"]+"""
name="""+'"'+toolSteelA2["name"]+'"'+"""
defaultDensity=""" + '"' + str(toolSteelA2["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(toolSteelA2["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(toolSteelA2["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(toolSteelA2["defaultYieldStrength"]) + '"' + """
/>
"""
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

maragingSteel300["materialString"]="<!--"+maragingSteel300["name"]+" parameterization of "+maragingSteel300["model"]+" model, version: "+str(maragingSteel300["version"])+"""-->
<"""+maragingSteel300["model"]+"""
name="""+'"'+maragingSteel300["name"]+'"'+"""
defaultDensity=""" + '"' + str(maragingSteel300["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(maragingSteel300["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(maragingSteel300["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(maragingSteel300["defaultYieldStrength"]) + '"' + """
/>
"""
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

grayCastIron["materialString"]="<!--"+grayCastIron["name"]+" parameterization of "+grayCastIron["model"]+" model, version: "+str(grayCastIron["version"])+"""-->
<"""+grayCastIron["model"]+"""
name="""+'"'+grayCastIron["name"]+'"'+"""
defaultDensity=""" + '"' + str(grayCastIron["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(grayCastIron["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(grayCastIron["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(grayCastIron["defaultYieldStrength"]) + '"' + """
/>
"""
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

ductileIron["materialString"]="<!--"+ductileIron["name"]+" parameterization of "+ductileIron["model"]+" model, version: "+str(ductileIron["version"])+"""-->
<"""+ductileIron["model"]+"""
name="""+'"'+ductileIron["name"]+'"'+"""
defaultDensity=""" + '"' + str(ductileIron["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(ductileIron["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(ductileIron["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(ductileIron["defaultYieldStrength"]) + '"' + """
/>
"""
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

titaniumGrade2["materialString"]="<!--"+titaniumGrade2["name"]+" parameterization of "+titaniumGrade2["model"]+" model, version: "+str(titaniumGrade2["version"])+"""-->
<"""+titaniumGrade2["model"]+"""
name="""+'"'+titaniumGrade2["name"]+'"'+"""
defaultDensity=""" + '"' + str(titaniumGrade2["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(titaniumGrade2["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(titaniumGrade2["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(titaniumGrade2["defaultYieldStrength"]) + '"' + """
/>
"""
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

ti64["materialString"]="<!--"+ti64["name"]+" parameterization of "+ti64["model"]+" model, version: "+str(ti64["version"])+"""-->
<"""+ti64["model"]+"""
name="""+'"'+ti64["name"]+'"'+"""
defaultDensity=""" + '"' + str(ti64["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(ti64["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(ti64["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(ti64["defaultYieldStrength"]) + '"' + """
/>
"""
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

copper["materialString"]="<!--"+copper["name"]+" parameterization of "+copper["model"]+" model, version: "+str(copper["version"])+"""-->
<"""+copper["model"]+"""
name="""+'"'+copper["name"]+'"'+"""
defaultDensity=""" + '"' + str(copper["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(copper["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(copper["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(copper["defaultYieldStrength"]) + '"' + """
/>
"""
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

brassC260["materialString"]="<!--"+brassC260["name"]+" parameterization of "+brassC260["model"]+" model, version: "+str(brassC260["version"])+"""-->
<"""+brassC260["model"]+"""
name="""+'"'+brassC260["name"]+'"'+"""
defaultDensity=""" + '"' + str(brassC260["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(brassC260["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(brassC260["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(brassC260["defaultYieldStrength"]) + '"' + """
/>
"""
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

phosphorBronzeC510["materialString"]="<!--"+phosphorBronzeC510["name"]+" parameterization of "+phosphorBronzeC510["model"]+" model, version: "+str(phosphorBronzeC510["version"])+"""-->
<"""+phosphorBronzeC510["model"]+"""
name="""+'"'+phosphorBronzeC510["name"]+'"'+"""
defaultDensity=""" + '"' + str(phosphorBronzeC510["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(phosphorBronzeC510["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(phosphorBronzeC510["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(phosphorBronzeC510["defaultYieldStrength"]) + '"' + """
/>
"""
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

nickel["materialString"]="<!--"+nickel["name"]+" parameterization of "+nickel["model"]+" model, version: "+str(nickel["version"])+"""-->
<"""+nickel["model"]+"""
name="""+'"'+nickel["name"]+'"'+"""
defaultDensity=""" + '"' + str(nickel["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(nickel["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(nickel["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(nickel["defaultYieldStrength"]) + '"' + """
/>
"""
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

inconel718["materialString"]="<!--"+inconel718["name"]+" parameterization of "+inconel718["model"]+" model, version: "+str(inconel718["version"])+"""-->
<"""+inconel718["model"]+"""
name="""+'"'+inconel718["name"]+'"'+"""
defaultDensity=""" + '"' + str(inconel718["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(inconel718["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(inconel718["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(inconel718["defaultYieldStrength"]) + '"' + """
/>
"""
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

cobaltChromeF75["materialString"]="<!--"+cobaltChromeF75["name"]+" parameterization of "+cobaltChromeF75["model"]+" model, version: "+str(cobaltChromeF75["version"])+"""-->
<"""+cobaltChromeF75["model"]+"""
name="""+'"'+cobaltChromeF75["name"]+'"'+"""
defaultDensity=""" + '"' + str(cobaltChromeF75["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(cobaltChromeF75["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(cobaltChromeF75["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(cobaltChromeF75["defaultYieldStrength"]) + '"' + """
/>
"""
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

magnesium["materialString"]="<!--"+magnesium["name"]+" parameterization of "+magnesium["model"]+" model, version: "+str(magnesium["version"])+"""-->
<"""+magnesium["model"]+"""
name="""+'"'+magnesium["name"]+'"'+"""
defaultDensity=""" + '"' + str(magnesium["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(magnesium["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(magnesium["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(magnesium["defaultYieldStrength"]) + '"' + """
/>
"""
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

magnesiumAZ31B["materialString"]="<!--"+magnesiumAZ31B["name"]+" parameterization of "+magnesiumAZ31B["model"]+" model, version: "+str(magnesiumAZ31B["version"])+"""-->
<"""+magnesiumAZ31B["model"]+"""
name="""+'"'+magnesiumAZ31B["name"]+'"'+"""
defaultDensity=""" + '"' + str(magnesiumAZ31B["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(magnesiumAZ31B["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(magnesiumAZ31B["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(magnesiumAZ31B["defaultYieldStrength"]) + '"' + """
/>
"""
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

zinc["materialString"]="<!--"+zinc["name"]+" parameterization of "+zinc["model"]+" model, version: "+str(zinc["version"])+"""-->
<"""+zinc["model"]+"""
name="""+'"'+zinc["name"]+'"'+"""
defaultDensity=""" + '"' + str(zinc["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(zinc["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(zinc["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(zinc["defaultYieldStrength"]) + '"' + """
/>
"""
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

zamak3["materialString"]="<!--"+zamak3["name"]+" parameterization of "+zamak3["model"]+" model, version: "+str(zamak3["version"])+"""-->
<"""+zamak3["model"]+"""
name="""+'"'+zamak3["name"]+'"'+"""
defaultDensity=""" + '"' + str(zamak3["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(zamak3["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(zamak3["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(zamak3["defaultYieldStrength"]) + '"' + """
/>
"""
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

tungsten["materialString"]="<!--"+tungsten["name"]+" parameterization of "+tungsten["model"]+" model, version: "+str(tungsten["version"])+"""-->
<"""+tungsten["model"]+"""
name="""+'"'+tungsten["name"]+'"'+"""
defaultDensity=""" + '"' + str(tungsten["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(tungsten["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(tungsten["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(tungsten["defaultYieldStrength"]) + '"' + """
/>
"""
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

molybdenum["materialString"]="<!--"+molybdenum["name"]+" parameterization of "+molybdenum["model"]+" model, version: "+str(molybdenum["version"])+"""-->
<"""+molybdenum["model"]+"""
name="""+'"'+molybdenum["name"]+'"'+"""
defaultDensity=""" + '"' + str(molybdenum["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(molybdenum["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(molybdenum["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(molybdenum["defaultYieldStrength"]) + '"' + """
/>
"""
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

tantalum["materialString"]="<!--"+tantalum["name"]+" parameterization of "+tantalum["model"]+" model, version: "+str(tantalum["version"])+"""-->
<"""+tantalum["model"]+"""
name="""+'"'+tantalum["name"]+'"'+"""
defaultDensity=""" + '"' + str(tantalum["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(tantalum["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(tantalum["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(tantalum["defaultYieldStrength"]) + '"' + """
/>
"""
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

niobium["materialString"]="<!--"+niobium["name"]+" parameterization of "+niobium["model"]+" model, version: "+str(niobium["version"])+"""-->
<"""+niobium["model"]+"""
name="""+'"'+niobium["name"]+'"'+"""
defaultDensity=""" + '"' + str(niobium["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(niobium["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(niobium["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(niobium["defaultYieldStrength"]) + '"' + """
/>
"""
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

lead["materialString"]="<!--"+lead["name"]+" parameterization of "+lead["model"]+" model, version: "+str(lead["version"])+"""-->
<"""+lead["model"]+"""
name="""+'"'+lead["name"]+'"'+"""
defaultDensity=""" + '"' + str(lead["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(lead["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(lead["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(lead["defaultYieldStrength"]) + '"' + """
/>
"""
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

tin["materialString"]="<!--"+tin["name"]+" parameterization of "+tin["model"]+" model, version: "+str(tin["version"])+"""-->
<"""+tin["model"]+"""
name="""+'"'+tin["name"]+'"'+"""
defaultDensity=""" + '"' + str(tin["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(tin["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(tin["defaultShearModulus"]) + '"' + """
defaultYieldStrength=""" + '"' + str(tin["defaultYieldStrength"]) + '"' + """
/>
"""
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

alumina995["materialString"]="<!--"+alumina995["name"]+" parameterization of "+alumina995["model"]+" model, version: "+str(alumina995["version"])+"""-->
<"""+alumina995["model"]+"""
name="""+'"'+alumina995["name"]+'"'+"""
defaultDensity=""" + '"' + str(alumina995["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(alumina995["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(alumina995["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(alumina995["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(alumina995["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(alumina995["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(alumina995["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(alumina995["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(alumina995["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(alumina995["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(alumina995["fractureToughness"]) + '"' + """
/>
"""
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

zirconiaYTZP["materialString"]="<!--"+zirconiaYTZP["name"]+" parameterization of "+zirconiaYTZP["model"]+" model, version: "+str(zirconiaYTZP["version"])+"""-->
<"""+zirconiaYTZP["model"]+"""
name="""+'"'+zirconiaYTZP["name"]+'"'+"""
defaultDensity=""" + '"' + str(zirconiaYTZP["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(zirconiaYTZP["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(zirconiaYTZP["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(zirconiaYTZP["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(zirconiaYTZP["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(zirconiaYTZP["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(zirconiaYTZP["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(zirconiaYTZP["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(zirconiaYTZP["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(zirconiaYTZP["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(zirconiaYTZP["fractureToughness"]) + '"' + """
/>
"""
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

siliconCarbide["materialString"]="<!--"+siliconCarbide["name"]+" parameterization of "+siliconCarbide["model"]+" model, version: "+str(siliconCarbide["version"])+"""-->
<"""+siliconCarbide["model"]+"""
name="""+'"'+siliconCarbide["name"]+'"'+"""
defaultDensity=""" + '"' + str(siliconCarbide["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(siliconCarbide["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(siliconCarbide["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(siliconCarbide["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(siliconCarbide["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(siliconCarbide["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(siliconCarbide["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(siliconCarbide["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(siliconCarbide["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(siliconCarbide["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(siliconCarbide["fractureToughness"]) + '"' + """
/>
"""
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

siliconNitride["materialString"]="<!--"+siliconNitride["name"]+" parameterization of "+siliconNitride["model"]+" model, version: "+str(siliconNitride["version"])+"""-->
<"""+siliconNitride["model"]+"""
name="""+'"'+siliconNitride["name"]+'"'+"""
defaultDensity=""" + '"' + str(siliconNitride["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(siliconNitride["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(siliconNitride["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(siliconNitride["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(siliconNitride["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(siliconNitride["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(siliconNitride["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(siliconNitride["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(siliconNitride["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(siliconNitride["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(siliconNitride["fractureToughness"]) + '"' + """
/>
"""
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

boronCarbide["materialString"]="<!--"+boronCarbide["name"]+" parameterization of "+boronCarbide["model"]+" model, version: "+str(boronCarbide["version"])+"""-->
<"""+boronCarbide["model"]+"""
name="""+'"'+boronCarbide["name"]+'"'+"""
defaultDensity=""" + '"' + str(boronCarbide["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(boronCarbide["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(boronCarbide["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(boronCarbide["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(boronCarbide["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(boronCarbide["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(boronCarbide["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(boronCarbide["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(boronCarbide["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(boronCarbide["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(boronCarbide["fractureToughness"]) + '"' + """
/>
"""
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

aluminumNitride["materialString"]="<!--"+aluminumNitride["name"]+" parameterization of "+aluminumNitride["model"]+" model, version: "+str(aluminumNitride["version"])+"""-->
<"""+aluminumNitride["model"]+"""
name="""+'"'+aluminumNitride["name"]+'"'+"""
defaultDensity=""" + '"' + str(aluminumNitride["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(aluminumNitride["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(aluminumNitride["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(aluminumNitride["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(aluminumNitride["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(aluminumNitride["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(aluminumNitride["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(aluminumNitride["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(aluminumNitride["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(aluminumNitride["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(aluminumNitride["fractureToughness"]) + '"' + """
/>
"""
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

boronNitride["materialString"]="<!--"+boronNitride["name"]+" parameterization of "+boronNitride["model"]+" model, version: "+str(boronNitride["version"])+"""-->
<"""+boronNitride["model"]+"""
name="""+'"'+boronNitride["name"]+'"'+"""
defaultDensity=""" + '"' + str(boronNitride["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(boronNitride["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(boronNitride["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(boronNitride["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(boronNitride["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(boronNitride["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(boronNitride["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(boronNitride["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(boronNitride["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(boronNitride["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(boronNitride["fractureToughness"]) + '"' + """
/>
"""
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

tungstenCarbide["materialString"]="<!--"+tungstenCarbide["name"]+" parameterization of "+tungstenCarbide["model"]+" model, version: "+str(tungstenCarbide["version"])+"""-->
<"""+tungstenCarbide["model"]+"""
name="""+'"'+tungstenCarbide["name"]+'"'+"""
defaultDensity=""" + '"' + str(tungstenCarbide["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(tungstenCarbide["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(tungstenCarbide["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(tungstenCarbide["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(tungstenCarbide["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(tungstenCarbide["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(tungstenCarbide["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(tungstenCarbide["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(tungstenCarbide["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(tungstenCarbide["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(tungstenCarbide["fractureToughness"]) + '"' + """
/>
"""
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

titaniumCarbide["materialString"]="<!--"+titaniumCarbide["name"]+" parameterization of "+titaniumCarbide["model"]+" model, version: "+str(titaniumCarbide["version"])+"""-->
<"""+titaniumCarbide["model"]+"""
name="""+'"'+titaniumCarbide["name"]+'"'+"""
defaultDensity=""" + '"' + str(titaniumCarbide["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(titaniumCarbide["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(titaniumCarbide["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(titaniumCarbide["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(titaniumCarbide["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(titaniumCarbide["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(titaniumCarbide["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(titaniumCarbide["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(titaniumCarbide["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(titaniumCarbide["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(titaniumCarbide["fractureToughness"]) + '"' + """
/>
"""
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

titaniumNitride["materialString"]="<!--"+titaniumNitride["name"]+" parameterization of "+titaniumNitride["model"]+" model, version: "+str(titaniumNitride["version"])+"""-->
<"""+titaniumNitride["model"]+"""
name="""+'"'+titaniumNitride["name"]+'"'+"""
defaultDensity=""" + '"' + str(titaniumNitride["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(titaniumNitride["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(titaniumNitride["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(titaniumNitride["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(titaniumNitride["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(titaniumNitride["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(titaniumNitride["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(titaniumNitride["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(titaniumNitride["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(titaniumNitride["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(titaniumNitride["fractureToughness"]) + '"' + """
/>
"""
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

fusedSilica["materialString"]="<!--"+fusedSilica["name"]+" parameterization of "+fusedSilica["model"]+" model, version: "+str(fusedSilica["version"])+"""-->
<"""+fusedSilica["model"]+"""
name="""+'"'+fusedSilica["name"]+'"'+"""
defaultDensity=""" + '"' + str(fusedSilica["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(fusedSilica["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(fusedSilica["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(fusedSilica["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(fusedSilica["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(fusedSilica["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(fusedSilica["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(fusedSilica["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(fusedSilica["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(fusedSilica["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(fusedSilica["fractureToughness"]) + '"' + """
/>
"""
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

quartz["materialString"]="<!--"+quartz["name"]+" parameterization of "+quartz["model"]+" model, version: "+str(quartz["version"])+"""-->
<"""+quartz["model"]+"""
name="""+'"'+quartz["name"]+'"'+"""
defaultDensity=""" + '"' + str(quartz["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(quartz["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(quartz["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(quartz["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(quartz["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(quartz["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(quartz["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(quartz["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(quartz["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(quartz["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(quartz["fractureToughness"]) + '"' + """
/>
"""
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

sapphire["materialString"]="<!--"+sapphire["name"]+" parameterization of "+sapphire["model"]+" model, version: "+str(sapphire["version"])+"""-->
<"""+sapphire["model"]+"""
name="""+'"'+sapphire["name"]+'"'+"""
defaultDensity=""" + '"' + str(sapphire["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(sapphire["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(sapphire["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(sapphire["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(sapphire["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(sapphire["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(sapphire["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(sapphire["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(sapphire["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(sapphire["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(sapphire["fractureToughness"]) + '"' + """
/>
"""
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

mullite["materialString"]="<!--"+mullite["name"]+" parameterization of "+mullite["model"]+" model, version: "+str(mullite["version"])+"""-->
<"""+mullite["model"]+"""
name="""+'"'+mullite["name"]+'"'+"""
defaultDensity=""" + '"' + str(mullite["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(mullite["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(mullite["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(mullite["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(mullite["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(mullite["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(mullite["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(mullite["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(mullite["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(mullite["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(mullite["fractureToughness"]) + '"' + """
/>
"""
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

cordierite["materialString"]="<!--"+cordierite["name"]+" parameterization of "+cordierite["model"]+" model, version: "+str(cordierite["version"])+"""-->
<"""+cordierite["model"]+"""
name="""+'"'+cordierite["name"]+'"'+"""
defaultDensity=""" + '"' + str(cordierite["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(cordierite["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(cordierite["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(cordierite["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(cordierite["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(cordierite["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(cordierite["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(cordierite["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(cordierite["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(cordierite["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(cordierite["fractureToughness"]) + '"' + """
/>
"""
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

steatite["materialString"]="<!--"+steatite["name"]+" parameterization of "+steatite["model"]+" model, version: "+str(steatite["version"])+"""-->
<"""+steatite["model"]+"""
name="""+'"'+steatite["name"]+'"'+"""
defaultDensity=""" + '"' + str(steatite["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(steatite["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(steatite["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(steatite["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(steatite["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(steatite["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(steatite["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(steatite["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(steatite["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(steatite["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(steatite["fractureToughness"]) + '"' + """
/>
"""
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

porcelain["materialString"]="<!--"+porcelain["name"]+" parameterization of "+porcelain["model"]+" model, version: "+str(porcelain["version"])+"""-->
<"""+porcelain["model"]+"""
name="""+'"'+porcelain["name"]+'"'+"""
defaultDensity=""" + '"' + str(porcelain["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(porcelain["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(porcelain["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(porcelain["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(porcelain["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(porcelain["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(porcelain["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(porcelain["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(porcelain["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(porcelain["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(porcelain["fractureToughness"]) + '"' + """
/>
"""
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

macorGlassCeramic["materialString"]="<!--"+macorGlassCeramic["name"]+" parameterization of "+macorGlassCeramic["model"]+" model, version: "+str(macorGlassCeramic["version"])+"""-->
<"""+macorGlassCeramic["model"]+"""
name="""+'"'+macorGlassCeramic["name"]+'"'+"""
defaultDensity=""" + '"' + str(macorGlassCeramic["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(macorGlassCeramic["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(macorGlassCeramic["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(macorGlassCeramic["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(macorGlassCeramic["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(macorGlassCeramic["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(macorGlassCeramic["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(macorGlassCeramic["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(macorGlassCeramic["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(macorGlassCeramic["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(macorGlassCeramic["fractureToughness"]) + '"' + """
/>
"""
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

magnesiumOxide["materialString"]="<!--"+magnesiumOxide["name"]+" parameterization of "+magnesiumOxide["model"]+" model, version: "+str(magnesiumOxide["version"])+"""-->
<"""+magnesiumOxide["model"]+"""
name="""+'"'+magnesiumOxide["name"]+'"'+"""
defaultDensity=""" + '"' + str(magnesiumOxide["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(magnesiumOxide["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(magnesiumOxide["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(magnesiumOxide["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(magnesiumOxide["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(magnesiumOxide["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(magnesiumOxide["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(magnesiumOxide["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(magnesiumOxide["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(magnesiumOxide["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(magnesiumOxide["fractureToughness"]) + '"' + """
/>
"""
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

stabilizedZirconia["materialString"]="<!--"+stabilizedZirconia["name"]+" parameterization of "+stabilizedZirconia["model"]+" model, version: "+str(stabilizedZirconia["version"])+"""-->
<"""+stabilizedZirconia["model"]+"""
name="""+'"'+stabilizedZirconia["name"]+'"'+"""
defaultDensity=""" + '"' + str(stabilizedZirconia["defaultDensity"]) + '"' + """
defaultBulkModulus=""" + '"' + str(stabilizedZirconia["defaultBulkModulus"]) + '"' + """
defaultShearModulus=""" + '"' + str(stabilizedZirconia["defaultShearModulus"]) + '"' + """
tensileStrength=""" + '"' + str(stabilizedZirconia["tensileStrength"]) + '"' + """
compressiveStrength=""" + '"' + str(stabilizedZirconia["compressiveStrength"]) + '"' + """
maximumStrength=""" + '"' + str(stabilizedZirconia["maximumStrength"]) + '"' + """
crackSpeed=""" + '"' + str(stabilizedZirconia["crackSpeed"]) + '"' + """
damagedMaterialFrictionSlope=""" + '"' + str(stabilizedZirconia["damagedMaterialFrictionSlope"]) + '"' + """
enableEnergyFailureCriterion=""" + '"' + str(stabilizedZirconia["enableEnergyFailureCriterion"]) + '"' + """
fractureEnergyReleaseRate=""" + '"' + str(stabilizedZirconia["fractureEnergyReleaseRate"]) + '"' + """
fractureToughness=""" + '"' + str(stabilizedZirconia["fractureToughness"]) + '"' + """
/>
"""
# #################################################################################################

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

engineeringMaterials = engineeringMetals + engineeringCeramics

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
elasticDemo["materialString"] = "<!--"+elasticDemo["name"]+" parameterization of "+elasticDemo["model"]+" model, version: "+str(elasticDemo["version"])+"""-->
<"""+elasticDemo["model"]+"""
name="""+'"'+elasticDemo["name"]+'"'+"""
defaultDensity="""+'"'+str(elasticDemo["defaultDensity"])+'"'+"""
defaultYoungModulus="""+'"'+str(elasticDemo["defaultYoungModulus"])+'"'+"""
defaultPoissonRatio="""+'"'+str(elasticDemo["defaultPoissonRatio"])+'"'+"""
/>
"""
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
elasticAluminumSI["materialString"] = "<!--"+elasticAluminumSI["name"]+" parameterization of "+elasticAluminumSI["model"]+" model, version: "+str(elasticAluminumSI["version"])+"""-->
<"""+elasticAluminumSI["model"]+"""
name="""+'"'+elasticAluminumSI["name"]+'"'+"""
defaultDensity="""+'"'+str(elasticAluminumSI["defaultDensity"])+'"'+"""
defaultYoungModulus="""+'"'+str(elasticAluminumSI["defaultYoungModulus"])+'"'+"""
defaultPoissonRatio="""+'"'+str(elasticAluminumSI["defaultPoissonRatio"])+'"'+"""
/>
"""
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

ghareb["materialString"] = "<!--"+ghareb["name"]+" parameterization of "+ghareb["model"]+" model, version: "+str(ghareb["version"])+"""-->
<"""+ghareb["model"]+"""
name="""+'"'+ghareb["name"]+'"'+"""
defaultDensity="""+'"'+str(ghareb["defaultDensity"])+'"'+"""
b0="""+'"'+str(ghareb["b0"])+'"'+"""
b1="""+'"'+str(ghareb["b1"])+'"'+"""
b2="""+'"'+str(ghareb["b2"])+'"'+"""
b3="""+'"'+str(ghareb["b3"])+'"'+"""
b4="""+'"'+str(ghareb["b4"])+'"'+"""
dstrendh="""+'"'+str(ghareb["dstrendh"])+'"'+"""
dfslopedh="""+'"'+str(ghareb["dfslopedh"])+'"'+"""
dpeakI1dh="""+'"'+str(ghareb["dpeakI1dh"])+'"'+"""
dcrdh="""+'"'+str(ghareb["dcrdh"])+'"'+"""
g0="""+'"'+str(ghareb["g0"])+'"'+"""
g1="""+'"'+str(ghareb["g1"])+'"'+"""
g2="""+'"'+str(ghareb["g2"])+'"'+"""
g3="""+'"'+str(ghareb["g3"])+'"'+"""
g4="""+'"'+str(ghareb["g4"])+'"'+"""
p0="""+'"'+str(ghareb["p0"])+'"'+"""
p1="""+'"'+str(ghareb["p1"])+'"'+"""
p2="""+'"'+str(ghareb["p2"])+'"'+"""
p3="""+'"'+str(ghareb["p3"])+'"'+"""
p4="""+'"'+str(ghareb["p4"])+'"'+"""
cr="""+'"'+str(ghareb["cr"])+'"'+"""
fluidBulkModulus="""+'"'+str(ghareb["fluidBulkModulus"])+'"'+"""
fluidInitialPressure="""+'"'+str(ghareb["fluidInitialPressure"])+'"'+"""
t1RateDependence="""+'"'+str(ghareb["t1RateDependence"])+'"'+"""
t2RateDependence="""+'"'+str(ghareb["t2RateDependence"])+'"'+"""
peakI1="""+'"'+str(ghareb["peakI1"])+'"'+"""
fSlope="""+'"'+str(ghareb["fSlope"])+'"'+"""
fSlopeFailed="""+'"'+str(ghareb["fSlopeFailed"])+'"'+"""
stren="""+'"'+str(ghareb["stren"])+'"'+"""
ySlope="""+'"'+str(ghareb["ySlope"])+'"'+"""
beta="""+'"'+str(ghareb["beta"])+'"'+"""
fractureEnergyReleaseRate="""+'"'+str(ghareb["fractureEnergyReleaseRate"])+'"'+"""
fractureSofteningExponent="""+'"'+str(ghareb["fractureSofteningExponent"])+'"'+"""
fractureStress="""+'"'+str(ghareb["fractureStress"])+'"'+"""
initialTemperature="""+'"'+str(ghareb["initialTemperature"])+'"'+"""
Q="""+'"'+str(ghareb["Q"])+'"'+"""
brittleDuctileTransition="""+'"'+str(ghareb["brittleDuctileTransition"])+'"'+"""
damageEvolutionCriterion="""+'"'+str(ghareb["damageEvolutionCriterion"])+'"'+"""
enableBuckling="""+'"'+str(ghareb["enableBuckling"])+'"'+"""
bucklingLength="""+'"'+str(ghareb["bucklingLength"])+'"'+"""
bucklingAmplitude="""+'"'+str(ghareb["bucklingAmplitude"])+'"'+"""
enableCreep="""+'"'+str(ghareb["enableCreep"])+'"'+"""
creepC0="""+'"'+str(ghareb["creepC0"])+'"'+"""
creepC1="""+'"'+str(ghareb["creepC1"])+'"'+"""
creepC2="""+'"'+str(ghareb["creepC2"])+'"'+"""
creepA="""+'"'+str(ghareb["creepA"])+'"'+"""
creepB="""+'"'+str(ghareb["creepB"])+'"'+"""
creepC="""+'"'+str(ghareb["creepC"])+'"'+"""
creepD="""+'"'+str(ghareb["creepD"])+'"'+"""
creepE="""+'"'+str(ghareb["creepE"])+'"'+"""
creepF="""+'"'+str(ghareb["creepF"])+'"'+"""
creepG="""+'"'+str(ghareb["creepG"])+'"'+"""
strainHardeningN="""+'"'+str(ghareb["strainHardeningN"])+'"'+"""
strainHardeningK="""+'"'+str(ghareb["strainHardeningK"])+'"'+"""
plasticStrainTolerance="""+'"'+str(ghareb["plasticStrainTolerance"])+'"'+"""
stressReturnTolerance="""+'"'+str(ghareb["stressReturnTolerance"])+'"'+"""
maxAllowedSubcycles="""+'"'+str(ghareb["maxAllowedSubcycles"])+'"'+"""
failedStepResponse="""+'"'+str(ghareb["failedStepResponse"])+'"'+"""
/>
"""
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
verificationElastic["materialString"] = (
    "<!--"+verificationElastic["name"]+" parameterization of "+verificationElastic["model"]+" model, version: "+str(verificationElastic["version"])+"-->\n"
    "<"+verificationElastic["model"]+"\n"
    "name=\""+verificationElastic["name"]+"\"\n"
    "defaultDensity=\""+str(verificationElastic["defaultDensity"])+"\"\n"
    "defaultYoungModulus=\""+str(verificationElastic["defaultYoungModulus"])+"\"\n"
    "defaultPoissonRatio=\""+str(verificationElastic["defaultPoissonRatio"])+"\"\n"
    "/>\n"
)
# #################################################################################################

###################################################################################################
# VERIFICATION QUARTZ DAMAGE MATERIAL:
# Deterministic quartz-like damage material for verification cases that should not rely on stochastic
# Weibull fields unless the input explicitly creates them.
verificationQuartz = {}
verificationQuartz["name"] = "verificationQuartz"
verificationQuartz["version"] = 2605180002
verificationQuartz["model"] = "Quartz"
verificationQuartz["defaultDensity"] = 2.65
verificationQuartz["defaultBulkModulus"] = 37.0
verificationQuartz["defaultShearModulus"] = 44.0
verificationQuartz["defaultTensileStrength"] = 0.030
verificationQuartz["defaultCompressiveStrength"] = 0.300
verificationQuartz["defaultShearStrength"] = 0.050
verificationQuartz["fractureEnergy"] = 1.0e-5
verificationQuartz["weibullReferenceVolume"] = 1.0
verificationQuartz["weibullModulus"] = 0.0
verificationQuartz["waveSpeed"] = float(np.sqrt((verificationQuartz["defaultBulkModulus"] + 4.0/3.0*verificationQuartz["defaultShearModulus"])/verificationQuartz["defaultDensity"]))
verificationQuartz["materialString"] = (
    "<!--"+verificationQuartz["name"]+" parameterization of "+verificationQuartz["model"]+" model, version: "+str(verificationQuartz["version"])+"-->\n"
    "<"+verificationQuartz["model"]+"\n"
    "name=\""+verificationQuartz["name"]+"\"\n"
    "defaultDensity=\""+str(verificationQuartz["defaultDensity"])+"\"\n"
    "defaultBulkModulus=\""+str(verificationQuartz["defaultBulkModulus"])+"\"\n"
    "defaultShearModulus=\""+str(verificationQuartz["defaultShearModulus"])+"\"\n"
    "defaultTensileStrength=\""+str(verificationQuartz["defaultTensileStrength"])+"\"\n"
    "defaultCompressiveStrength=\""+str(verificationQuartz["defaultCompressiveStrength"])+"\"\n"
    "defaultShearStrength=\""+str(verificationQuartz["defaultShearStrength"])+"\"\n"
    "fractureEnergy=\""+str(verificationQuartz["fractureEnergy"])+"\"\n"
    "/>\n"
)
# #################################################################################################

