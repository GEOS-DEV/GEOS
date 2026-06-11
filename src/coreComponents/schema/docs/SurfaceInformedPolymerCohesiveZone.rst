
================================= ========= ======== ===============================================================================
Name                              Type      Default  Description
================================= ========= ======== ===============================================================================
bulkModulus                       real64    required Reference bulk modulus at glassTransitionTemperature
crystallinity                     real64    0        Current crystallinity measure used by the optional crystallinity multipliers
compressivePressureStrengtheningCap             real64    -1       Magnitude of the compressive mean-stress cap applied to the pressure-asymmetry term
defaultYieldStrength              real64    required Reference yield strength at glassTransitionTemperature
elasticCrystallinityCoeff         real64    0        Linear crystallinity coefficient for bulk and shear moduli
glassTransitionTemperature        real64    300      Reference transition temperature for the normalized thermal scale
hardeningScaleExponent            real64    1        Exponent p_H in H(T)=H_ref*S_T(T)^p_H
maximumStretch                    real64    DBL_MAX  Stretch at which the cohesive law returns zero traction
name                              groupName required A name is required for any non-unique nodes
pressureAsymmetryAmplitude        real64    0        Amplitude of the pressure-sensitive yield-strength asymmetry
pressureAsymmetryWidth            real64    1        Temperature width of the pressure-sensitive yield-strength asymmetry
referenceCrystallinity            real64    0        Reference crystallinity for the tabulated material parameters
shearModulus                      real64    required Reference shear modulus at glassTransitionTemperature
shearSofteningMagnitude           real64    required Reference magnitude of the decaying plastic softening term
shearSofteningShapeParameter1     real64    required Equivalent plastic strain scale for the softening term
shearSofteningShapeParameter2     real64    required Shape exponent for the softening term
strainHardeningSlope              real64    required Reference stretch-hardening slope
temperatureColdSlope              real64    0        Slope of log temperature scale below glassTransitionTemperature
temperatureHotSlope               real64    0        Slope of log temperature scale above glassTransitionTemperature
temperatureTransitionMagnitude    real64    0        Centered logistic drop magnitude in log temperature scale
temperatureTransitionWidth        real64    1        Width of the smooth temperature transition
thickness                         real64    required Cohesive film thickness
yieldStrengthCrystallinityCoeff   real64    0        Linear crystallinity coefficient for yield strength
================================= ========= ======== ===============================================================================
