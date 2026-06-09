
================================= ========= ======== ===============================================================================
Name                              Type      Default  Description
================================= ========= ======== ===============================================================================
crystallinity                     real64    0        Current crystallinity measure used by the optional crystallinity multipliers
defaultBulkModulus                real64    -1       Reference bulk modulus at glassTransitionTemperature
defaultDensity                    real64    required Default material density
defaultDrainedLinearTEC           real64    0        Default linear thermal expansion coefficient
defaultPoissonRatio               real64    -1       Default Poisson ratio
defaultShearModulus               real64    -1       Reference shear modulus at glassTransitionTemperature
defaultYieldStrength              real64    required Reference yield strength at glassTransitionTemperature
defaultYoungModulus               real64    -1       Default Young modulus
elasticCrystallinityCoeff         real64    0        Linear crystallinity coefficient for bulk and shear moduli
glassTransitionTemperature        real64    300      Reference transition temperature for the normalized thermal scale
hardeningScaleExponent            real64    1        Exponent p_H in H(T)=H_ref*S_T(T)^p_H
maximumStretch                    real64    DBL_MAX  Maximum principal stretch used to flag damage
name                              groupName required A name is required for any non-unique nodes
pressureAsymmetryAmplitude        real64    0        Amplitude of the pressure-sensitive yield-strength asymmetry
pressureAsymmetryWidth            real64    1        Temperature width of the pressure-sensitive yield-strength asymmetry
referenceCrystallinity            real64    0        Reference crystallinity for the tabulated material parameters
shearSofteningMagnitude           real64    required Reference magnitude of the decaying plastic softening term
shearSofteningShapeParameter1     real64    required Equivalent plastic strain scale for the softening term
shearSofteningShapeParameter2     real64    required Shape exponent for the softening term
strainHardeningSlope              real64    required Reference stretch-hardening slope
temperatureColdSlope              real64    0        Slope of log temperature scale below glassTransitionTemperature
temperatureHotSlope               real64    0        Slope of log temperature scale above glassTransitionTemperature
temperatureTransitionMagnitude    real64    0        Centered logistic drop magnitude in log temperature scale
temperatureTransitionWidth        real64    1        Width of the smooth temperature transition
yieldStrengthCrystallinityCoeff   real64    0        Linear crystallinity coefficient for yield strength
================================= ========= ======== ===============================================================================
