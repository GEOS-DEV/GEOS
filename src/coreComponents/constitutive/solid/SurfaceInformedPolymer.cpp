/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SurfaceInformedPolymer.cpp
 */

#include "SurfaceInformedPolymer.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

SurfaceInformedPolymer::SurfaceInformedPolymer( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_deformationGradient(),
  m_plasticStrain(),
  m_equivalentPlasticStrain(),
  m_damage(),
  m_temperature(),
  m_jacobian(),
  m_yieldStrength(),
  m_defaultYieldStrength( DBL_MAX ),
  m_shearSofteningMagnitude( 0.0 ),
  m_shearSofteningShapeParameter1( 1.0 ),
  m_shearSofteningShapeParameter2( 1.0 ),
  m_strainHardeningSlope( 0.0 ),
  m_hardeningScaleExponent( 1.0 ),
  m_maximumStretch( DBL_MAX ),
  m_fractureStretchLambdaMin( 1.0 ),
  m_fractureStretchLambda0( 0.0 ),
  m_fractureStretchT0( 300.0 ),
  m_fractureStretchTemperatureScale( 1.0 ),
  m_glassTransitionTemperature( 300.0 ),
  m_temperatureColdSlope( 0.0 ),
  m_temperatureHotSlope( 0.0 ),
  m_temperatureTransitionMagnitude( 0.0 ),
  m_temperatureTransitionWidth( 1.0 ),
  m_crystallinity( 0.0 ),
  m_referenceCrystallinity( 0.0 ),
  m_crystallinityTransitionWidth( 1.0 ),
  m_elasticCrystallinityCoeff( 0.0 ),
  m_yieldStrengthCrystallinityCoeff( 0.0 ),
  m_pressureAsymmetryAmplitude( 0.0 ),
  m_pressureAsymmetryWidth( 1.0 ),
  m_compressivePressureStrengtheningCap( -1.0 )
{
  registerWrapper( viewKeyStruct::defaultYieldStrengthString(), &m_defaultYieldStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference yield strength at glassTransitionTemperature before softening and hardening contributions" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeString(), &m_shearSofteningMagnitude ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference magnitude of the decaying plastic softening term" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter1String(), &m_shearSofteningShapeParameter1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Equivalent plastic strain scale for the decaying softening term" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter2String(), &m_shearSofteningShapeParameter2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shape exponent for the decaying softening term" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeString(), &m_strainHardeningSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference stretch-hardening slope" );

  registerWrapper( viewKeyStruct::hardeningScaleExponentString(), &m_hardeningScaleExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent p_H in H(T)=H_ref*S_T(T)^p_H" );

  registerWrapper( viewKeyStruct::maximumStretchString(), &m_maximumStretch ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Constant chain stretch used to flag damage when fractureStretchLambda0 is non-positive" );

  registerWrapper( viewKeyStruct::fractureStretchLambdaMinString(), &m_fractureStretchLambdaMin ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum/asymptotic stretch lambda_min in the optional exponential fracture-stretch law" );

  registerWrapper( viewKeyStruct::fractureStretchLambda0String(), &m_fractureStretchLambda0 ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Scale lambda_0 in the optional exponential fracture-stretch law; non-positive values disable the law" );

  registerWrapper( viewKeyStruct::fractureStretchT0String(), &m_fractureStretchT0 ).
    setApplyDefaultValue( 300.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature T_0 in the optional exponential fracture-stretch law" );

  registerWrapper( viewKeyStruct::fractureStretchTemperatureScaleString(), &m_fractureStretchTemperatureScale ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature scale a in the optional exponential fracture-stretch law" );

  registerWrapper( viewKeyStruct::glassTransitionTemperatureString(), &m_glassTransitionTemperature ).
    setApplyDefaultValue( 300.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference transition temperature for the normalized thermal scale" );

  registerWrapper( viewKeyStruct::temperatureColdSlopeString(), &m_temperatureColdSlope ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Slope of log temperature scale below glassTransitionTemperature" );

  registerWrapper( viewKeyStruct::temperatureHotSlopeString(), &m_temperatureHotSlope ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Slope of log temperature scale above glassTransitionTemperature" );

  registerWrapper( viewKeyStruct::temperatureTransitionMagnitudeString(), &m_temperatureTransitionMagnitude ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Centered logistic drop magnitude in log temperature scale" );

  registerWrapper( viewKeyStruct::temperatureTransitionWidthString(), &m_temperatureTransitionWidth ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Width of the smooth temperature transition" );

  registerWrapper( viewKeyStruct::crystallinityString(), &m_crystallinity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Current crystallinity measure; units are user-defined but must match the coefficients" );

  registerWrapper( viewKeyStruct::referenceCrystallinityString(), &m_referenceCrystallinity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference crystallinity for tabulated model parameters" );

  registerWrapper( viewKeyStruct::crystallinityTransitionWidthString(), &m_crystallinityTransitionWidth ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature width used to activate crystallinity corrections" );

  registerWrapper( viewKeyStruct::elasticCrystallinityCoeffString(), &m_elasticCrystallinityCoeff ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Linear crystallinity coefficient for bulk and shear moduli" );

  registerWrapper( viewKeyStruct::yieldStrengthCrystallinityCoeffString(), &m_yieldStrengthCrystallinityCoeff ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Linear crystallinity coefficient for yield strength" );

  registerWrapper( viewKeyStruct::pressureAsymmetryAmplitudeString(), &m_pressureAsymmetryAmplitude ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Amplitude of the pressure-sensitive yield-strength asymmetry" );

  registerWrapper( viewKeyStruct::pressureAsymmetryWidthString(), &m_pressureAsymmetryWidth ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature width of the pressure-sensitive yield-strength asymmetry" );

  registerWrapper( viewKeyStruct::compressivePressureStrengtheningCapString(), &m_compressivePressureStrengtheningCap ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Magnitude of the compressive mean-stress cap applied only to the pressure-asymmetry term; a negative value disables the cap" );

  registerWrapper( viewKeyStruct::deformationGradientString(), &m_deformationGradient ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Particle deformation gradient" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Plastic strain in engineering Voigt notation" );

  registerWrapper( viewKeyStruct::equivalentPlasticStrainString(), &m_equivalentPlasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Equivalent plastic strain-like history variable" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Stretch-based damage flag" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 300.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Particle temperature used by the polymer scale functions" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Jacobian history updated from volumetric strain increments" );

  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Current scalar flow strength" );
}

SurfaceInformedPolymer::~SurfaceInformedPolymer()
{}

void SurfaceInformedPolymer::allocateConstitutiveData( dataRepository::Group & parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_deformationGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_equivalentPlasticStrain.resize( 0, numConstitutivePointsPerParentIndex );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_temperature.resize( 0 );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_yieldStrength.resize( 0 );
}

void SurfaceInformedPolymer::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_defaultBulkModulus <= 0.0, "Reference bulk modulus must be positive.", InputError );
  GEOS_THROW_IF( m_defaultShearModulus <= 0.0, "Reference shear modulus must be positive.", InputError );
  GEOS_THROW_IF( m_defaultYieldStrength < 0.0, "Reference yield strength must be nonnegative.", InputError );
  GEOS_THROW_IF( m_shearSofteningMagnitude < 0.0, "Shear softening magnitude must be nonnegative.", InputError );
  GEOS_THROW_IF( m_shearSofteningShapeParameter1 <= 0.0, "Shear softening shape parameter 1 must be positive.", InputError );
  GEOS_THROW_IF( m_shearSofteningShapeParameter2 <= 0.0, "Shear softening shape parameter 2 must be positive.", InputError );
  GEOS_THROW_IF( m_strainHardeningSlope < 0.0, "Strain hardening slope must be nonnegative.", InputError );
  GEOS_THROW_IF( m_maximumStretch <= 1.0, "Maximum stretch must be greater than one.", InputError );
  GEOS_THROW_IF( m_fractureStretchLambda0 < 0.0, "Fracture stretch lambda_0 must be nonnegative.", InputError );
  if( m_fractureStretchLambda0 > 0.0 )
  {
    GEOS_THROW_IF( m_fractureStretchLambdaMin <= 1.0, "Fracture stretch lambda_min must be greater than one when the exponential law is active.", InputError );
    GEOS_THROW_IF( m_fractureStretchTemperatureScale <= 0.0, "Fracture stretch temperature scale must be positive when the exponential law is active.", InputError );
  }
  GEOS_THROW_IF( m_temperatureTransitionWidth < 0.0, "Temperature transition width must be nonnegative.", InputError );
  GEOS_THROW_IF( m_crystallinityTransitionWidth < 0.0, "Crystallinity transition width must be nonnegative.", InputError );
  GEOS_THROW_IF( m_pressureAsymmetryWidth < 0.0, "Pressure asymmetry width must be nonnegative.", InputError );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::temperatureString() ).setApplyDefaultValue( m_glassTransitionTemperature );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::yieldStrengthString() ).setApplyDefaultValue( m_defaultYieldStrength );
}

void SurfaceInformedPolymer::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SurfaceInformedPolymer, std::string const &, Group * const )

} // namespace constitutive
} // namespace geos
