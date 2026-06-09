/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file Graphite.cpp
 */

#include "Graphite.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Graphite::Graphite( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultYoungModulusTransverse(),
  m_defaultYoungModulusAxial(),
  m_defaultPoissonRatioTransverse(),
  m_defaultPoissonRatioAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_defaultYoungModulusTransversePressureDerivative(),
  m_defaultYoungModulusAxialPressureDerivative(),
  m_defaultShearModulusAxialTransversePressureDerivative(),
  m_velocityGradient(),
  m_plasticStrain(),
  m_relaxation(),
  m_enableCrackTipStressConcentration(),
  m_crackTipStressConcentration(),
  m_distanceToCrackTip(),
  m_basalPlanePlasticWork(),
  m_plasticWork(),
  m_alphaL(),
  m_alphaT(),
  m_damage(),
  m_basalPlaneDamage(),
  m_comminutionDamage(),
  m_temperature(),
  m_temperatureRate(),
  m_jacobian(),
  m_lengthScale(),
  m_strengthScale(),
  m_failureStrength(),
  m_maximumPrincipalStressDamage(),
  m_crackSpeed(),
  m_scaleFractureEnergyReleaseRate(),
  m_fractureEnergyStrengthScaleExponent(),
  m_basalPlaneFractureEnergyReleaseRate(),
  m_totalFractureEnergyReleaseRate(),
  m_damageEvolutionExponent(),
  m_damagedMaterialFrictionalSlope(),
  m_distortionShearResponseX2(),
  m_distortionShearResponseY1(),
  m_distortionShearResponseY2(),
  m_distortionShearResponseM1(),
  m_positiveDistortionShearResponseX2(),
  m_positiveDistortionShearResponseY1(),
  m_positiveDistortionShearResponseY2(),
  m_positiveDistortionShearResponseM1(),
  m_inPlaneShearResponseX2(),
  m_inPlaneShearResponseY1(),
  m_inPlaneShearResponseY2(),
  m_inPlaneShearResponseM1(),
  m_coupledShearResponseX2(),
  m_coupledShearResponseY1(),
  m_coupledShearResponseY2(),
  m_coupledShearResponseM1(),
  m_distortionStrainHardeningC0(),
  m_inPlaneStrainHardeningC0(),
  m_coupledStrainHardeningC0(),
  m_maximumPlasticStrain()
{
  // register default values
  registerWrapper( viewKeyStruct::defaultYoungModulusTransverseString(), &m_defaultYoungModulusTransverse ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Transverse Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialString(), &m_defaultYoungModulusAxial ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Axial Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverseString(), &m_defaultPoissonRatioTransverse ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverseString(), &m_defaultPoissonRatioAxialTransverse ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Axial-Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverseString(), &m_defaultShearModulusAxialTransverse ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Axial-Transverse Shear Modulus" );

  registerWrapper( viewKeyStruct::defaultYoungModulusTransversePressureDerivativeString(), &m_defaultYoungModulusTransversePressureDerivative ).
    setApplyDefaultValue( 0. ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Transverse Young's modulus pressure derivative" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialPressureDerivativeString(), &m_defaultYoungModulusAxialPressureDerivative ).
    setApplyDefaultValue( 0. ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Axial Young's modulus pressure derivative" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransversePressureDerivativeString(), &m_defaultShearModulusAxialTransversePressureDerivative ).
    setApplyDefaultValue( 0. ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Axial transverse shear modulus pressure derivative" );

  registerWrapper( viewKeyStruct::failureStrengthString(), &m_failureStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum theoretical strength" );

  registerWrapper( viewKeyStruct::maximumPrincipalStressDamageString(), &m_maximumPrincipalStressDamage ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to apply Max principal Stress Failure / Damage criterion" );

  registerWrapper( viewKeyStruct::crackSpeedString(), &m_crackSpeed ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Crack speed for brittle failure" );

  registerWrapper( viewKeyStruct::scaleFractureEnergyReleaseRateString(), &m_scaleFractureEnergyReleaseRate ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Legacy flag used to initialize fractureEnergyStrengthScaleExponent when that exponent is not supplied" );

  registerWrapper( viewKeyStruct::fractureEnergyStrengthScaleExponentString(), &m_fractureEnergyStrengthScaleExponent ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent eta_G in effective fracture energy G_f^eff = G_f * strengthScale^eta_G; negative default derives eta_G from scaleFractureEnergyReleaseRate" );

  registerWrapper( viewKeyStruct::basalPlaneFractureEnergyReleaseRateString(), &m_basalPlaneFractureEnergyReleaseRate ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fracture energy release rates for basal shear and normal modes" );

  registerWrapper( viewKeyStruct::totalFractureEnergyReleaseRateString(), &m_totalFractureEnergyReleaseRate ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fracture Energy Release Rate for non-basal shear and normal modes" );

  registerWrapper( viewKeyStruct::damageEvolutionExponentString(), &m_damageEvolutionExponent ).
    setApplyDefaultValue( 32.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent used in work-based damage evolution D = min(1, W/(G_f/l))^n" );

  registerWrapper( viewKeyStruct::damagedMaterialFrictionalSlopeString(), &m_damagedMaterialFrictionalSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Damaged material frictional slope" );

  registerWrapper( viewKeyStruct::distortionShearResponseX2String(), &m_distortionShearResponseX2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Negative signed-distortion response X2 for I_d^- = max(-I_d,0)" );

  registerWrapper( viewKeyStruct::distortionShearResponseY1String(), &m_distortionShearResponseY1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Negative signed-distortion response Y1 for I_d^- = max(-I_d,0)" );

  registerWrapper( viewKeyStruct::distortionShearResponseY2String(), &m_distortionShearResponseY2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Negative signed-distortion response Y2 for I_d^- = max(-I_d,0)" );

  registerWrapper( viewKeyStruct::distortionShearResponseM1String(), &m_distortionShearResponseM1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Negative signed-distortion response M1 for I_d^- = max(-I_d,0)" );

  registerWrapper( viewKeyStruct::positiveDistortionShearResponseX2String(), &m_positiveDistortionShearResponseX2 ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional positive signed-distortion response X2 for I_d^+ = max(I_d,0); negative default copies distortionShearResponseX2" );

  registerWrapper( viewKeyStruct::positiveDistortionShearResponseY1String(), &m_positiveDistortionShearResponseY1 ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional positive signed-distortion response Y1 for I_d^+ = max(I_d,0); negative default copies distortionShearResponseY1" );

  registerWrapper( viewKeyStruct::positiveDistortionShearResponseY2String(), &m_positiveDistortionShearResponseY2 ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional positive signed-distortion response Y2 for I_d^+ = max(I_d,0); negative default copies distortionShearResponseY2" );

  registerWrapper( viewKeyStruct::positiveDistortionShearResponseM1String(), &m_positiveDistortionShearResponseM1 ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional positive signed-distortion response M1 for I_d^+ = max(I_d,0); negative default copies distortionShearResponseM1" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseX2String(), &m_inPlaneShearResponseX2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response X2" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseY1String(), &m_inPlaneShearResponseY1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response Y1" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseY2String(), &m_inPlaneShearResponseY2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response Y2" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseM1String(), &m_inPlaneShearResponseM1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response M1" );

  registerWrapper( viewKeyStruct::coupledShearResponseX2String(), &m_coupledShearResponseX2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response X2" );

  registerWrapper( viewKeyStruct::coupledShearResponseY1String(), &m_coupledShearResponseY1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response Y1" );

  registerWrapper( viewKeyStruct::coupledShearResponseY2String(), &m_coupledShearResponseY2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response Y2" );

  registerWrapper( viewKeyStruct::coupledShearResponseM1String(), &m_coupledShearResponseM1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response M1" );

  registerWrapper( viewKeyStruct::distortionStrainHardeningC0(), &m_distortionStrainHardeningC0 ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Distortion Strain Hardening Multiplier C0" );

  registerWrapper( viewKeyStruct::inPlaneStrainHardeningC0(), &m_inPlaneStrainHardeningC0 ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "In Plane train Hardening Multiplier C0" );

  registerWrapper( viewKeyStruct::coupledStrainHardeningC0(), &m_coupledStrainHardeningC0 ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coupled Strain Hardening Multiplier C0" );

  registerWrapper( viewKeyStruct::maximumPlasticStrainString(), &m_maximumPlasticStrain ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum plastic strain" );

  // register fields
  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Plastic strain" );

  registerWrapper( viewKeyStruct::relaxationString(), &m_relaxation ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Relaxation" );

  registerWrapper( viewKeyStruct::enableCrackTipStressConcentrationString(), &m_enableCrackTipStressConcentration ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Use crack-tip stress concentration" );

  registerWrapper( viewKeyStruct::crackTipStressConcentrationString(), &m_crackTipStressConcentration ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Crack tip stress concentration" );

  registerWrapper( viewKeyStruct::distanceToCrackTipString(), &m_distanceToCrackTip ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Distance to crack tip" );

  registerWrapper( viewKeyStruct::basalPlanePlasticWorkString(), &m_basalPlanePlasticWork ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Basal Plane Plastic Work" );

  registerWrapper( viewKeyStruct::plasticWorkString(), &m_plasticWork ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Plastic Work" );

  registerWrapper( viewKeyStruct::alphaLString(), &m_alphaL ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "constant for thermal expansion longitudinal to symmetry axis" );

  registerWrapper( viewKeyStruct::alphaTString(), &m_alphaT ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "constant for thermal expansion transverse to symmetry axis" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Host-visible DFG damage/failure indicator" );

  registerWrapper( viewKeyStruct::basalPlaneDamageString(), &m_basalPlaneDamage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Internal basal-plane fracture damage used to degrade basal opening/sliding strength" );

  registerWrapper( viewKeyStruct::comminutionDamageString(), &m_comminutionDamage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Internal comminution/powder damage used to transition strengths to residual powder response" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 300.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point temperature values" );

  registerWrapper( viewKeyStruct::temperatureRateString(), &m_temperatureRate ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point temperatureRate values" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point jacobian values" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point length scale values" );

  registerWrapper( viewKeyStruct::strengthScaleString(), &m_strengthScale ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point strength scale values" );

  registerWrapper( viewKeyStruct::effectiveBulkModulusString(), &m_effectiveBulkModulus ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Effective bulk modulus for stress control and wavespeed calculations" );

  registerWrapper( viewKeyStruct::effectiveShearModulusString(), &m_effectiveShearModulus ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Effective shear modulus for stress control and wavespeed calculations" );

  registerWrapper( viewKeyStruct::materialDirectionString(), &m_materialDirection ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Material direction - first row is used as graphite basal plane normal" );
}


Graphite::~Graphite()
{}


void Graphite::allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_effectiveBulkModulus.resize( 0 );
  m_effectiveShearModulus.resize( 0 );
  m_materialDirection.resize( 0, 3, 3 );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_relaxation.resize( 0, numConstitutivePointsPerParentIndex );
  m_crackTipStressConcentration.resize( 0 );
  m_distanceToCrackTip.resize( 0 );
  m_basalPlanePlasticWork.resize( 0, numConstitutivePointsPerParentIndex );
  m_plasticWork.resize( 0, numConstitutivePointsPerParentIndex );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_basalPlaneDamage.resize( 0, numConstitutivePointsPerParentIndex );
  m_comminutionDamage.resize( 0, numConstitutivePointsPerParentIndex );
  m_temperature.resize( 0 );
  m_temperatureRate.resize( 0 );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_lengthScale.resize( 0 );
  m_strengthScale.resize( 0 );
}


void Graphite::postInputInitialization()
{
  SolidBase::postInputInitialization();

  // TODO: initialize m_effectiveBulkModulus here. it might be needed by stress control before first updateStress.

  // GEOS_LOG_RANK_0( "Ez: " << m_defaultYoungModulusAxial << "\n" <<
  //                  "Ep: " << m_defaultYoungModulusTransverse << "\n" <<
  //                  "Nup: " << m_defaultPoissonRatioTransverse << "\n" <<
  //                  "Nuzp: " << m_defaultPoissonRatioAxialTransverse << "\n" <<
  //                  "Gzp: " << m_defaultShearModulusAxialTransverse << "\n" <<
  //                  "dEzdp: " << m_defaultYoungModulusAxialPressureDerivative << "\n" <<
  //                  "dEpdp: " << m_defaultYoungModulusTransversePressureDerivative << "\n" <<
  //                  "dGzpdp: " << m_defaultShearModulusAxialTransversePressureDerivative << "\n" <<
  //                  "sigmaFail: " << m_failureStrength << "\n" <<
  //                  "basalPlaneFractureEnergyReleaseRate: " << m_basalPlaneFractureEnergyReleaseRate << "\n" <<
  //                  "ds X2: " << m_distortionShearResponseX2 << "\n" <<
  //                  "ds Y1: " << m_distortionShearResponseY1 << "\n" <<
  //                  "ds Y2: " << m_distortionShearResponseY2 << "\n" <<
  //                  "ds M1: " << m_distortionShearResponseM1 << "\n" <<
  //                  "ips X2: " << m_inPlaneShearResponseX2 << "\n" <<
  //                  "ips Y1: " << m_inPlaneShearResponseY1 << "\n" <<
  //                  "ips Y2: " << m_inPlaneShearResponseY2 << "\n" <<
  //                  "ips M1: " << m_inPlaneShearResponseM1 << "\n" <<
  //                  "cs X2: " << m_coupledShearResponseX2 << "\n" <<
  //                  "cs Y1: " << m_coupledShearResponseY1 << "\n" <<
  //                  "cs Y2: " << m_coupledShearResponseY2 << "\n" <<
  //                  "cs M1: " << m_coupledShearResponseM1 << "\n" <<
  //                  "max ep: " << m_maximumPlasticStrain );

  // Add elastic constants check
  GEOS_THROW_IF( m_defaultYoungModulusAxial <= 0.0, "defaultYoungModulusAxial must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultYoungModulusTransverse <= 0.0, "defaultYoungModulusTransverse must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultShearModulusAxialTransverse <= 0.0, "defaultShearModulusAxialTransverse must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultYoungModulusAxialPressureDerivative < 0.0, "defaultYoungModulusAxialPressureDerivative must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultYoungModulusTransversePressureDerivative < 0.0, "defaultYoungModulusTransversePressureDerivative must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultShearModulusAxialTransversePressureDerivative < 0.0, "defaultShearModulusAxialTransversePressureDerivative must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultPoissonRatioAxialTransverse < -0.499999,  "defaultPoissonRatioAxialTransverse must be > -0.5 ", InputError );
  GEOS_THROW_IF( m_defaultPoissonRatioAxialTransverse > 0.499999,  "defaultPoissonRatioAxialTransverse must be < 0.5 ", InputError );
  GEOS_THROW_IF( m_defaultPoissonRatioTransverse < -0.499999,  "defaultPoissonRatioTransverse must be > -0.5 ", InputError );
  GEOS_THROW_IF( m_defaultPoissonRatioTransverse > 0.499999,  "defaultPoissonRatioTransverse must be < 0.5 ", InputError );

  real64 stability =
  1.0 - m_defaultPoissonRatioTransverse
      - 2.0 * m_defaultYoungModulusTransverse
              / m_defaultYoungModulusAxial
              * m_defaultPoissonRatioAxialTransverse
              * m_defaultPoissonRatioAxialTransverse;

  GEOS_THROW_IF( stability <= 1.e-12, "Transversely isotropic elastic constants are not positive definite.", InputError );

  GEOS_THROW_IF( m_failureStrength <= 0.0, "Maximum theoretical strength must be greater than 0", InputError );
  GEOS_THROW_IF( m_maximumPrincipalStressDamage != 0 &&
               m_maximumPrincipalStressDamage != 1, "Max Princ. Stress Damage flag should be 0 or 1", InputError );
  GEOS_THROW_IF( m_crackSpeed <= 0, "Crack speed should be positive", InputError );
  GEOS_THROW_IF( m_scaleFractureEnergyReleaseRate != 0 && m_scaleFractureEnergyReleaseRate != 1,  "Fracture energy scale flag should be 0 or 1", InputError );
  if( m_fractureEnergyStrengthScaleExponent < 0.0 )
  {
    m_fractureEnergyStrengthScaleExponent = ( m_scaleFractureEnergyReleaseRate == 1 ) ? 1.0 : 0.0;
  }
  GEOS_THROW_IF( m_fractureEnergyStrengthScaleExponent < 0.0, "fractureEnergyStrengthScaleExponent must be non-negative.", InputError );

  GEOS_THROW_IF( m_enableCrackTipStressConcentration != 0 &&  m_enableCrackTipStressConcentration != 1,  "Crack tip flag should be 0 or 1", InputError );

  GEOS_THROW_IF( m_basalPlaneFractureEnergyReleaseRate <= 0.0, "Basal plane fracture energy release rate must be a positive number.", InputError );
  GEOS_THROW_IF( m_totalFractureEnergyReleaseRate <= 0.0, "Total Fracture Energy Release Rate must be a positive number.", InputError );
  GEOS_THROW_IF( m_damageEvolutionExponent <= 0.0, "damageEvolutionExponent must be a positive number.", InputError );

  GEOS_THROW_IF( m_damagedMaterialFrictionalSlope < 0.0, "Damaged material frictional slope must be greater than 0", InputError );

  GEOS_THROW_IF( m_distortionShearResponseX2 <= 0.0, "Distortion shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_distortionShearResponseY1 < 0.0, "Distortion shear response y1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_distortionShearResponseY2 <= m_distortionShearResponseY1, "Distortion shear response y2 must > y1.", InputError );
  GEOS_THROW_IF( m_distortionShearResponseM1 < 0.0, "Distortion shear response m1 must be a positive number.", InputError );

  if( m_positiveDistortionShearResponseX2 < 0.0 )
  {
    m_positiveDistortionShearResponseX2 = m_distortionShearResponseX2;
  }
  if( m_positiveDistortionShearResponseY1 < 0.0 )
  {
    m_positiveDistortionShearResponseY1 = m_distortionShearResponseY1;
  }
  if( m_positiveDistortionShearResponseY2 < 0.0 )
  {
    m_positiveDistortionShearResponseY2 = m_distortionShearResponseY2;
  }
  if( m_positiveDistortionShearResponseM1 < 0.0 )
  {
    m_positiveDistortionShearResponseM1 = m_distortionShearResponseM1;
  }
  GEOS_THROW_IF( m_positiveDistortionShearResponseX2 <= 0.0, "Positive distortion shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_positiveDistortionShearResponseY1 < 0.0, "Positive distortion shear response y1 must be non-negative.", InputError );
  GEOS_THROW_IF( m_positiveDistortionShearResponseY2 <= m_positiveDistortionShearResponseY1, "Positive distortion shear response y2 must be greater than y1.", InputError );
  GEOS_THROW_IF( m_positiveDistortionShearResponseM1 < 0.0, "Positive distortion shear response m1 must be non-negative.", InputError );

  GEOS_THROW_IF( m_inPlaneShearResponseX2 <= 0.0, "In plane shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_inPlaneShearResponseY1 < 0.0, "In plane shear response y1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_inPlaneShearResponseY2 <= m_inPlaneShearResponseY1, "In plane shear response y2 must > y1.", InputError );
  GEOS_THROW_IF( m_inPlaneShearResponseM1 < 0.0, "In plane shear response m1 must be a positive number.", InputError );

  GEOS_THROW_IF( m_coupledShearResponseX2 <= 0.0, "Coupled shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_coupledShearResponseY1 < 0.0, "Coupled shear response y1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_coupledShearResponseY2 <= m_coupledShearResponseY1, "Coupled shear response y2 must > y1.", InputError );
  GEOS_THROW_IF( m_coupledShearResponseM1 < 0.0, "Coupled shear response m1 must be a positive number.", InputError );

  GEOS_THROW_IF( m_maximumPlasticStrain <= 0.0, "Maximum plastic strain must be a positive number.", InputError );
}


void Graphite::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Graphite, std::string const &, Group * const )
}
} /* namespace geos */
