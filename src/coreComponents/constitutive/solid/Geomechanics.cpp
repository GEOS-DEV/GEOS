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
 *  @file Geomechanics.cpp
 */

#include "Geomechanics.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{
Geomechanics::Geomechanics( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_b0( 0.0 ),
  m_b1( 0.0 ),
  m_b2( 0.0 ),
  m_b3( 0.0 ),
  m_b4( 0.0 ),
  m_dstrendh( 0.0 ),
  m_dfslopedh( 0.0 ),
  m_dpeakI1dh( 0.0 ),
  m_dcrdh( 0.0 ),
  m_g0( 0.0 ),
  m_g1( 0.0 ),
  m_g2( 0.0 ),
  m_g3( 0.0 ),
  m_g4( 0.0 ),
  m_p0( 0.0 ),
  m_p1( 0.0 ),
  m_p2( 0.0 ),
  m_p3( 0.0 ),
  m_p4( 0.0 ),
  m_peakI1( 0.0 ),
  m_fSlope( 0.0 ),
  m_fSlopeFailed( 0.0 ),
  m_stren( 0.0 ),
  m_ySlope( 0.0 ),
  m_beta( 1.0 ),
  m_t1RateDependence( 0.0 ),
  m_t2RateDependence( 0.0 ),
  m_fractureEnergyReleaseRate( 0.0 ),
  m_fractureSofteningExponent( 1.0 ),
  m_fractureStress( 0.0 ),
  m_initialTemperature( 0.0 ),
  m_Q( 0.0 ),
  m_brittleDuctileTransition( 0.0 ),
  m_damageEvolutionCriterion( 0 ),
  m_cr( 0.0 ),
  m_fluidBulkModulus(0.0 ),
  m_fluidInitialPressure( 0.0 ),
  m_enableBuckling( 0 ),
  m_bucklingLength( 1. ),
  m_bucklingAmplitude( 0. ),
  m_enableCreep( 0 ),
  m_creepC0( 0.0),
  m_creepC1( 0.0 ),
  m_creepC2( 0.0 ),
  m_creepA( 0.0 ),
  m_creepB( 0.0 ),
  m_creepC( 0.0 ),
  m_creepD( 1.0 ),
  m_creepE( 0.0 ),
  m_creepF( 1.0 ),
  m_creepG( 1.0 ),
  m_strainHardeningN( 0.0 ),
  m_strainHardeningK( 0.0 ),
  m_plasticStrainTolerance( 1.0e-10 ),
  m_stressReturnTolerance( 1.0e-6 ),
  m_maxAllowedSubcycles( 256 ),
  m_failedStepResponse( 2 ),
  m_bulkModulus(),
  m_shearModulus(),
  m_velocityGradient(),
  m_deformationGradient(),
  m_plasticStrain(),
  m_damage(),
  m_constitutiveUpdateFlag(),
  m_temperature(),
  m_lengthScale(),
  m_strengthScale()
{
  // register default values
  registerWrapper( viewKeyStruct::b0String(), &m_b0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 0" );

  registerWrapper( viewKeyStruct::b1String(), &m_b1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 1" );

  registerWrapper( viewKeyStruct::b2String(), &m_b2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 2" );

  registerWrapper( viewKeyStruct::b3String(), &m_b3 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 3" );

  registerWrapper( viewKeyStruct::b4String(), &m_b4 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 4" );

  registerWrapper( viewKeyStruct::dstrendhString(), &m_dstrendh ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Constant for STREN hardening rate" );

  registerWrapper( viewKeyStruct::dfslopedhString(), &m_dfslopedh ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Constant for FSLOPE hardening rate" );

  registerWrapper( viewKeyStruct::dpeakI1dhString(), &m_dpeakI1dh ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Constant for PEAKI1 hardening rate" );

  registerWrapper( viewKeyStruct::dcrdhString(), &m_dcrdh ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Constant for hardened CR" );

  registerWrapper( viewKeyStruct::g0String(), &m_g0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent shear shear modulus parameter 0" );

  registerWrapper( viewKeyStruct::g1String(), &m_g1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 1" );

  registerWrapper( viewKeyStruct::g2String(), &m_g2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 2" );

  registerWrapper( viewKeyStruct::g3String(), &m_g3 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 3" );

  registerWrapper( viewKeyStruct::g4String(), &m_g4 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 4" );

  registerWrapper( viewKeyStruct::p0String(), &m_p0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 0" );

  registerWrapper( viewKeyStruct::p1String(), &m_p1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 1" );

  registerWrapper( viewKeyStruct::p2String(), &m_p2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Crush curve parameter 2. Currently unused by this implementation." );

  registerWrapper( viewKeyStruct::p3String(), &m_p3 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 3" );

  registerWrapper( viewKeyStruct::p4String(), &m_p4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Crush curve parameter 4. Currently unused by this implementation." );

  registerWrapper( viewKeyStruct::crString(), &m_cr ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cap shape parameter" );

  registerWrapper( viewKeyStruct::fluidBulkModulusString(), &m_fluidBulkModulus ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid bulk modulus. Fluid effects are currently disabled; leave as zero." );

  registerWrapper( viewKeyStruct::fluidInitialPressureString(), &m_fluidInitialPressure ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid initial pressure. Fluid effects are currently disabled; leave as zero." );

  registerWrapper( viewKeyStruct::t1RateDependenceString(), &m_t1RateDependence ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rate dependence parameter 1. Currently unused by this implementation." );

  registerWrapper( viewKeyStruct::t2RateDependenceString(), &m_t2RateDependence ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rate dependence parameter 2. Currently unused by this implementation." );

  registerWrapper( viewKeyStruct::fractureEnergyReleaseRateString(), &m_fractureEnergyReleaseRate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fracture energy release rate parameter" );

  registerWrapper( viewKeyStruct::brittleDuctileTransitionString(), &m_brittleDuctileTransition ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "brittleDuctileTransition parameter" );

  registerWrapper( viewKeyStruct::fractureSofteningExponentString(), &m_fractureSofteningExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fracture softening exponent parameter" );

  registerWrapper( viewKeyStruct::fractureStressString(), &m_fractureStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fracture stress" );

  registerWrapper( viewKeyStruct::initialTemperatureString(), &m_initialTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "initial Temperature" );

  registerWrapper( viewKeyStruct::QString(), &m_Q ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "activation energy" );

  registerWrapper( viewKeyStruct::damageEvolutionCriterionString(), &m_damageEvolutionCriterion ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "damageEvolutionCriterion" );

  registerWrapper( viewKeyStruct::peakI1String(), &m_peakI1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Peak I1 shear limit parameter" );

  registerWrapper( viewKeyStruct::fSlopeString(), &m_fSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "F slope shear limit parameter" );

  registerWrapper( viewKeyStruct::fSlopeFailedString(), &m_fSlopeFailed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "F slope shear limit parameter after damage" );

  registerWrapper( viewKeyStruct::strenString(), &m_stren ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Stren shear limit parameter" );

  registerWrapper( viewKeyStruct::ySlopeString(), &m_ySlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Y slope shear limit parameter" );

  registerWrapper( viewKeyStruct::betaString(), &m_beta ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Nonassociativity parameter" );

  registerWrapper( viewKeyStruct::enableBucklingString(), &m_enableBuckling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Buckling flag" );

  registerWrapper( viewKeyStruct::bucklingLengthString(), &m_bucklingLength ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Buckling Length" );    

  registerWrapper( viewKeyStruct::bucklingAmplitudeString(), &m_bucklingAmplitude ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Buckling Amplitude" );

  registerWrapper( viewKeyStruct::creepString(), &m_enableCreep ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep flag" );

  registerWrapper( viewKeyStruct::creepC0String(), &m_creepC0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C0 parameter" );

  registerWrapper( viewKeyStruct::creepC1String(), &m_creepC1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C1 parameter" );

  registerWrapper( viewKeyStruct::creepC2String(), &m_creepC2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C2 parameter" );

  registerWrapper( viewKeyStruct::creepAString(), &m_creepA ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep A parameter" );
  
  registerWrapper( viewKeyStruct::creepBString(), &m_creepB ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep B parameter" );

  registerWrapper( viewKeyStruct::creepCString(), &m_creepC ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C parameter" );
  
  registerWrapper( viewKeyStruct::creepDString(), &m_creepD ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep D parameter" );

  registerWrapper( viewKeyStruct::creepEString(), &m_creepE ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep E parameter" );

  registerWrapper( viewKeyStruct::creepFString(), &m_creepF ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep F parameter" );

  registerWrapper( viewKeyStruct::creepGString(), &m_creepG ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep G parameter" );

  registerWrapper( viewKeyStruct::strainHardeningNString(), &m_strainHardeningN ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain Hardening n parameter" );

  registerWrapper( viewKeyStruct::strainHardeningKString(), &m_strainHardeningK ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain Hardening K parameter" );

  registerWrapper( viewKeyStruct::plasticStrainToleranceString(), &m_plasticStrainTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Tolerance in the volumetric plastic strain consistency bisection - units of strain" );

  registerWrapper( viewKeyStruct::stressReturnToleranceString(), &m_stressReturnTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Dimensionless precision of the non-hardening return relative to the trial stress increment" );

  registerWrapper( viewKeyStruct::maxAllowedSubcyclesString(), &m_maxAllowedSubcycles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Number of subcycles allowed before step is deemed to fail" );

  registerWrapper( viewKeyStruct::failedStepResponseString(), &m_failedStepResponse ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Treatment if requested subcycles exceed maxAllowedSubcycles: 0/1 cap and try, 2/3 mark failed. Failed convergence always sets constitutiveUpdateFlag=-1 for solver deletion." );

  // register fields
  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Bulk modulus" );
  
  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setInputFlag( InputFlags::FALSE).
    setDescription( "Shear modulus");

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient" );

  registerWrapper( viewKeyStruct::materialDirectionString(), &m_materialDirection ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Material direction" );

  registerWrapper( viewKeyStruct::deformationGradientString(), &m_deformationGradient ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Array of element/particle deformation gradient values" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Plastic strain" );

  registerWrapper( viewKeyStruct::porosityString(), &m_porosity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Porosity" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::constitutiveUpdateFlagString(), &m_constitutiveUpdateFlag ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Constitutive update status: 0=good, 1=subcycling warning, -1=failed update/delete flag" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 300.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point temperature values" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setApplyDefaultValue( DBL_MIN ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point length scale values" );

  registerWrapper( viewKeyStruct::strengthScaleString(), &m_strengthScale ).
    setApplyDefaultValue( DBL_MIN ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point strength scale values" );
}

Geomechanics::~Geomechanics()
{}

void Geomechanics::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  m_bulkModulus.resize( 0 );
  m_shearModulus.resize( 0 );
  m_velocityGradient.resize( 0, 3, 3 );
  m_materialDirection.resize( 0, 3 );
  m_deformationGradient.resize( 0, 3, 3);
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_porosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_constitutiveUpdateFlag.resize( 0, numConstitutivePointsPerParentIndex );
  //m_temperature.resize( 0, numConstitutivePointsPerParentIndex );
}

void Geomechanics::postInputInitialization()
{
  SolidBase::postInputInitialization();

  GEOS_THROW_IF( m_b0 <= 0.0, "b0 must be greater than 0", InputError );
  GEOS_THROW_IF( m_b0 + m_b1 <= 0.0, "b0 + b1 must be greater than 0", InputError );
  GEOS_THROW_IF( m_g0 <= 0.0, "g0 must be greater than 0", InputError );

  GEOS_THROW_IF( m_p0 >= 0.0, "p0 must be less than 0", InputError );
  GEOS_THROW_IF( m_p1 <= 0.0, "p1 must be greater than 0", InputError );
  GEOS_THROW_IF( m_p3 <= 0.0, "p3 must be greater than 0", InputError );

  GEOS_THROW_IF( m_fSlope < 0.0, "fSlope must be greater than or equal to 0", InputError );
  GEOS_THROW_IF( m_fSlopeFailed < 0.0, "fSlopeFailed must be greater than or equal to 0", InputError );
  GEOS_THROW_IF( m_fSlopeFailed > m_fSlope, "fSlopeFailed must be less than or equal to fSlope", InputError );
  GEOS_THROW_IF( m_peakI1 < 0.0, "peakI1 must be greater than or equal to 0", InputError );
  GEOS_THROW_IF( m_stren < 0.0, "stren must be greater than or equal to 0", InputError );
  GEOS_THROW_IF( m_ySlope < 0.0, "ySlope must be greater than or equal to 0", InputError );
  GEOS_THROW_IF( m_beta <= 0.0, "beta must be greater than 0", InputError );
  GEOS_THROW_IF( m_fractureSofteningExponent <= 0.0, "fracture softening exponent must be greater than 0", InputError );
  GEOS_THROW_IF( m_cr <= 0.0 || m_cr >= 1.0, "cr must satisfy 0 < cr < 1", InputError );

  auto const isNearlyZero = []( real64 const value ) -> bool { return std::abs( value ) < 1.0e-12; };
  real64 const ySlopeInitial = std::min( 0.99999 * m_fSlope, m_ySlope );
  bool const validLimitSurface =
    ( m_fSlope > 0.0 && m_peakI1 >= 0.0 && isNearlyZero( m_stren ) && isNearlyZero( ySlopeInitial ) ) ||
    ( isNearlyZero( m_fSlope ) && isNearlyZero( m_peakI1 ) && m_stren > 0.0 && isNearlyZero( ySlopeInitial ) ) ||
    ( m_fSlope > 0.0 && isNearlyZero( ySlopeInitial ) && m_stren > 0.0 && isNearlyZero( m_peakI1 ) ) ||
    ( m_fSlope > ySlopeInitial && ySlopeInitial > 0.0 && m_stren > ySlopeInitial * m_peakI1 && m_peakI1 >= 0.0 );
  GEOS_THROW_IF( !validLimitSurface,
                 "Invalid shear limit surface parameters. Expected one of: linear Drucker-Prager, Von Mises, zero-PEAKI1 transition, or nonlinear Drucker-Prager.",
                 InputError );

  GEOS_THROW_IF( !isNearlyZero( m_fluidBulkModulus ) || !isNearlyZero( m_fluidInitialPressure ),
                 "Fluid effects are currently disabled in Geomechanics. Set fluidBulkModulus and fluidInitialPressure to zero.",
                 InputError );

  GEOS_THROW_IF( m_enableBuckling < 0 || m_enableBuckling > 2,
                 "enableBuckling must be 0, 1, or 2", InputError );
  if( m_enableBuckling > 0 )
  {
    GEOS_THROW_IF( m_bucklingLength <= 0.0, "bucklingLength must be greater than 0 when buckling is enabled", InputError );
    GEOS_THROW_IF( m_bucklingAmplitude < 0.0 || m_bucklingAmplitude > 1.0,
                   "bucklingAmplitude must satisfy 0 <= bucklingAmplitude <= 1", InputError );
  }

  GEOS_THROW_IF( m_enableCreep < 0 || m_enableCreep > 1, "enableCreep must be 0 or 1", InputError );
  if( m_enableCreep == 1 )
  {
    GEOS_THROW_IF( m_creepB <= 0.0, "creepB must be greater than 0 when creep is enabled", InputError );
    GEOS_THROW_IF( m_creepG <= 0.0, "creepG must be greater than 0 when creep is enabled", InputError );
    GEOS_THROW_IF( !isNearlyZero( m_Q ) && m_initialTemperature <= 0.0,
                   "initialTemperature must be greater than 0 when nonzero creep activation energy Q is used",
                   InputError );
  }

  GEOS_THROW_IF( m_plasticStrainTolerance <= 0.0, "plasticStrainTolerance must be greater than 0", InputError );
  GEOS_THROW_IF( m_stressReturnTolerance <= 0.0, "stressReturnTolerance must be greater than 0", InputError );
  GEOS_THROW_IF( m_maxAllowedSubcycles <= 0, "maxAllowedSubcycles must be greater than 0", InputError );
  GEOS_THROW_IF( m_failedStepResponse < 0 || m_failedStepResponse > 3,
                 "failedStepResponse must be 0, 1, 2, or 3", InputError );
}



void Geomechanics::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Geomechanics, std::string const &, Group * const )
}
} /* namespace geos */
