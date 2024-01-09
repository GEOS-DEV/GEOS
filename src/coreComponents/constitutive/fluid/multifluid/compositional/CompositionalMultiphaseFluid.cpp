/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseFluid.cpp
 */

#include "CompositionalMultiphaseFluid.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace constitutive
{

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::
CompositionalMultiphaseFluid( string const & name, Group * const parent )
  : MultiFluidBase( name, parent )
{
  using InputFlags = dataRepository::InputFlags;

  getWrapperBase( viewKeyStruct::componentNamesString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::componentCriticalPressureString(), &m_componentCriticalPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical pressures" );

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString(), &m_componentCriticalTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical temperatures" );

  registerWrapper( viewKeyStruct::componentCriticalVolumeString(), &m_componentCriticalVolume ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component critical volumnes" );

  registerWrapper( viewKeyStruct::componentAcentricFactorString(), &m_componentAcentricFactor ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component acentric factors" );

  registerWrapper( viewKeyStruct::componentVolumeShiftString(), &m_componentVolumeShift ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component volume shifts" );

  registerWrapper( viewKeyStruct::componentBinaryCoeffString(), &m_componentBinaryCoeff ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Table of binary interaction coefficients" );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
integer CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

// Naming conventions
namespace compositional
{
template< int NP > struct PhaseName {};
template<> struct PhaseName< 2 > { static constexpr char const * catalogName() { return "TwoPhase"; } };
template<> struct PhaseName< 3 > { static constexpr char const * catalogName() { return "ThreePhase"; } };
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
string CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::catalogName()
{
  return GEOS_FMT( "Compositional{}Fluid{}{}",
                   compositional::PhaseName< FLASH::KernelWrapper::getNumberOfPhases() >::catalogName(),
                   FLASH::catalogName(),
                   PHASE1::Viscosity::catalogName() );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

  GEOS_THROW_IF_NE_MSG( NP, NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases in '{}'. There should be {} phases",
                                  getFullName(), viewKeyStruct::phaseNamesString(), NUM_PHASES ),
                        InputError );

  auto const checkInputSize = [&]( auto const & array, integer const expected, string const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), expected,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );

  };
  checkInputSize( m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString() );
  checkInputSize( m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString() );
  checkInputSize( m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString() );

  if( m_componentCriticalVolume.empty() )
  {
    m_componentCriticalVolume.resize( NC );
    calculateCriticalVolume( m_componentCriticalPressure,
                             m_componentCriticalTemperature,
                             m_componentCriticalVolume );
  }
  checkInputSize( m_componentCriticalVolume, NC, viewKeyStruct::componentCriticalVolumeString() );

  if( m_componentVolumeShift.empty() )
  {
    m_componentVolumeShift.resize( NC );
    m_componentVolumeShift.zero();
  }
  checkInputSize( m_componentVolumeShift, NC, viewKeyStruct::componentVolumeShiftString() );

  if( m_componentBinaryCoeff.empty() )
  {
    m_componentBinaryCoeff.resize( NC, NC );
    m_componentBinaryCoeff.zero();
  }
  checkInputSize( m_componentBinaryCoeff, NC * NC, viewKeyStruct::componentBinaryCoeffString() );

  // Binary interaction coefficients should be symmetric and have zero diagonal
  GEOS_THROW_IF_NE_MSG( m_componentBinaryCoeff.size( 0 ), NC,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentBinaryCoeffString() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_componentBinaryCoeff.size( 1 ), NC,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentBinaryCoeffString() ),
                        InputError );
  for( integer ic = 0; ic < NC; ++ic )
  {
    GEOS_THROW_IF_GT_MSG( LvArray::math::abs( m_componentBinaryCoeff( ic, ic )), MultiFluidConstants::epsilon,
                          GEOS_FMT( "{}: {} entry at ({},{}) is {}: should be zero", getFullName(), viewKeyStruct::componentBinaryCoeffString(),
                                    ic, ic, m_componentBinaryCoeff( ic, ic ) ),
                          InputError );
    for( integer jc = ic + 1; jc < NC; ++jc )
    {
      real64 const difference = LvArray::math::abs( m_componentBinaryCoeff( ic, jc )-m_componentBinaryCoeff( jc, ic ));
      GEOS_THROW_IF_GT_MSG( difference, MultiFluidConstants::epsilon,
                            GEOS_FMT( "{}: {} entry at ({},{}) is {} and is different from entry at ({},{}) which is {}",
                                      getFullName(), viewKeyStruct::componentBinaryCoeffString(),
                                      ic, jc, m_componentBinaryCoeff( ic, jc ), jc, ic, m_componentBinaryCoeff( jc, ic ) ),
                            InputError );
    }
  }
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();

  // Create the fluid models
  createModels();
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::deliverClone( string const & name,
                                                                             Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  return clone;
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
typename CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::KernelWrapper
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createKernelWrapper()
{
  return KernelWrapper( *m_componentProperties,
                        *m_flash,
                        *m_phase1,
                        *m_phase2,
                        *m_phase3,
                        m_componentMolarWeight,
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

// Create the fluid models
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createModels()
{
  m_componentProperties = std::make_unique< compositional::ComponentProperties >(
    m_componentNames,
    m_componentMolarWeight,
    m_componentCriticalPressure,
    m_componentCriticalTemperature,
    m_componentCriticalVolume,
    m_componentAcentricFactor,
    m_componentVolumeShift,
    m_componentBinaryCoeff );

  m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                       *m_componentProperties );

  m_phase1 = std::make_unique< PHASE1 >( GEOS_FMT( "{}_PhaseModel1", getName() ),
                                         *m_componentProperties );

  m_phase2 = std::make_unique< PHASE2 >( GEOS_FMT( "{}_PhaseModel2", getName() ),
                                         *m_componentProperties );

  m_phase3 = std::make_unique< PHASE3 >( GEOS_FMT( "{}_PhaseModel3", getName() ),
                                         *m_componentProperties );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::calculateCriticalVolume(
  arrayView1d< const real64 > const criticalPressure,
  arrayView1d< const real64 > const criticalTemperature,
  arrayView1d< real64 > const criticalVolume ) const
{
  integer const numComponents = criticalPressure.size( 0 );
  for( integer ic=0; ic<numComponents; ++ic )
  {
    criticalVolume[ic] = 2.215e-6 * criticalTemperature[ic] / (0.025 + 1e-6*criticalPressure[ic] );   // m^3/mol
  }
}

// Explicit instantiation of the model template.
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashPRPR,
    compositional::PhaseModel< compositional::CompositionalDensity< compositional::CubicEOSPR >, compositional::ConstantViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity< compositional::CubicEOSPR >, compositional::ConstantViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashSRKSRK,
    compositional::PhaseModel< compositional::CompositionalDensity< compositional::CubicEOSSRK >, compositional::ConstantViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity< compositional::CubicEOSSRK >, compositional::ConstantViscosity, compositional::NullModel > >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhasePengRobinsonConstantViscosity,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity,
                        string const &,
                        dataRepository::Group * const )

} // namespace constitutive

} // namespace geos
