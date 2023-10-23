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

#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace constitutive
{

template< typename FLASH, typename ... PHASES >
CompositionalMultiphaseFluid< FLASH, PHASES... >::
CompositionalMultiphaseFluid( string const & name, Group * const parent )
  : MultiFluidBase( name, parent )
{
  using InputFlags = dataRepository::InputFlags;

  getWrapperBase( viewKeyStruct::componentNamesString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::equationsOfStateString(), &m_equationsOfState ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of equation of state types for each phase. Available options are:\n"
                    + EnumStrings< EquationOfStateType >::concat( "\n *" ) );

  registerWrapper( viewKeyStruct::componentCriticalPressureString(), &m_componentCriticalPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical pressures" );

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString(), &m_componentCriticalTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical temperatures" );

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

template< typename FLASH, typename ... PHASES >
integer CompositionalMultiphaseFluid< FLASH, PHASES... >::getWaterPhaseIndex() const
{
  return -1;
}

template< typename FLASH, typename ... PHASES >
void CompositionalMultiphaseFluid< FLASH, PHASES... >::postProcessInput()
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
  checkInputSize( m_equationsOfState, NP, viewKeyStruct::equationsOfStateString() );
  checkInputSize( m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString() );
  checkInputSize( m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString() );
  checkInputSize( m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString() );
  checkInputSize( m_equationsOfState, NP, viewKeyStruct::equationsOfStateString() );

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
}

template< typename FLASH, typename ... PHASES >
void CompositionalMultiphaseFluid< FLASH, PHASES... >::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();

  // Create the fluid models
  createModels();
}

template< typename FLASH, typename ... PHASES >
std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid< FLASH, PHASES... >::deliverClone( string const & name,
                                                                Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  return clone;
}

template< typename FLASH, typename ... PHASES >
typename CompositionalMultiphaseFluid< FLASH, PHASES... >::KernelWrapper
CompositionalMultiphaseFluid< FLASH, PHASES... >::createKernelWrapper()
{
  return KernelWrapper( *m_componentProperties,
                        *m_flash,
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
template< typename FLASH, typename ... PHASES >
void CompositionalMultiphaseFluid< FLASH, PHASES... >::createModels()
{
  m_componentProperties = std::make_unique< compositional::ComponentProperties >(
    m_componentCriticalPressure,
    m_componentCriticalTemperature,
    m_componentAcentricFactor,
    m_componentVolumeShift,
    m_componentBinaryCoeff );

  m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                       m_componentNames,
                                       m_componentMolarWeight,
                                       *m_componentProperties );

  std::apply( [&]( std::unique_ptr< PHASES > & ... phaseModel ) {
    integer phaseIndex = 0;
    (void(phaseModel = std::make_unique< PHASES >(
            GEOS_FMT( "{}_PhaseModel{}", getName(), phaseIndex++ ),
            m_componentNames,
            m_componentMolarWeight,
            *m_componentProperties )), ...);
  }, m_phases );
}

// Explicit instantiation of the model template.
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashPRPR,
    compositional::PhaseModel< compositional::CubicEOSDensityPR, compositional::ConstantViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CubicEOSDensityPR, compositional::ConstantViscosity, compositional::NullModel > >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalTwoPhaseFluidPengRobinson, string const &, dataRepository::Group * const )

} // namespace constitutive

} // namespace geos
