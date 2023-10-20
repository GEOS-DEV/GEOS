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

using namespace dataRepository;

namespace constitutive
{

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid( string const & name, Group * const parent )
  : MultiFluidBase( name, parent )
{
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

integer CompositionalMultiphaseFluid::getWaterPhaseIndex() const
{
  return m_aqueousPhaseIndex;
}

void CompositionalMultiphaseFluid::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

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

  m_phaseTypes.resize( NP );
  std::transform( m_phaseNames.begin(), m_phaseNames.end(), m_phaseTypes.begin(), [this]( string const & name ){ return this->getPhaseType( name ); } );
  m_eosTypes.resize( NP );
  std::transform( m_equationsOfState.begin(), m_equationsOfState.end(), m_eosTypes.begin(),
    []( string const & name ){ return EnumStrings<EquationOfStateType>::fromString(name); } );

  // Reorder phases: liquid, vapour, aqueous
  std::multimap< PhaseType, EquationOfStateType > ordering;
  for( integer ip = 0; ip < NP; ++ip )
  {
    ordering.insert( {m_phaseTypes[ip], m_eosTypes[ip]} );
  }
  integer ip = 0;
  for( const auto [phase, eos] : ordering )
  {
    m_phaseTypes[ip] = phase;
    m_eosTypes[ip] = eos;
    if( phase == PhaseType::aqueous )
    {
      m_aqueousPhaseIndex = ip;
    }
    ++ip;
  }
}

void CompositionalMultiphaseFluid::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createFluid();
}

void CompositionalMultiphaseFluid::createFluid()
{
  m_fluid = std::make_unique< IFluid >();
}

std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid::deliverClone( string const & name,
                                            Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  CompositionalMultiphaseFluid & fluid = dynamicCast< CompositionalMultiphaseFluid & >( *clone );
  fluid.m_phaseTypes = m_phaseTypes;
  fluid.m_eosTypes = m_eosTypes;
  fluid.createFluid();
  return clone;
}

CompositionalMultiphaseFluid::KernelWrapper::
  KernelWrapper( CompositionalMultiphaseFluid::IFluid & fluid,
                 arrayView1d< PhaseType > const & phaseTypes,
                 arrayView1d< geos::real64 const > const & componentMolarWeight,
                 bool useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( componentMolarWeight,
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_fluid( fluid ),
  m_phaseTypes( phaseTypes )
{}

CompositionalMultiphaseFluid::KernelWrapper
CompositionalMultiphaseFluid::createKernelWrapper()
{
  return KernelWrapper( *m_fluid,
                        m_phaseTypes,
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

PhaseType
CompositionalMultiphaseFluid::getPhaseType( string const & name ) const
{
  map< string, PhaseType > const phaseTypes =
  {
    { "oil", PhaseType::liquid },
    { "liquid", PhaseType::liquid },
    { "gas", PhaseType::vapour },
    { "vapour", PhaseType::vapour },
    { "wat", PhaseType::aqueous },
    { "water", PhaseType::aqueous },
    { "aqueous", PhaseType::aqueous }
  };
  return findOption( phaseTypes, name, viewKeyStruct::phaseNamesString(), getFullName() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, string const &, Group * const )

} // namespace constitutive

} // namespace geos
