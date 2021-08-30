/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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

#include "pvt/pvt.hpp"

#include <map>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

namespace
{

pvt::EOS_TYPE getCompositionalEosType( string const & name )
{
  static std::map< string, pvt::EOS_TYPE > const eosTypes =
  {
    { "PR", pvt::EOS_TYPE::PENG_ROBINSON },
    { "SRK", pvt::EOS_TYPE::REDLICH_KWONG_SOAVE }
  };
  auto const it = eosTypes.find( name );
  GEOSX_ERROR_IF( it == eosTypes.end(), "Compositional EOS type not supported by PVTPackage: " << name );
  return it->second;
}

} // namespace

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid( string const & name, Group * const parent )
  : MultiFluidBase( name, parent )
{
  getWrapperBase( viewKeyStruct::componentNamesString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::equationsOfStateString(), &m_equationsOfState ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of equation of state types for each phase" );

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

namespace
{

template< typename ARRAY >
void checkInputSize( ARRAY const & array, integer const expected, string const & attr, string const & name )
{
  GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                         name << ": invalid number of entries in " << attr << " attribute",
                         InputError );

}

pvt::PHASE_TYPE getPVTPackagePhaseType( string const & name, string const & groupName )
{
  static map< string, pvt::PHASE_TYPE > const phaseTypes
  {
    { "gas", pvt::PHASE_TYPE::GAS },
    { "oil", pvt::PHASE_TYPE::OIL },
    { "water", pvt::PHASE_TYPE::LIQUID_WATER_RICH }
  };
  return findOption( phaseTypes, name, MultiFluidBase::viewKeyStruct::phaseNamesString(), groupName );
}

}

void CompositionalMultiphaseFluid::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  m_phaseTypes.resize( numFluidPhases() );
  std::transform( m_phaseNames.begin(), m_phaseNames.end(), m_phaseTypes.begin(),
                  [&]( string const & name ){ return getPVTPackagePhaseType( name, getFullName() ); } );

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

  checkInputSize( m_equationsOfState, NP, viewKeyStruct::equationsOfStateString(), getFullName() );
  checkInputSize( m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString(), getFullName() );
  checkInputSize( m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString(), getFullName() );
  checkInputSize( m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString(), getFullName() );
  checkInputSize( m_equationsOfState, NP, viewKeyStruct::equationsOfStateString(), getFullName() );

  if( m_componentVolumeShift.empty() )
  {
    m_componentVolumeShift.resize( NC );
    m_componentVolumeShift.zero();
  }
  checkInputSize( m_componentVolumeShift, NC, viewKeyStruct::componentVolumeShiftString(), getFullName() );

  if( m_componentBinaryCoeff.empty() )
  {
    m_componentBinaryCoeff.resize( NC, NC );
    m_componentBinaryCoeff.zero();
  }
  checkInputSize( m_componentBinaryCoeff, NC * NC, viewKeyStruct::componentBinaryCoeffString(), getFullName() );
}

void CompositionalMultiphaseFluid::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createFluid();
}

void CompositionalMultiphaseFluid::createFluid()
{
  std::vector< pvt::EOS_TYPE > eos( numFluidPhases() );
  std::transform( m_equationsOfState.begin(), m_equationsOfState.end(), eos.begin(), getCompositionalEosType );

  std::vector< pvt::PHASE_TYPE > phases( m_phaseTypes.begin(), m_phaseTypes.end() );
  std::vector< string > const components( m_componentNames.begin(), m_componentNames.end() );
  std::vector< double > const Mw( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );
  std::vector< double > const Tc( m_componentCriticalTemperature.begin(), m_componentCriticalTemperature.end() );
  std::vector< double > const Pc( m_componentCriticalPressure.begin(), m_componentCriticalPressure.end() );
  std::vector< double > const Omega( m_componentAcentricFactor.begin(), m_componentAcentricFactor.end() );

  m_fluid = pvt::MultiphaseSystemBuilder::buildCompositional( pvt::COMPOSITIONAL_FLASH_TYPE::NEGATIVE_OIL_GAS, phases, eos,
                                                              components, Mw, Tc, Pc, Omega );
}

std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid::deliverClone( string const & name,
                                            Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  CompositionalMultiphaseFluid & fluid = dynamicCast< CompositionalMultiphaseFluid & >( *clone );
  fluid.m_phaseTypes = m_phaseTypes;
  fluid.createFluid();
  return clone;
}

CompositionalMultiphaseFluid::KernelWrapper::
  KernelWrapper( pvt::MultiphaseSystem & fluid,
                 arrayView1d< pvt::PHASE_TYPE > const & phaseTypes,
                 arrayView1d< geosx::real64 const > const & componentMolarWeight,
                 bool useMass,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseFraction,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseDensity,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseMassDensity,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseViscosity,
                 MultiFluidBase::KernelWrapper::PhaseCompViews const & phaseCompFraction,
                 MultiFluidBase::KernelWrapper::FluidPropViews const & totalDensity )
  : MultiFluidBase::KernelWrapper( componentMolarWeight,
                                   useMass,
                                   phaseFraction,
                                   phaseDensity,
                                   phaseMassDensity,
                                   phaseViscosity,
                                   phaseCompFraction,
                                   totalDensity ),
  m_fluid( fluid ),
  m_phaseTypes( phaseTypes )
{}

CompositionalMultiphaseFluid::KernelWrapper
CompositionalMultiphaseFluid::createKernelWrapper() const
{
  return KernelWrapper( *m_fluid,
                        m_phaseTypes,
                        m_componentMolarWeight,
                        m_useMass,
                        { m_phaseFraction,
                          m_dPhaseFraction_dPressure,
                          m_dPhaseFraction_dTemperature,
                          m_dPhaseFraction_dGlobalCompFraction },
                        { m_phaseDensity,
                          m_dPhaseDensity_dPressure,
                          m_dPhaseDensity_dTemperature,
                          m_dPhaseDensity_dGlobalCompFraction },
                        { m_phaseMassDensity,
                          m_dPhaseMassDensity_dPressure,
                          m_dPhaseMassDensity_dTemperature,
                          m_dPhaseMassDensity_dGlobalCompFraction },
                        { m_phaseViscosity,
                          m_dPhaseViscosity_dPressure,
                          m_dPhaseViscosity_dTemperature,
                          m_dPhaseViscosity_dGlobalCompFraction },
                        { m_phaseCompFraction,
                          m_dPhaseCompFraction_dPressure,
                          m_dPhaseCompFraction_dTemperature,
                          m_dPhaseCompFraction_dGlobalCompFraction },
                        { m_totalDensity,
                          m_dTotalDensity_dPressure,
                          m_dTotalDensity_dTemperature,
                          m_dTotalDensity_dGlobalCompFraction } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
