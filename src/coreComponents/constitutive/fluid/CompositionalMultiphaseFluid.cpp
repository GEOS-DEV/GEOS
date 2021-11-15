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
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"

#include "pvt/pvt.hpp"

#include <map>
#include <utility>

namespace geosx
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

integer CompositionalMultiphaseFluid::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

void CompositionalMultiphaseFluid::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  auto const getPVTPackagePhaseType = [&]( string const & phaseName )
  {
    static map< string, pvt::PHASE_TYPE > const phaseTypes
    {
      { "gas", pvt::PHASE_TYPE::GAS },
      { "oil", pvt::PHASE_TYPE::OIL },
      { "water", pvt::PHASE_TYPE::LIQUID_WATER_RICH }
    };
    return findOption( phaseTypes, phaseName, viewKeyStruct::phaseNamesString(), getFullName() );
  };

  m_phaseTypes.resize( numFluidPhases() );
  std::transform( m_phaseNames.begin(), m_phaseNames.end(), m_phaseTypes.begin(), getPVTPackagePhaseType );

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

  auto const checkInputSize = [&]( auto const & array, integer const expected, string const & attribute )
  {
    GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                           GEOSX_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
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

void CompositionalMultiphaseFluid::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createFluid();
}

void CompositionalMultiphaseFluid::createFluid()
{
  auto const getCompositionalEosType = [&]( string const & name )
  {
    static map< string, pvt::EOS_TYPE > const eosTypes =
    {
      { "PR", pvt::EOS_TYPE::PENG_ROBINSON },
      { "SRK", pvt::EOS_TYPE::REDLICH_KWONG_SOAVE }
    };
    return findOption( eosTypes, name, viewKeyStruct::phaseNamesString(), getFullName() );
  };

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
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( componentMolarWeight,
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
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
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
