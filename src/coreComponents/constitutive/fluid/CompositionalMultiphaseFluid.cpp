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
  : MultiFluidPVTPackageWrapper( name, parent )
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

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid() = default;

std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid::deliverClone( string const & name,
                                            Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidPVTPackageWrapper::deliverClone( name, parent );
  CompositionalMultiphaseFluid & fluid = dynamicCast< CompositionalMultiphaseFluid & >( *clone );

  fluid.createFluid();
  return clone;
}

namespace
{

template< typename ARRAY >
void checkInputSize( ARRAY const & array, localIndex const expected, string const & attr )
{
  GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                         "CompositionalMultiphaseFluid: invalid number of entries in " << attr << " attribute",
                         InputError );

}

}

void CompositionalMultiphaseFluid::postProcessInput()
{
  MultiFluidPVTPackageWrapper::postProcessInput();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

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

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
