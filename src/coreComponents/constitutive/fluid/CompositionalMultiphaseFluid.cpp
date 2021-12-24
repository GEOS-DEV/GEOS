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

template< typename VISCOSITY_MODEL >
CompositionalMultiphaseFluid< VISCOSITY_MODEL >::CompositionalMultiphaseFluid( string const & name, Group * const parent )
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

  registerWrapper( viewKeyStruct::componentCriticalVolumeString(), &m_componentCriticalVolume ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component critical volumes" );

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

template< typename VISCOSITY_MODEL >
integer CompositionalMultiphaseFluid< VISCOSITY_MODEL >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}


template< typename VISCOSITY_MODEL >
void CompositionalMultiphaseFluid< VISCOSITY_MODEL >::postProcessInput()
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

  // if unspecified, estimate critical volume using Ihmels' (2010) correlation 
  // reference: http://dx.doi.org/10.1021/je100167w

  if( m_componentCriticalVolume.empty() )
  {
    m_componentCriticalVolume.resize( NC );
    for(localIndex c=0; c<NC; ++c)
    {
      m_componentCriticalVolume[c] = 2.215e-6 * m_componentCriticalTemperature[c] / (0.025 + 1e-6*m_componentCriticalPressure[c] ); // m^3/mol
    }
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
}

template< typename VISCOSITY_MODEL >
void CompositionalMultiphaseFluid< VISCOSITY_MODEL >::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createFluid();
}

template< typename VISCOSITY_MODEL >
void CompositionalMultiphaseFluid< VISCOSITY_MODEL >::createFluid()
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

  m_viscosity = std::make_unique< VISCOSITY_MODEL >( getName() + '_' + VISCOSITY_MODEL::catalogName(), 
                                                     m_componentNames, 
                                                     m_componentMolarWeight,
                                                     m_componentCriticalPressure,
                                                     m_componentCriticalTemperature,
                                                     m_componentCriticalVolume,
                                                     m_componentAcentricFactor,
                                                     m_componentVolumeShift,
                                                     m_componentBinaryCoeff );
}

template< typename VISCOSITY_MODEL >
std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid< VISCOSITY_MODEL >::deliverClone( string const & name,
                                                               Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  CompositionalMultiphaseFluid & fluid = dynamicCast< CompositionalMultiphaseFluid & >( *clone );
  fluid.m_phaseTypes = m_phaseTypes;
  fluid.createFluid();
  return clone;
}

template< typename VISCOSITY_MODEL >
CompositionalMultiphaseFluid< VISCOSITY_MODEL >::KernelWrapper::
  KernelWrapper( pvt::MultiphaseSystem & fluid,
                 VISCOSITY_MODEL const & viscosity,
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
  m_viscosity( viscosity.createKernelWrapper() ),
  m_phaseTypes( phaseTypes )
{}

template< typename VISCOSITY_MODEL >
typename CompositionalMultiphaseFluid< VISCOSITY_MODEL >::KernelWrapper
CompositionalMultiphaseFluid< VISCOSITY_MODEL >::createKernelWrapper()
{
  return KernelWrapper( *m_fluid,
                        *m_viscosity,
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

template class CompositionalMultiphaseFluid< PVTProps::LohrenzBrayClarkViscosity >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluidLBC, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
