/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalEnthalpy.cpp
 */

#include "CompositionalEnthalpy.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/HeatCapacityCoefficients.hpp"
#include "common/PhysicsConstants.hpp"

namespace geos
{
namespace constitutive
{
namespace compositional
{

CompositionalEnthalpyUpdate::CompositionalEnthalpyUpdate( integer const phaseIndex,
                                                          integer const vapourIndex,
                                                          PhaseType const phaseType,
                                                          EquationOfStateType const equationOfState,
                                                          real64 const refTemperature,
                                                          arrayView2d< real64 const > const & referenceEnthalpy,
                                                          arrayView2d< real64 const > const & coefficients ):
  m_phaseIndex( phaseIndex ),
  m_vapourIndex( vapourIndex ),
  m_phaseType( phaseType ),
  m_equationOfState( equationOfState ),
  m_refTemperature( refTemperature ),
  m_referenceEnthalpy( referenceEnthalpy ),
  m_coefficients( coefficients )
{}

CompositionalEnthalpy::CompositionalEnthalpy( string const & name,
                                              ComponentProperties const & componentProperties,
                                              integer const phaseIndex,
                                              ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties ),
  m_phaseIndex( phaseIndex ),
  m_heatCapacityCoefficients ( modelParameters.get< HeatCapacityCoefficients >() )
{
  EquationOfState const * equationOfState = modelParameters.get< EquationOfState >();
  string const eosName = equationOfState->m_equationsOfStateNames[phaseIndex];
  m_equationOfState = EnumStrings< EquationOfStateType >::fromString( eosName );

  integer const numPhases = m_heatCapacityCoefficients->m_phaseTypes.size();

  m_phaseType = static_cast< PhaseType >(m_heatCapacityCoefficients->m_phaseTypes[phaseIndex]);
  for( integer ip = 0; ip < numPhases; ip++ )
  {
    if( isPhaseType( m_heatCapacityCoefficients->m_phaseTypes[ip], PhaseType::VAPOUR ) )
    {
      m_vapourIndex = ip;
    }
  }
}

CompositionalEnthalpy::KernelWrapper
CompositionalEnthalpy::createKernelWrapper() const
{
  return KernelWrapper( m_phaseIndex,
                        m_vapourIndex,
                        m_phaseType,
                        m_equationOfState,
                        m_heatCapacityCoefficients->m_referenceTemperature,
                        m_heatCapacityCoefficients->m_referenceEnthalpy.toViewConst(),
                        m_heatCapacityCoefficients->m_coefficients.toViewConst() );
}

std::unique_ptr< ModelParameters >
CompositionalEnthalpy::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  auto params = EquationOfState::create( std::move( parameters ) );
  params = HeatCapacityCoefficients::create( std::move( params ) );
  return params;
}

} // namespace compositional
} // namespace constitutive
} // namespace geos
