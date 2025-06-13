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

CompositionalEnthalpyUpdate::CompositionalEnthalpyUpdate( EquationOfStateType const equationOfState,
                                                          arrayView1d< real64 const > const & referenceEnthalpy,
                                                          arrayView2d< real64 const > const & coefficients )
  : m_equationOfState( equationOfState ),
  m_referenceEnthalpy( referenceEnthalpy ),
  m_coefficients( coefficients )
{}

CompositionalEnthalpy::CompositionalEnthalpy( string const & name,
                                              ComponentProperties const & componentProperties,
                                              integer const phaseIndex,
                                              ModelParameters const & modelParameters )
  : FunctionBase( name, componentProperties ),
  m_heatCapacityCoefficients ( modelParameters.get< HeatCapacityCoefficients >() )
{
  EquationOfState const * equationOfState = modelParameters.get< EquationOfState >();
  string const eosName = equationOfState->m_equationsOfStateNames[phaseIndex];
  m_equationOfState = EnumStrings< EquationOfStateType >::fromString( eosName );

  const auto & heatCapacityCoefficients = m_heatCapacityCoefficients->m_coefficients;
  integer const numComps = heatCapacityCoefficients.size( 1 );
  integer const numCoeffs = heatCapacityCoefficients.size( 2 );
  m_referenceEnthalpy.resize( numComps );
  m_coefficients.resize( numComps, numCoeffs );
  constexpr real64 R = constants::gasConstant;
  for( integer ic = 0; ic < numComps; ic++ )
  {
    for( integer j = 0; j < numCoeffs; j++ )
    {
      m_coefficients( ic, j ) = R * heatCapacityCoefficients( phaseIndex, ic, j );
    }
    // Calculate the enthalpy at the reference temperature
    real64 refEnthalpy = 0.0;
    real64 refHeatCapacity = 0.0;
    KernelWrapper::evaluatePolynomial( m_heatCapacityCoefficients->m_referenceTemperature[ic],
                                       m_coefficients[ic],
                                       refEnthalpy,
                                       refHeatCapacity );
    m_referenceEnthalpy[ic] = m_heatCapacityCoefficients->m_referenceEnthalpy( phaseIndex, ic ) - refEnthalpy;
  }
}

CompositionalEnthalpy::KernelWrapper
CompositionalEnthalpy::createKernelWrapper() const
{
  return KernelWrapper( m_equationOfState,
                        m_referenceEnthalpy.toViewConst(),
                        m_coefficients.toViewConst());
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
