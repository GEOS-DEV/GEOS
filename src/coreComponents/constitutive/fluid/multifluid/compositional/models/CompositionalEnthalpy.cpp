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

namespace geos
{
namespace constitutive
{
namespace compositional
{

CompositionalEnthalpyUpdate::CompositionalEnthalpyUpdate( EquationOfStateType const equationOfState )
  : m_equationOfState( equationOfState )
{}

CompositionalEnthalpy::CompositionalEnthalpy( string const & name,
                                              ComponentProperties const & componentProperties,
                                              integer const phaseIndex,
                                              ModelParameters const & modelParameters )
  : FunctionBase( name, componentProperties )
{
  EquationOfState const * equationOfState = modelParameters.get< EquationOfState >();
  string const eosName = equationOfState->m_equationsOfStateNames[phaseIndex];
  m_equationOfState = EnumStrings< EquationOfStateType >::fromString( eosName );
}

CompositionalEnthalpy::KernelWrapper
CompositionalEnthalpy::createKernelWrapper() const
{
  return KernelWrapper( m_equationOfState );
}

std::unique_ptr< ModelParameters >
CompositionalEnthalpy::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return EquationOfState::create( std::move( parameters ) );
}

} // namespace compositional
} // namespace constitutive
} // namespace geos
