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
 * @file EquationOfState.cpp
 */

#include "EquationOfState.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

EquationOfState::EquationOfState( string_array const & names ):
  m_names( names )
{
  convertNames( m_names, m_types );
}

EquationOfState::KernelWrapper::KernelWrapper( arrayView1d< integer const > const & types )
  : m_types( types )
{}

EquationOfState::KernelWrapper
EquationOfState::createKernelWrapper() const
{
  return KernelWrapper( m_types );
}

integer EquationOfState::getEquationOfStateType( integer const & phaseIndex ) const
{
  return m_types[phaseIndex];
}

bool EquationOfState::validateNames( string_array const & names )
{
  for( string const & name : names )
  {
    EnumStrings< EquationOfStateType >::fromString( name );
  }
  return true;
}

void EquationOfState::convertNames( string_array const & names, array1d< integer > & types )
{
  types.resize( names.size());
  for( integer i = 0; i < names.size(); ++i )
  {
    EquationOfStateType const e = EnumStrings< EquationOfStateType >::fromString( names[i] );
    types[i] = static_cast< int >(e);
  }
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
