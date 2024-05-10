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
 * @file CompositionalDensity.cpp
 */

#include "CompositionalDensity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

CompositionalDensity::CompositionalDensity( string const & name,
                                            ComponentProperties const & componentProperties,
                                            EquationOfState const & equationOfState,
                                            integer const phaseIndex ):
  FunctionBase( name, componentProperties, equationOfState, phaseIndex ),
  m_phaseIndex( phaseIndex )
{
  // Calculate the dimensional volume shift
  m_componentDimensionalVolumeShift.resize( componentProperties.getNumberOfComponents());
  calculateDimensionalVolumeShift( componentProperties,
                                   equationOfState,
                                   m_componentDimensionalVolumeShift );
}

CompositionalDensityUpdate::CompositionalDensityUpdate( arrayView1d< real64 const > const & volumeShift,
                                                        integer const eosType ):
  m_componentDimensionalVolumeShift( volumeShift ),
  m_eosType( eosType )
{}

CompositionalDensity::KernelWrapper CompositionalDensity::createKernelWrapper() const
{
  integer const eosType = m_equationOfState.getEquationOfStateType( m_phaseIndex );
  return KernelWrapper( m_componentDimensionalVolumeShift, eosType );
}

void CompositionalDensity::calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                                            EquationOfState const & equationOfState,
                                                            arraySlice1d< real64 > componentDimensionalVolumeShift )
{
  integer const eosType = equationOfState.getEquationOfStateType( m_phaseIndex );
  if( EquationOfState::equals( eosType, EquationOfStateType::PengRobinson ))
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                            componentDimensionalVolumeShift );
  }
  else if( EquationOfState::equals( eosType, EquationOfStateType::SoaveRedlichKwong ))
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                                 componentDimensionalVolumeShift );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
