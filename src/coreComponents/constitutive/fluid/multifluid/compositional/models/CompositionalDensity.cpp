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

CompositionalDensity::KernelWrapper
CompositionalDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentDimensionalVolumeShift, m_equationOfState );
}

std::unique_ptr< ModelParameters >
CompositionalDensity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return EquationOfState::create( std::move( parameters ) );
}

void CompositionalDensity::calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                                            EquationOfStateType const & equationOfState,
                                                            arraySlice1d< real64 > componentDimensionalVolumeShift )
{
  if( equationOfState == EquationOfStateType::PengRobinson )
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                            componentDimensionalVolumeShift );
  }
  else if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                                 componentDimensionalVolumeShift );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
