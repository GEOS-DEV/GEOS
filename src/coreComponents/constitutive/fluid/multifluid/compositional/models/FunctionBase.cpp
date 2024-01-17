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
 * @file FunctionBase.cpp
 */

#include "FunctionBase.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

GEOS_HOST_DEVICE
void
FunctionBaseUpdate::convertDerivativesToTotalMoleFraction( integer const numComps,
                                                           arraySlice2d< real64 const > const & dPhaseComposition,
                                                           arraySlice1d< real64 > const & dProperty,
                                                           arraySlice1d< real64 > const & workSpace )
{
  using Deriv = multifluid::DerivativeOffset;
  integer const numDofs = numComps + 2;
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    workSpace[kc] = dProperty[kc];
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dProperty[Deriv::dC+ic] = 0.0;
  }
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dProperty[kc] += (dPhaseComposition( ic, kc ) * workSpace[Deriv::dC+ic]);
    }
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
