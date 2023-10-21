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
 * @file CubicEOSDensityImpl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CUBICEOSDENSITYIMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CUBICEOSDENSITYIMPL_HPP_

#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE >
template< int USD1 >
GEOS_HOST_DEVICE
void CubicEOSDensityUpdate< EOS_TYPE >::compute( real64 const & pressure,
                                                 real64 const & temperature,
                                                 arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                                 real64 & value,
                                                 bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure,
                   temperature,
                   phaseComposition,
                   useMass );

  value = 100.0;
}

template< typename EOS_TYPE >
template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
void CubicEOSDensityUpdate< EOS_TYPE >::compute( real64 const & pressure,
                                                 real64 const & temperature,
                                                 arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                                 arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                                 real64 & value,
                                                 arraySlice1d< real64, USD3 > const & dValue,
                                                 bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure,
                   temperature,
                   phaseComposition,
                   dPhaseComposition,
                   useMass );

  value = 100.0;

  LvArray::forValuesInSlice( dValue, []( real64 & val ){ val = 0.0; } );
}

using CubicEOSDensityPR = CubicEOSDensity< PengRobinsonEOS >;
using CubicEOSDensitySRK = CubicEOSDensity< SoaveRedlichKwongEOS >;

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CUBICEOSDENSITYIMPL_HPP_
