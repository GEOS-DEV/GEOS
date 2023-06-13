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
 * @file CompositionalProperties.hpp
 */

#include "CompositionalProperties.hpp"

namespace geos
{

namespace constitutive
{

GEOS_HOST_DEVICE
void CompositionalProperties::computeMolarDensity( integer const numComps,
                                                   real64 const pressure,
                                                   real64 const temperature,
                                                   arraySlice1d< real64 const > const & composition,
                                                   arraySlice1d< real64 const > const & volumeShift,
                                                   real64 const compressibilityFactor,
                                                   real64 & molarDensity )
{

  real64 vEos = gasConstant * temperature * compressibilityFactor / pressure;
  real64 vCorrected = vEos;

  for( integer ic = 0; ic < numComps; ++ic )
  {
    vCorrected = vCorrected + composition[ic] * ( volumeShift[ic] * temperature + volumeShift[ic] );
  }

  if( epsilon < vCorrected )
  {
    molarDensity = 1.0 / vCorrected;
  }
  else
  {
    molarDensity = 0.0;
  }
}

GEOS_HOST_DEVICE
void CompositionalProperties::computeMolarDensity( integer const numComps,
                                                   real64 const pressure,
                                                   real64 const temperature,
                                                   arraySlice1d< real64 const > const & composition,
                                                   arraySlice1d< real64 const > const & volumeShift,
                                                   real64 const compressibilityFactor,
                                                   real64 const dCompressibilityFactor_dp,
                                                   real64 const dCompressibilityFactor_dt,
                                                   arraySlice1d< real64 const > const & dCompressibilityFactor_dz,
                                                   real64 const molarDensity,
                                                   real64 & dMolarDensity_dp,
                                                   real64 & dMolarDensity_dt,
                                                   arraySlice1d< real64 > const & dMolarDensity_dz )
{
  if( molarDensity < epsilon )
  {
    dMolarDensity_dp = 0.0;
    dMolarDensity_dt = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dMolarDensity_dz[ic] = 0.0;
    }
    return;
  }

  real64 dvCorrected_dx = 0.0;

  // Pressure derivative
  dvCorrected_dx = gasConstant * temperature * (dCompressibilityFactor_dp - compressibilityFactor / pressure) / pressure;
  dMolarDensity_dp = -molarDensity * molarDensity * dvCorrected_dx;

  // Temperature derivative
  dvCorrected_dx = gasConstant * (temperature * dCompressibilityFactor_dt + compressibilityFactor) / pressure;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dvCorrected_dx += composition[ic] * volumeShift[ic];
  }
  dMolarDensity_dt = -molarDensity * molarDensity * dvCorrected_dx;

  // Composition derivative
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dvCorrected_dx = gasConstant * temperature * dCompressibilityFactor_dz[ic] / pressure
                     + ( volumeShift[ic] * temperature + volumeShift[ic] );
    dMolarDensity_dz[ic] = -molarDensity * molarDensity * dvCorrected_dx;
  }
}

GEOS_HOST_DEVICE
void CompositionalProperties::computeMassDensity( integer const numComps,
                                                  arraySlice1d< real64 const > const & composition,
                                                  arraySlice1d< real64 const > const & molecularWeight,
                                                  real64 const molarDensity,
                                                  real64 & massDensity )
{
  massDensity = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    massDensity += molecularWeight[ic] * composition[ic] * molarDensity;
  }
}

GEOS_HOST_DEVICE
void CompositionalProperties::computeMassDensity( integer const numComps,
                                                  arraySlice1d< real64 const > const & molecularWeight,
                                                  real64 const molarDensity,
                                                  real64 const dMolarDensity_dp,
                                                  real64 const dMolarDensity_dt,
                                                  arraySlice1d< real64 const > const dMolarDensity_dz,
                                                  real64 const massDensity,
                                                  real64 & dMassDensity_dp,
                                                  real64 & dMassDensity_dt,
                                                  arraySlice1d< real64 > const & dMassDensity_dz )
{
  // Pressure derivative
  dMassDensity_dp = massDensity * dMolarDensity_dp / molarDensity;

  // Temperature derivative
  dMassDensity_dt = massDensity * dMolarDensity_dt / molarDensity;

  // Composition derivative
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dMassDensity_dz[ic] = massDensity * dMolarDensity_dz[ic] / molarDensity + molecularWeight[ic] * molarDensity;
  }
}

} // namespace constitutive

} // namespace geos
