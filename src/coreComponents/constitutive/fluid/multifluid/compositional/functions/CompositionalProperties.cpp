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
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/*
 * Calculate the molar volume and apply the Peneloux shift parameters. The parameters should be in
 * dimensional form.
 *   Peneloux, A et al. 1982. Fluid phase equilibria, 8(1):7–23.
 *   https://doi.org/10.1016/0378-3812(82)80002-2
 */
GEOS_HOST_DEVICE
void CompositionalProperties::computeMolarDensity( integer const numComps,
                                                   real64 const pressure,
                                                   real64 const temperature,
                                                   arraySlice1d< real64 const > const & composition,
                                                   arraySlice1d< real64 const > const & volumeShift,
                                                   real64 const compressibilityFactor,
                                                   real64 & molarDensity )
{

  real64 vEos = MultiFluidConstants::gasConstant * temperature * compressibilityFactor / pressure;
  real64 vCorrected = vEos;

  for( integer ic = 0; ic < numComps; ++ic )
  {
    vCorrected += composition[ic] * volumeShift[ic];
  }

  if( MultiFluidConstants::epsilon < vCorrected )
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
                                                   arraySlice1d< real64 const > const & GEOS_UNUSED_PARAM ( composition ),
                                                   arraySlice1d< real64 const > const & volumeShift,
                                                   real64 const compressibilityFactor,
                                                   arraySlice1d< real64 const > const & compressibilityFactorDerivs,
                                                   real64 const molarDensity,
                                                   arraySlice1d< real64 > const & molarDensityDerivs )
{
  if( molarDensity < MultiFluidConstants::epsilon )
  {
    auto const setZero = []( real64 & val ){ val = 0.0; };
    LvArray::forValuesInSlice( molarDensityDerivs, setZero );
    return;
  }

  real64 dvCorrected_dx = 0.0;

  // Pressure derivative
  dvCorrected_dx = MultiFluidConstants::gasConstant * temperature * (compressibilityFactorDerivs[Deriv::dP] - compressibilityFactor / pressure) / pressure;
  molarDensityDerivs[Deriv::dP] = -molarDensity * molarDensity * dvCorrected_dx;

  // Temperature derivative
  dvCorrected_dx = MultiFluidConstants::gasConstant * (temperature * compressibilityFactorDerivs[Deriv::dT] + compressibilityFactor) / pressure;
  molarDensityDerivs[Deriv::dT] = -molarDensity * molarDensity * dvCorrected_dx;

  // Composition derivative
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const kc = Deriv::dC + ic;
    dvCorrected_dx = MultiFluidConstants::gasConstant * temperature * compressibilityFactorDerivs[kc] / pressure + volumeShift[ic];
    molarDensityDerivs[kc] = -molarDensity * molarDensity * dvCorrected_dx;
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
                                                  arraySlice1d< real64 const > const & molarDensityDerivs,
                                                  real64 const massDensity,
                                                  arraySlice1d< real64 > const & massDensityDerivs )
{
  // Pressure and temperature derivatives
  for( integer const kc : {Deriv::dP, Deriv::dT} )
  {
    massDensityDerivs[kc] = massDensity * molarDensityDerivs[kc] / molarDensity;
  }

  // Composition derivative
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const kc = Deriv::dC + ic;
    massDensityDerivs[kc] = massDensity * molarDensityDerivs[kc] / molarDensity + molecularWeight[ic] * molarDensity;
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
