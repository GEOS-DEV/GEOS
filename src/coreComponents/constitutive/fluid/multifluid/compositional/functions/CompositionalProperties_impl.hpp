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
 * @file CompositionalPropertiesImpl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_IMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_IMPL_HPP_

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
template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalProperties::computeMolarDensity( integer const numComps,
                                                   real64 const pressure,
                                                   real64 const temperature,
                                                   arraySlice1d< real64 const, USD1 > const & composition,
                                                   arraySlice1d< real64 const > const & volumeShift,
                                                   real64 const compressibilityFactor,
                                                   arraySlice1d< real64 const > const & compressibilityFactorDerivs,
                                                   real64 & molarDensity,
                                                   arraySlice1d< real64, USD2 > const & molarDensityDerivs )
{
  integer const numDofs = numComps + 2;

  real64 const vEos = constants::gasConstant * temperature * compressibilityFactor / pressure;
  auto const & vEosDerivs = molarDensityDerivs;
  for( integer dof = 0; dof < numDofs; ++dof )
  {
    vEosDerivs[dof] = constants::gasConstant * temperature * compressibilityFactorDerivs[dof] / pressure;
  }
  vEosDerivs[Deriv::dP] -= constants::gasConstant * temperature * compressibilityFactor / pressure / pressure;
  vEosDerivs[Deriv::dT] += constants::gasConstant * compressibilityFactor / pressure;

  real64 vCorrected = vEos;
  auto const & vCorrectedDerivs = molarDensityDerivs;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    vCorrected -= composition[ic] * volumeShift[ic];
    vCorrectedDerivs[Deriv::dC+ic] -= volumeShift[ic];
  }

  molarDensity = ( MultiFluidConstants::epsilon < vCorrected ) ? 1.0 / vCorrected : 0.0;
  for( integer dof = 0; dof < numDofs; ++dof )
  {
    molarDensityDerivs[dof] = -vCorrectedDerivs[dof]*molarDensity*molarDensity;
  }
}

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalProperties::computeMassDensity( integer const numComps,
                                                  arraySlice1d< real64 const, USD1 > const & composition,
                                                  arraySlice1d< real64 const > const & molecularWeight,
                                                  real64 const molarDensity,
                                                  arraySlice1d< real64 const, USD2 > const & molarDensityDerivs,
                                                  real64 & massDensity,
                                                  arraySlice1d< real64, USD2 > const & massDensityDerivs )
{
  massDensity = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    massDensity += molecularWeight[ic] * composition[ic] * molarDensity;
  }

  if( massDensity < MultiFluidConstants::epsilon )
  {
    auto const setZero = []( real64 & val ){ val = 0.0; };
    LvArray::forValuesInSlice( massDensityDerivs, setZero );
    return;
  }

  real64 const oneOverMolarDensity = 1.0 / molarDensity;

  // Pressure and temperature derivatives
  for( integer const kc : {Deriv::dP, Deriv::dT} )
  {
    massDensityDerivs[kc] = massDensity * molarDensityDerivs[kc] * oneOverMolarDensity;
  }

  // Composition derivative
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const kc = Deriv::dC + ic;
    massDensityDerivs[kc] = massDensity * molarDensityDerivs[kc] * oneOverMolarDensity + molecularWeight[ic] * molarDensity;
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_IMPL_HPP_
