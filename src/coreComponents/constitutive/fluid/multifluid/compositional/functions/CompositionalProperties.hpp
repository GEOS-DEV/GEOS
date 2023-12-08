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

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct CompositionalProperties
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Compute the molar density of a mixture from the composition and the compressibility factor
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] volumeShift dimensional volume shift parameters
   * @param[in] compressibilityFactor compressibility factor (z-factor)
   * @param[out] molarDensity the calculated molar density
   * @note The volume shifts can result in a negative molar density which will be truncated to zero
   */
  GEOS_HOST_DEVICE
  static void computeMolarDensity( integer const numComps,
                                   real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const > const & composition,
                                   arraySlice1d< real64 const > const & volumeShift,
                                   real64 const compressibilityFactor,
                                   real64 & molarDensity );

  /**
   * @brief Compute the molar density derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] volumeShift dimensional volume shift parameters
   * @param[in] compressibilityFactor compressibility factor (z-factor)
   * @param[in] compressibilityFactorDerivs derivatives of the compressibility factor (z-factor)
   * @param[in] molarDensity the calculated molar density
   * @param[out] molarDensityDerivs derivatives of the molar density
   */
  GEOS_HOST_DEVICE
  static void computeMolarDensity( integer const numComps,
                                   real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const > const & composition,
                                   arraySlice1d< real64 const > const & volumeShift,
                                   real64 const compressibilityFactor,
                                   arraySlice1d< real64 const > const & compressibilityFactorDerivs,
                                   real64 const molarDensity,
                                   arraySlice1d< real64 > const & molarDensityDerivs );

  /**
   * @brief Compute the mass density of a mixture from the composition and the molar density
   * @param[in] numComps number of components
   * @param[in] composition composition of the mixture
   * @param[in] molecularWeight the component molecular weights
   * @param[in] molarDensity the mixture molar density
   * @param[out] massDensity the calculated mass density
   */
  GEOS_HOST_DEVICE
  static void computeMassDensity( integer const numComps,
                                  arraySlice1d< real64 const > const & composition,
                                  arraySlice1d< real64 const > const & molecularWeight,
                                  real64 const molarDensity,
                                  real64 & massDensity );

  /**
   * @brief Compute the mass density derivatives
   * @param[in] numComps number of components
   * @param[in] molecularWeight the component molecular weights
   * @param[in] molarDensity the mixture molar density
   * @param[in] molarDensityDerivs derivatives of the molar density
   * @param[in] massDensity mass density
   * @param[out] massDensityDerivs derivatives of the mass density
   */
  GEOS_HOST_DEVICE
  static void computeMassDensity( integer const numComps,
                                  arraySlice1d< real64 const > const & molecularWeight,
                                  real64 const molarDensity,
                                  arraySlice1d< real64 const > const & molarDensityDerivs,
                                  real64 const massDensity,
                                  arraySlice1d< real64 > const & massDensityDerivs );
};

} //namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_
