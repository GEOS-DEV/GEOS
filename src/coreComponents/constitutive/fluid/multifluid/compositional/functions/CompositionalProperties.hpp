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

namespace geos
{

namespace constitutive
{

struct CompositionalProperties
{
public:
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
  /// Universal gas constant
  static constexpr real64 gasConstant = 8.31446261815324;

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
   * @param[in] dCompressibilityFactor_dp derivative of the compressibility factor (z-factor) wrt pressure
   * @param[in] dCompressibilityFactor_dp derivative of the compressibility factor (z-factor) wrt temperature
   * @param[in] dCompressibilityFactor_dz derivative of the compressibility factor (z-factor) wrt composition
   * @param[in] molarDensity the calculated molar density
   * @param[out] dMolarDensity_dp derivative of the molar density wrt pressure
   * @param[out] dMolarDensity_dt derivative of the molar density wrt temperature
   * @param[out] dMolarDensity_dz derivative of the molar density wrt composition
   */
  GEOS_HOST_DEVICE
  static void computeMolarDensity( integer const numComps,
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
                                   arraySlice1d< real64 > const & dMolarDensity_dz );

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
   * @param[in] dMolarDensity_dp derivative of the molar density wrt pressure
   * @param[in] dMolarDensity_dt derivative of the molar density wrt temperature
   * @param[in] dMolarDensity_dz derivative of the molar density wrt composition
   * @param[in] massDensity mass density
   * @param[out] dMassDensity_dp derivative of the mass density wrt pressure
   * @param[out] dMassDensity_dt derivative of the mass density wrt temperature
   * @param[out] dMassDensity_dz derivative of the mass density wrt composition
   */
  GEOS_HOST_DEVICE
  static void computeMassDensity( integer const numComps,
                                  arraySlice1d< real64 const > const & molecularWeight,
                                  real64 const molarDensity,
                                  real64 const dMolarDensity_dp,
                                  real64 const dMolarDensity_dt,
                                  arraySlice1d< real64 const > const dMolarDensity_dz,
                                  real64 const massDensity,
                                  real64 & dMassDensity_dp,
                                  real64 & dMassDensity_dt,
                                  arraySlice1d< real64 > const & dMassDensity_dz );
};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_
