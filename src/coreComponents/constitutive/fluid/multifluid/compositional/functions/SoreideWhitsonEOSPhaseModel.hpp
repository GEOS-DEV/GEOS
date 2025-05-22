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
 * @file SoreideWhitsonEOSPhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_HPP_

#include "CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE >
struct SoreideWhitsonEOSPhaseModel
{
private:
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;

public:
  using CubicModel = CubicEOSPhaseModel< EOS_TYPE >;
  using Deriv = typename CubicModel::Deriv;

  /**
   * @brief Main entry point of the Soreide-Whitson EOS model
   * @details Computes the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  real64 const & salinity,
                                  arraySlice1d< real64 > const & logFugacityCoefficients );

  /**
   * @brief Secondary entry point of the Soreide-Whitson EOS model
   * @details Computes the derivatives of the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] logFugacityCoefficients log of the fugacity coefficients
   * @param[out] logFugacityCoefficientDerivs derivatives of the log of the fugacity coefficients
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeLogFugacityCoefficientDerivs( integer const numComps,
                                       real64 const & pressure,
                                       real64 const & temperature,
                                       arraySlice1d< real64 const, USD > const & composition,
                                       ComponentProperties::KernelWrapper const & componentProperties,
                                       real64 const & salinity,
                                       arraySlice1d< real64 const > const & logFugacityCoefficients,
                                       arraySlice2d< real64 > const & logFugacityCoefficientDerivs );

  /**
   * @brief Calculate the pure coefficients
   * @details Computes the pure coefficients
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   */
  GEOS_HOST_DEVICE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 const & salinity,
                           real64 & aCoefficient,
                           real64 & bCoefficient );

  /**
   * @brief Calculate the pure coefficients derivatives
   * @details Computes the pure coefficients derivatives
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   * @param[out] daCoefficient_dp pure coefficient (A) derivative w.r.t. pressure
   * @param[out] dbCoefficient_dp pure coefficient (B) derivative w.r.t. pressure
   * @param[out] daCoefficient_dt pure coefficient (A) derivative w.r.t. temperature
   * @param[out] dbCoefficient_dt pure coefficient (B) derivative w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  static void
  computePureCoefficientsAndDerivs( integer const ic,
                                    real64 const & pressure,
                                    real64 const & temperature,
                                    ComponentProperties::KernelWrapper const & componentProperties,
                                    real64 const & salinity,
                                    real64 & aCoefficient,
                                    real64 & bCoefficient,
                                    real64 & daCoefficient_dp,
                                    real64 & dbCoefficient_dp,
                                    real64 & daCoefficient_dt,
                                    real64 & dbCoefficient_dt );

  /**
   * @brief Compute the mixture coefficients using pressure, temperature, composition and input
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] aPureCoefficient pure coefficient (A)
   * @param[out] bPureCoefficient pure coefficient (B)
   * @param[out] aMixtureCoefficient mixture coefficient (A)
   * @param[out] bMixtureCoefficient mixture coefficient (B)
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 const & salinity,
                              arraySlice1d< real64 > const & aPureCoefficient,
                              arraySlice1d< real64 > const & bPureCoefficient,
                              real64 & aMixtureCoefficient,
                              real64 & bMixtureCoefficient );

  /**
   * @brief Compute the mixture coefficients derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] aMixtureCoefficientDerivs derivatives of mixture coefficient (A)
   * @param[out] bMixtureCoefficientDerivs derivatives of mixture coefficient (B)
   * @note Assumes that pressure and temperature are strictly positive
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeMixtureCoefficientDerivs( integer const numComps,
                                   real64 const & pressure,
                                   real64 const & temperature,
                                   arraySlice1d< real64 const, USD > const & composition,
                                   ComponentProperties::KernelWrapper const & componentProperties,
                                   real64 const & salinity,
                                   arraySlice1d< real64 const > const & aPureCoefficient,
                                   arraySlice1d< real64 const > const & bPureCoefficient,
                                   real64 const aMixtureCoefficient,
                                   real64 const bMixtureCoefficient,
                                   arraySlice1d< real64 > const & aMixtureCoefficientDerivs,
                                   arraySlice1d< real64 > const & bMixtureCoefficientDerivs );

  /**
   * @brief Compute compressibility factor
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] compressibilityFactor the current compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 const & salinity,
                                real64 & compressibilityFactor );

  /**
   * @brief Compute compressibility factor and derivatives
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] compressibilityFactor the current compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactorAndDerivs( integer const numComps,
                                         real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD1 > const & composition,
                                         ComponentProperties::KernelWrapper const & componentProperties,
                                         real64 const & salinity,
                                         real64 & compressibilityFactor,
                                         arraySlice1d< real64, USD2 > const & compressibilityFactorDerivs );

  /**
   * @brief Get the binary interaction coefficient between two components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[in] i index of first component
   * @param[in] j index of second component
   * @param[out] kij the binary interaction coefficient
   * @param[out] dkij_dT derivative of the binary interaction coefficient w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  static void
  getBinaryInteractionCoefficient( real64 const & pressure,
                                   real64 const & temperature,
                                   ComponentProperties::KernelWrapper const & componentProperties,
                                   real64 const & salinity,
                                   integer const i,
                                   integer const j,
                                   real64 & kij,
                                   real64 & dkij_dT );

  /**
   * @brief Safe pow function
   * @param[in] a The value
   * @param[in] b The exponent
   * @return returns a**b if a is positive else 0
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 power( real64 const a, real64 const b )
  {
    return a < MultiFluidConstants::minForSpeciesPresence ? 0.0 : std::pow( a, b );
  }
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

// Include the implementation
#include "SoreideWhitsonEOSPhaseModel_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_HPP_
