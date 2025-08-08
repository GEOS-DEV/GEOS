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
  using CubicModel = CubicEOSPhaseModel< EOS_TYPE >;
  using Deriv = typename CubicModel::Deriv;

  template< bool DERIVATIVES = false >
  using StackVariables = SalinityStackVariables_Impl< void, DERIVATIVES >;

  /**
   * @brief Allocate and initialise composition independent data
   * @details Will allocate and initialise the data that is independent of the composition. This can be used in subsequent calls
   * @tparam DERIVATIVES a flag to indicate if derivatives (wrt p and/or t) should be allocated and calculated
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] data The composition data
   */
  template< bool DERIVATIVES = false >
  GEOS_HOST_DEVICE
  static void
  initialiseStack( integer const numComps, real64 const & pressure,
                   real64 const & temperature,
                   ComponentProperties::KernelWrapper const & componentProperties,
                   real64 const & salinity,
                   StackVariables< DERIVATIVES > & data );

  /**
   * @brief Main entry point of the Soreide-Whitson EOS model
   * @details Computes the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity the salinity
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
   * @param[in] salinity the salinity
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   * @param[out] logFugacityCoefficientDerivs derivatives of the log of the fugacity coefficients
   */
  template< integer USD1, integer USD2 = 0 >
  GEOS_HOST_DEVICE
  static void
  computeLogFugacityCoefficientsAndDerivs( integer const numComps,
                                           real64 const & pressure,
                                           real64 const & temperature,
                                           arraySlice1d< real64 const, USD1 > const & composition,
                                           ComponentProperties::KernelWrapper const & componentProperties,
                                           real64 const & salinity,
                                           arraySlice1d< real64 > const & logFugacityCoefficients,
                                           arraySlice2d< real64, USD2 > const & logFugacityCoefficientDerivs );

  /**
   * @brief Compute compressibility factor
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity the salinity
   * @param[out] compressibilityFactor the calculated compressibility factor
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
   * @param[in] salinity the salinity
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
   * @brief Calculate the pure coefficients derivatives for the water component
   * @tparam DERIVATIVES a flag to indicate if derivatives (wrt p and/or t) should be calculated
   * @param[in] h2o_index Water component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity the salinity
   * @param[out] data data The component mixture properties
   */
  template< bool DERIVATIVES = false >
  GEOS_HOST_DEVICE
  static void
  computeWaterCoefficients( integer const h2o_index,
                            real64 const & pressure,
                            real64 const & temperature,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            real64 const & salinity,
                            StackVariables< DERIVATIVES > & data );

  /**
   * @brief Calculate the pure component coefficients for the water component
   *
   */
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
