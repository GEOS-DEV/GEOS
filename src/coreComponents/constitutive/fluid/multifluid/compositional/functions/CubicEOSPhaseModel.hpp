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
 * @file CubicEOSPhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_

#include "EOSStackVariables.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct PengRobinsonEOS
{
  static constexpr real64 omegaA = 0.457235529;
  static constexpr real64 omegaB = 0.077796074;
  static constexpr real64 delta1 = 2.4142135624; // 1 + sqrt( 2 )
  static constexpr real64 delta2 = -0.4142135624; // 1 - sqrt( 2 )

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64
  kappa( real64 const & omega )
  {
    return ( omega < 0.49 )
      ? 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
      : 0.3796 + 1.485 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
  }
};

struct SoaveRedlichKwongEOS
{
  static constexpr real64 omegaA = 0.42748;
  static constexpr real64 omegaB = 0.08664;
  static constexpr real64 delta1 = 0.0;
  static constexpr real64 delta2 = 1.0;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64
  kappa( real64 const & omega )
  {
    return 0.480 + 1.574 * omega - 0.176 * omega * omega;
  }
};

template< typename EOS_TYPE >
struct CubicEOSPhaseModel
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using CubicModel = CubicEOSPhaseModel< EOS_TYPE >;

  // Enumeration for selected root
  enum class SelectedRoot : int
  {
    AUTO = 0,
    MINIMUM = 1,
    MAXIMUM = 2
  };

  template< bool DERIVATIVES >
  using StackVariables = EOSStackVariables_Impl< void, DERIVATIVES >;

  template< integer DIM, bool DERIVATIVES, integer USD = 0 >
  using StackDerivativeType = typename StackVariables< DERIVATIVES >::template DerivativeType< DIM, USD >;

  template< integer DIM, bool DERIVATIVES, integer USD = 0 >
  using StackConstDerivativeType = typename StackVariables< DERIVATIVES >::template ConstDerivativeType< DIM, USD >;

public:
  /**
   * @brief Allocate and initialise composition independent data
   * @details Will allocate and initialise the data that is independent of the composition. This can be used in subsequent calls
   * @tparam DERIVATIVES a flag to indicate if derivatives (wrt p and/or t) should be allocated and calculated
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] data The composition data
   */
  template< bool DERIVATIVES = false >
  GEOS_HOST_DEVICE
  static void
  initialiseStack( integer const numComps,
                   real64 const & pressure,
                   real64 const & temperature,
                   ComponentProperties::KernelWrapper const & componentProperties,
                   StackVariables< DERIVATIVES > & data );

  /**
   * @brief Main entry point of the cubic EOS model
   * @details Computes the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
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
                                  arraySlice1d< real64 > const & logFugacityCoefficients );

  /**
   * @brief Secondary entry point of the cubic EOS model
   * @details Computes the derivatives of the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   * @param[out] logFugacityCoefficientDerivs derivatives of the log of the fugacity coefficients
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static void
  computeLogFugacityCoefficientsAndDerivs( integer const numComps,
                                           real64 const & pressure,
                                           real64 const & temperature,
                                           arraySlice1d< real64 const, USD1 > const & composition,
                                           ComponentProperties::KernelWrapper const & componentProperties,
                                           arraySlice1d< real64 > const & logFugacityCoefficients,
                                           arraySlice2d< real64, USD2 > const & logFugacityCoefficientDerivs );

  /**
   * @brief Compute compressibility factor for the cubic EOS model
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] compressibilityFactor the current compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 & compressibilityFactor );

  /**
   * @brief Compute compressibility factor for the cubic EOS model
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] compressibilityFactor the current compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactorAndDerivs( integer const numComps,
                                         real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD > const & composition,
                                         ComponentProperties::KernelWrapper const & componentProperties,
                                         real64 & compressibilityFactor,
                                         arraySlice1d< real64 > const & compressibilityFactorDerivs );

  /**
   * @brief Calculate the dimensional volume shift
   * @details Computes the dimensional form of the volume shifts given the user defined non-dimensional form.
   * @param[in] numComps The number of components
   * @param[in] componentProperties The compositional model properties
   * @param[out] dimensionalVolumeShift The calculated dimensional volume shifts
   */
  GEOS_FORCE_INLINE
  static void calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                               arraySlice1d< real64 > const & dimensionalVolumeShift );

  /**
   * @brief Calculate the pure coefficients derivatives
   * @tparam DERIVATIVES a flag to indicate if derivatives (wrt p and/or t) should be calculated
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] data data The component mixture properties
   */
  template< bool DERIVATIVES = false >
  GEOS_HOST_DEVICE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           StackVariables< DERIVATIVES > & data );

  /**
   * @brief Compute the mixture coefficients using pressure, temperature, composition and input
   * @tparam DERIVATIVES a flag to indicate if derivatives should be calculated
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in/out] data The component mixture properties
   */
  template< integer USD, bool DERIVATIVES = false >
  GEOS_HOST_DEVICE
  static void
  computeMixtureCoefficients( integer const numComps,
                              arraySlice1d< real64 const, USD > const & composition,
                              StackVariables< DERIVATIVES > & data );

  /**
   * @brief Compute the compressibility factor using compositions, BICs, and mixture coefficients
   * @tparam DERIVATIVES a flag to indicate if derivatives should be calculated
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in] data The component mixture properties
   * @param[out] compressibilityFactor compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   * @param[in] selectedRoot Indicates which root is to be preferred in the case of 3 roots
   */
  template< integer USD, bool DERIVATIVES = false >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactor( integer const numComps,
                                arraySlice1d< real64 const, USD > const & composition,
                                StackVariables< DERIVATIVES > const & data,
                                real64 & compressibilityFactor,
                                StackDerivativeType< 1, DERIVATIVES > const & compressibilityFactorDerivs,
                                SelectedRoot const selectedRoot = SelectedRoot::AUTO );

  /**
   * @brief Compute the log of the fugacity coefficients using compositions, BICs, compressibility factor and mixture coefficients
   * @tparam DERIVATIVES a flag to indicate if derivatives should be calculated
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in] data The component mixture properties
   * @param[in] compressibilityFactor compressibility factor
   * @param[in] compressibilityFactorDerivs derivatives of the compressibility factor
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   * @param[out] logFugacityCoefficientDerivs derivatives of the log of the fugacity coefficients
   */
  template< integer USD1, bool DERIVATIVES = false, integer USD2 = 0 >
  GEOS_HOST_DEVICE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  arraySlice1d< real64 const, USD1 > const & composition,
                                  StackVariables< DERIVATIVES > const & data,
                                  real64 const & compressibilityFactor,
                                  StackConstDerivativeType< 1, DERIVATIVES > const & compressibilityFactorDerivs,
                                  arraySlice1d< real64 > const & logFugacityCoefficients,
                                  StackDerivativeType< 2, DERIVATIVES, USD2 > const & logFugacityCoefficientDerivs );

  /**
   * @brief Helper functions solving a cubic equation using trigonometry
   *        m3 * x^3 + m2 * x^2 + m1 *x + m0  = 0
   * @param[in] m3 first coefficient (in front of x^3)
   * @param[in] m2 second coefficient (in front of x^2)
   * @param[in] m1 third coefficient (in front of x)
   * @param[in] m0 fourth coefficient
   * @param[out] roots the roots of the polynomial
   * @param[out] numRoots the number of roots
   */
  GEOS_HOST_DEVICE
  static void solveCubicPolynomial( real64 const & m3,
                                    real64 const & m2,
                                    real64 const & m1,
                                    real64 const & m0,
                                    real64 ( &roots )[3],
                                    integer & numRoots );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setZero( real64 & val ){ val = 0.0; }
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

// Include implementation
#include "CubicEOSPhaseModel_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
