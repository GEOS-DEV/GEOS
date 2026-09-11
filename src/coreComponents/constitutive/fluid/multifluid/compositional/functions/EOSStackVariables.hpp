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
 * @file EOSStackVariables.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_EOSSTACKVARIABLES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_EOSSTACKVARIABLES_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief Temporary storage structure for Equation of State (EOS) computations.
 *
 * This templated structure holds intermediate variables used during EOS evaluations,
 * including mixture parameters and component-specific coefficients. It supports optional
 * derivative storage based on the `DERIVATIVES` flag.
 *
 * @tparam T Type used for numerical computations (e.g., float, double).
 * @tparam DERIVATIVES Boolean flag indicating whether derivative-related storage is enabled.
 */
template< typename T, bool DERIVATIVES >
struct EOSStackVariables_Impl
{
  /// Maximum number of components supported in the mixture.
  static constexpr integer maxNumComp = MultiFluidConstants::MAX_NUM_COMPONENTS;

  /**
   * @brief Type alias for derivative storage.
   *
   * @tparam DIM Dimensionality of the derivative (default is 1).
   */
  template< integer DIM = 1, integer USD = 0 >
  using DerivativeType = T *;

  /**
   * @brief Type alias for constant derivative storage.
   *
   * @tparam DIM Dimensionality of the derivative (default is 1).
   */
  template< integer DIM = 1, integer USD = 0 >
  using ConstDerivativeType = const T *;

  /**
   * @brief Constructor for initializing EOS stack variables.
   *
   * @param numComps Number of components in the mixture.
   * @param bip Binary interaction parameter matrix (k_ij).
   */
  GEOS_HOST_DEVICE
  EOSStackVariables_Impl( integer const numComps,
                          arraySlice2d< real64 const > const bip ):
    kij( bip ),
    m_data( 2, numComps ),
    aic( m_data[0] ),
    bic( m_data[1] )
  {}

  /// Binary interaction parameters (k_ij) between components.
  arraySlice2d< real64 const > const kij;

  /// Mixture parameter 'a' in cubic EOS
  real64 aMixture{0.0};

  /// Mixture parameter 'b' in cubic EOS.
  real64 bMixture{0.0};

  /// Internal storage for temporary EOS-related data.
  StackArray< real64, 2, 2 * maxNumComp > m_data;

  /// Component-specific 'a_i' coefficients.
  arraySlice1d< real64 > const aic;

  /// Component-specific 'b_i' coefficients.
  arraySlice1d< real64 > const bic;
};

/**
 * @brief Extension of EOS stack variables to include derivative information.
 *
 * This specialization of EOSStackVariables_Impl enables storage and access to
 * temperature and pressure derivatives of EOS parameters and component-specific
 * coefficients. It is used when derivative computations are required.
 *
 * @tparam T Type used for numerical computations (e.g., float, double).
 */
template< typename T >
struct EOSStackVariables_Impl< T, true > : public EOSStackVariables_Impl< T, false >
{
  /// Inherit the maximum number of components from the base class.
  using EOSStackVariables_Impl< T, false >::maxNumComp;

  /// Maximum number of degrees of freedom (components + temperature + pressure).
  static constexpr integer maxNumDof = maxNumComp + 2;

  /**
   * @brief Type alias for writable derivative arrays.
   *
   * @tparam DIM Dimensionality of the derivative (default is 1).
   */
  template< integer DIM = 1, integer USD = 0 >
  using DerivativeType = ArraySlice< real64, DIM, USD >;

  /**
   * @brief Type alias for read-only derivative arrays.
   *
   * @tparam DIM Dimensionality of the derivative (default is 1).
   */
  template< integer DIM = 1, integer USD = 0 >
  using ConstDerivativeType = ArraySlice< real64 const, DIM, USD >;

  /**
   * @brief Constructor for initializing EOS stack variables with derivatives.
   *
   * @param numComps Number of components in the mixture.
   * @param bip Binary interaction parameter matrix (k_ij).
   * @param dbip_dT Temperature derivatives of binary interaction parameters (dk_ij/dT).
   */
  GEOS_HOST_DEVICE
  EOSStackVariables_Impl( integer const numComps,
                          arraySlice2d< real64 const > const bip,
                          arraySlice2d< real64 const > const dbip_dT ):
    EOSStackVariables_Impl< T, false >( numComps, bip ),
    dkij_dT( dbip_dT ),
    m_derivativeData( 9, numComps+2 ),
    daic_dp( m_derivativeData[0] ),
    dbic_dp( m_derivativeData[1] ),
    daic_dt( m_derivativeData[2] ),
    dbic_dt( m_derivativeData[3] ),
    d2aic_dt2( m_derivativeData[4] ),
    d2bic_dt2( m_derivativeData[5] ),
    daMixture( m_derivativeData[6] ),
    dbMixture( m_derivativeData[7] ),
    ddaMixture( m_derivativeData[8] )
  {}

  /// Temperature derivatives of binary interaction parameters (dk_ij/dT).
  arraySlice2d< real64 const > const dkij_dT;

  /// Internal storage for all derivative-related temporary data.
  StackArray< real64, 2, 9 * maxNumDof > m_derivativeData;

  /// Pressure derivatives of component-specific 'a_i' coefficients (da_i/dp).
  arraySlice1d< real64 > const daic_dp;

  /// Pressure derivatives of component-specific 'b_i' coefficients (db_i/dp).
  arraySlice1d< real64 > const dbic_dp;

  /// Temperature derivatives of component-specific 'a_i' coefficients (da_i/dT).
  arraySlice1d< real64 > const daic_dt;

  /// Temperature derivatives of component-specific 'b_i' coefficients (db_i/dT).
  arraySlice1d< real64 > const dbic_dt;

  /// Second temperature derivatives of 'a_i' coefficients (d2a_i/dT2).
  arraySlice1d< real64 > const d2aic_dt2;

  /// Second temperature derivatives of 'b_i' coefficients (d2b_i/dT2).
  arraySlice1d< real64 > const d2bic_dt2;

  /// Derivatives of mixture parameter 'a' with respect to all degrees of freedom.
  DerivativeType<> const daMixture;

  /// Derivatives of mixture parameter 'b' with respect to all degrees of freedom.
  DerivativeType<> const dbMixture;

  /// Derivatives of temperature derivatives of mixture parameter 'a' with respect to all degrees of freedom.
  DerivativeType<> const ddaMixture;
};

/**
 * @brief Container for storing binary interaction coefficients (BICs) used in EOS models.
 *
 * This structure holds temperature- and pressure-dependent binary interaction parameters
 * (e.g., \( k_{ij} \)) for a given number of components. These parameters are typically used
 * in cubic equations of state to model non-ideal interactions between fluid components.
 *
 * @tparam T Type used for numerical computations (e.g., float, double).
 * @tparam DERIVATIVES Boolean flag indicating whether derivative support is enabled.
 */
template< typename T, bool DERIVATIVES >
struct BICStackVariables_Impl
{
  /// Maximum number of components supported by the EOS model.
  static constexpr integer maxNumComp = EOSStackVariables_Impl< T, DERIVATIVES >::maxNumComp;

  /**
   * @brief Constructor for initializing the BIC storage structure.
   *
   * @param numComps Number of components in the fluid mixture.
   */
  GEOS_HOST_DEVICE
  BICStackVariables_Impl( integer numComps ):
    kij_data( numComps, numComps )
  {}

  /// Salinity value used in BIC calculations (if applicable).
  real64 salinity{0.0};

  /**
   * @brief Binary interaction coefficients \( k_{ij} \), stored in a flat 2D array.
   *
   * The array is sized to hold all pairwise interactions between components.
   * These coefficients may vary with temperature and pressure.
   */
  StackArray< real64, 2, maxNumComp * maxNumComp > kij_data;
};

/**
 * @brief Specialization of BICStackVariables_Impl to include temperature derivatives.
 *
 * This structure extends the base BIC container to store the first-order temperature
 * derivatives of binary interaction coefficients \( \frac{\partial k_{ij}}{\partial T} \),
 * which are used in thermodynamic models requiring sensitivity to temperature.
 *
 * @tparam T Type used for numerical computations (e.g., float, double).
 */
template< typename T >
struct BICStackVariables_Impl< T, true > : public BICStackVariables_Impl< T, false >
{
  /// Inherit the maximum number of components from the base class.
  using BICStackVariables_Impl< T, false >::maxNumComp;

  /**
   * @brief Constructor for initializing BIC storage with temperature derivatives.
   *
   * @param numComps Number of components in the fluid mixture.
   */
  GEOS_HOST_DEVICE
  BICStackVariables_Impl( integer numComps ):
    BICStackVariables_Impl< T, false >( numComps ),
    dkij_dT_data( numComps, numComps )
  {}

  /**
   * @brief Temperature derivatives of binary interaction coefficients \( \frac{\partial k_{ij}}{\partial T} \).
   *
   * Stored in a flat 2D array sized for all component pairs.
   */
  StackArray< real64, 2, maxNumComp * maxNumComp > dkij_dT_data;
};

/**
 * @brief Generic template declaration for salinity-aware stack variables.
 *
 * This template is specialized for `DERIVATIVES = false` and `DERIVATIVES = true`.
 * It combines EOS and BIC temporary variables for thermodynamic models that
 * incorporate salinity effects.
 *
 * @tparam T Type used for numerical computations (e.g., float, double).
 * @tparam DERIVATIVES Boolean flag indicating whether derivative support is enabled.
 */
template< typename T, bool DERIVATIVES >
struct SalinityStackVariables_Impl {};

/**
 * @brief Specialization for salinity-aware stack variables without derivatives.
 *
 * This version includes EOS and BIC temporary variables for salinity-dependent
 * thermodynamic models, but does not store or compute derivatives.
 *
 * @tparam T Type used for numerical computations.
 */
template< typename T >
struct SalinityStackVariables_Impl< T, false > : public BICStackVariables_Impl< T, false >,
  public EOSStackVariables_Impl< T, false >
{
  /**
   * @brief Constructor for initializing salinity-aware variables (non-derivative version).
   *
   * @param numComps Number of components in the fluid mixture.
   */
  GEOS_HOST_DEVICE
  SalinityStackVariables_Impl( integer numComps ):
    BICStackVariables_Impl< T, false >( numComps ),
    EOSStackVariables_Impl< T, false >( numComps,
                                        BICStackVariables_Impl< T, false >::kij_data.toSliceConst() )
  {}
};

/**
 * @brief Specialization for salinity-aware stack variables with derivative support.
 *
 * This version includes EOS and BIC temporary variables along with their
 * temperature and pressure derivatives, used in sensitivity analysis or
 * Jacobian evaluations in salinity-influenced EOS models.
 *
 * @tparam T Type used for numerical computations.
 */
template< typename T >
struct SalinityStackVariables_Impl< T, true > : public BICStackVariables_Impl< T, true >,
  public EOSStackVariables_Impl< T, true >
{
  /**
   * @brief Constructor for initializing salinity-aware variables (derivative version).
   *
   * @param numComps Number of components in the fluid mixture.
   */
  GEOS_HOST_DEVICE
  SalinityStackVariables_Impl( integer numComps ):
    BICStackVariables_Impl< T, true >( numComps ),
    EOSStackVariables_Impl< T, true >( numComps,
                                       BICStackVariables_Impl< T, true >::kij_data.toSliceConst(),
                                       BICStackVariables_Impl< T, true >::dkij_dT_data.toSliceConst() )
  {}
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_EOSSTACKVARIABLES_HPP_
