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
 * @file SoreideWhitsonEOSPhaseModel_impl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_IMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_IMPL_HPP_

#include "SoreideWhitsonEOSPhaseModel.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentType.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE >
template< bool DERIVATIVES >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
initialiseStack( integer const numComps,
                 real64 const & pressure,
                 real64 const & temperature,
                 ComponentProperties::KernelWrapper const & componentProperties,
                 real64 const & salinity,
                 StackVariables< DERIVATIVES > & stack )
{
  typename CubicModel::template StackVariables< DERIVATIVES > & cubicStack = stack;
  CubicModel::initialiseStack( numComps,
                               pressure,
                               temperature,
                               componentProperties,
                               cubicStack );

  stack.salinity = salinity;
  arraySlice1d< integer const > const & componentType = componentProperties.m_componentType;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    if( isComponentType( componentType[ic], ComponentType::Water ))
    {
      computeWaterCoefficients( ic,
                                pressure,
                                temperature,
                                componentProperties,
                                salinity,
                                stack );
    }
  }

  // Calculate the binary interaction coefficients
  double kij = 0.0;
  double dkij_dT = 0.0;
  if constexpr (DERIVATIVES)
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      stack.kij_data( ic, ic ) = 0.0;
      stack.dkij_dT_data( ic, ic ) = 0.0;
      for( integer jc = 0; jc < ic; ++jc )
      {
        getBinaryInteractionCoefficient( pressure,
                                         temperature,
                                         componentProperties,
                                         salinity,
                                         ic,
                                         jc,
                                         kij,
                                         dkij_dT );
        stack.kij_data( ic, jc ) = kij;
        stack.kij_data( jc, ic ) = kij;
        stack.dkij_dT_data( ic, jc ) = dkij_dT;
        stack.dkij_dT_data( jc, ic ) = dkij_dT;
      }
    }
  }
  else
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      stack.kij_data( ic, ic ) = 0.0;
      for( integer jc = 0; jc < ic; ++jc )
      {
        getBinaryInteractionCoefficient( pressure,
                                         temperature,
                                         componentProperties,
                                         salinity,
                                         ic,
                                         jc,
                                         kij,
                                         dkij_dT );
        stack.kij_data( ic, jc ) = kij;
        stack.kij_data( jc, ic ) = kij;
      }
    }
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 const & salinity,
                                arraySlice1d< real64 > const & logFugacityCoefficients )
{
  real64 compressibilityFactor = 0.0;

  // Step 1: Allocate the stack memory needed for the update
  StackVariables< false > stack( numComps );
  initialiseStack( numComps,
                   pressure,
                   temperature,
                   componentProperties,
                   salinity,
                   stack );

  // Step 2: Compute the mixture coefficients
  CubicModel::computeMixtureCoefficients( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  CubicModel::computeCompressibilityFactor( numComps,
                                            composition,
                                            stack,
                                            compressibilityFactor,
                                            nullptr /* No derivatives */ );

  // Step 4: Use mixture coefficients and compressibility factor to update fugacity coefficients
  CubicModel::computeLogFugacityCoefficients( numComps,
                                              composition,
                                              stack,
                                              compressibilityFactor,
                                              nullptr, /* No derivatives */
                                              logFugacityCoefficients,
                                              nullptr /* No derivatives */ );
}

template< typename EOS_TYPE >
template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficientsAndDerivs( integer const numComps,
                                         real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD1 > const & composition,
                                         ComponentProperties::KernelWrapper const & componentProperties,
                                         real64 const & salinity,
                                         arraySlice1d< real64 > const & logFugacityCoefficients,
                                         arraySlice2d< real64, USD2 > const & logFugacityCoefficientDerivs )
{
  integer constexpr numMaxDofs = StackVariables< true >::maxNumDof;
  integer const numDofs = 2 + numComps;
  // Allocate space for the compressibility derivatives
  real64 compressibilityFactor = 0.0;
  StackArray< real64, 1, numMaxDofs > compressibilityFactorDerivs( numDofs );

  // Step 1: Allocate the stack memory needed for the update
  StackVariables< true > stack( numComps );
  initialiseStack< true >( numComps,
                           pressure,
                           temperature,
                           componentProperties,
                           salinity,
                           stack );

  // Step 2: Compute the mixture coefficients
  CubicModel::template computeMixtureCoefficients< USD1, true >( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  CubicModel::template computeCompressibilityFactor< USD1, true >(
    numComps,
    composition,
    stack,
    compressibilityFactor,
    compressibilityFactorDerivs.toSlice() );

  // Step 4: Use mixture coefficients and compressibility factor to update fugacity coefficients
  CubicModel::template computeLogFugacityCoefficients< USD1, true, USD2 >(
    numComps,
    composition,
    stack,
    compressibilityFactor,
    compressibilityFactorDerivs.toSliceConst(),
    logFugacityCoefficients,
    logFugacityCoefficientDerivs );
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 const & salinity,
                              real64 & compressibilityFactor )
{
  // Step 1: Allocate the stack memory needed for the update
  StackVariables< false > stack( numComps );
  initialiseStack( numComps,
                   pressure,
                   temperature,
                   componentProperties,
                   salinity,
                   stack );

  // Step 2: Compute the mixture coefficients
  CubicModel::computeMixtureCoefficients( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  CubicModel::computeCompressibilityFactor( numComps,
                                            composition,
                                            stack,
                                            compressibilityFactor,
                                            nullptr /* No derivatives */ );
}

template< typename EOS_TYPE >
template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactorAndDerivs( integer const numComps,
                                       real64 const & pressure,
                                       real64 const & temperature,
                                       arraySlice1d< real64 const, USD1 > const & composition,
                                       ComponentProperties::KernelWrapper const & componentProperties,
                                       real64 const & salinity,
                                       real64 & compressibilityFactor,
                                       arraySlice1d< real64, USD2 > const & compressibilityFactorDerivs )
{
  // Step 1: Allocate the stack memory needed for the update
  StackVariables< true > stack( numComps );
  initialiseStack< true >( numComps,
                           pressure,
                           temperature,
                           componentProperties,
                           salinity,
                           stack );

  // Step 2: Compute the mixture coefficients
  CubicModel::template computeMixtureCoefficients< USD1, true >( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  CubicModel::template computeCompressibilityFactor< USD1, true >(
    numComps,
    composition,
    stack,
    compressibilityFactor,
    compressibilityFactorDerivs );
}

template< typename EOS_TYPE >
template< bool DERIVATIVES >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeWaterCoefficients( integer const h2o_index,
                          real64 const & pressure,
                          real64 const & temperature,
                          ComponentProperties::KernelWrapper const & componentProperties,
                          real64 const & salinity,
                          StackVariables< DERIVATIVES > & stack )
{
  real64 const Pc = componentProperties.m_componentCriticalPressure[h2o_index];
  real64 const Tc = componentProperties.m_componentCriticalTemperature[h2o_index];

  // Reduced properties
  real64 const Pr = pressure / Pc;
  real64 const Tr = temperature / Tc;

  // Inverse powers of Tr
  real64 const invTr2 = 1.0 / (Tr*Tr);
  real64 const invTr3 = invTr2 / Tr;
  real64 const invTr4 = invTr2 * invTr2;
  real64 const invTr5 = invTr2 * invTr3;

  // Salinity term
  real64 const m_s = power( salinity, 1.1 );

  // sqrtAlpha and alpha
  real64 const sqrtAlpha = 1.0 + 0.4530 * (1.0 - Tr * (1.0 - 0.0103 * m_s)) + 0.0034 * (invTr3 - 1.0);
  real64 const alpha = sqrtAlpha * sqrtAlpha;
  stack.aic[h2o_index] = EOS_TYPE::omegaA * Pr * invTr2 * alpha;

  if constexpr (DERIVATIVES)
  {
    // Derivatives of sqrtAlpha
    real64 const dsqrtAlpha_dT = -0.4530 * (1.0 - 0.0103 * m_s) / Tc - 0.0034 * 3.0 * invTr4 / Tc;
    real64 const d2sqrtAlpha_dT2 = 0.0034 * 12.0 * invTr5 / (Tc * Tc);

    // Derivatives of alpha
    real64 const dalpha_dT = 2.0 * sqrtAlpha * dsqrtAlpha_dT;
    real64 const d2alpha_dT2 = 2.0 * (dsqrtAlpha_dT * dsqrtAlpha_dT + sqrtAlpha * d2sqrtAlpha_dT2);

    // Derivatives of a
    stack.daic_dp[h2o_index] = EOS_TYPE::omegaA * invTr2 * alpha / Pc;
    stack.daic_dt[h2o_index] = EOS_TYPE::omegaA * Pr * (dalpha_dT * invTr2 - 2.0 * alpha * invTr3 / Tc);
    stack.d2aic_dt2[h2o_index] = EOS_TYPE::omegaA * Pr * ( d2alpha_dT2 * invTr2 - 4.0 * dalpha_dT * invTr3 / Tc + 6.0 * alpha * invTr4 / (Tc*Tc));
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
getBinaryInteractionCoefficient( real64 const & pressure,
                                 real64 const & temperature,
                                 ComponentProperties::KernelWrapper const & componentProperties,
                                 real64 const & salinity,
                                 integer const ic,
                                 integer const jc,
                                 real64 & kij,
                                 real64 & dkij_dT )
{
  GEOS_UNUSED_VAR( pressure );

  // Initialise default values
  kij = componentProperties.m_componentBinaryCoeff( ic, jc );
  dkij_dT = 0.0;

  integer i = ic;
  integer j = jc;
  integer type_i = componentProperties.m_componentType[i];
  integer type_j = componentProperties.m_componentType[j];

  if( isComponentType( type_i, ComponentType::Water ) &&
      !isComponentType( type_j, ComponentType::Water ) )
  {
    integer t = type_i;
    type_i = type_j;
    type_j = t;
    t = i;
    i = j;
    j = t;
  }

  if( isComponentType( type_i, ComponentType::Water ) ||
      !isComponentType( type_j, ComponentType::Water ))
  {
    return;
  }

  real64 const Tci = componentProperties.m_componentCriticalTemperature( i );
  real64 const dTr_dT = 1.0 / Tci;
  real64 const Tr = temperature * dTr_dT;

  // References:
  // - Soreide & Whitson (1992) https://doi.org/10.1016/0378-3812(92)85105-H
  // - Yan et al. (2011) https://doi.org/10.1016/j.ijggc.2011.08.004
  // - Chabab et al. (2019) https://doi.org/10.1016/j.ijggc.2019.102825
  // - Chabab et al. (2024) https://doi.org/10.1016/j.ijhydene.2023.10.290
  if( isComponentType( type_i, ComponentType::CarbonDioxide ))
  {
    // We have options here:
    // Original Soreide & Whitson (1992); Yan et al. (2011); Chabab et al. (2019)
    // Equation (14) Soreide & Whitson (1992)
    real64 constexpr A0 = -0.31092;
    real64 constexpr A1 = 0.23580;
    real64 constexpr A2 = -21.2566;
    real64 const a0 = 1.0 + 0.15587 * power( salinity, 0.7505 );
    real64 const a1 = 1.0 + 0.17837 * power( salinity, 0.979 );
    real64 constexpr e2 = -6.7222;
    real64 const a2 = LvArray::math::exp( e2*Tr - salinity );
    kij = A0*a0 + A1*a1*Tr + A2*a2;
    dkij_dT = (A1*a1 + e2*A2*a2) * dTr_dT;
  }
  else if( isComponentType( type_i, ComponentType::HydrogenSulphide ))
  {
    // Equation (15) Soreide & Whitson (1992)
    real64 constexpr A0 = -0.20441;
    real64 constexpr A1 = 0.23426;
    kij = A0 + A1*Tr;
    dkij_dT = A1*dTr_dT;
  }
  else if( isComponentType( type_i, ComponentType::Nitrogen ))
  {
    // Equation (13) Soreide & Whitson (1992)
    real64 constexpr A0 = -1.70235;
    real64 constexpr A1 = 0.44338;
    real64 const csw = power( salinity, 0.75 );
    real64 const a0 = 1.0 + 0.025587*csw;
    real64 const a1 = 1.0 + 0.08126*csw;
    kij = A0*a0 + A1*a1*Tr;
    dkij_dT = A1*a1*dTr_dT;
  }
  else if( isComponentType( type_i, ComponentType::Hydrogen ))
  {
    // Equation (12) and Table 5. Chabab et al. (2024)
    real64 constexpr D0 = -2.11917;
    real64 constexpr D1 = 0.14888;
    real64 constexpr D2 = -13.01835;
    real64 constexpr D3 = -0.43946;
    real64 constexpr a0 = -2.26322e-2;
    real64 constexpr a1 = -4.4736e-3;
    kij = D0*(1.0 + a0*salinity) + D1*Tr*(1.0 + a1*salinity) + D2*exp( D3*Tr );
    dkij_dT = dTr_dT*( D1*(1.0 + a1*salinity) + D2*D3*exp( D3*Tr ) );
  }
  else
  {
    // Hydrocarbon-water interaction
    real64 const omega = componentProperties.m_componentAcentricFactor( i );

    // Table 2 from Soreide & Whitson (1992)
    // See also equations (11) and (12)
    // For methane we will use the tuned parameters
    if( 0.02 < omega )
    {
      // Original parameters from Soreide & Whitson (1992)
      real64 const A0 = 1.1120 - 1.7369 * power( omega, -0.1 );
      real64 const A1 = 1.1001 + 0.8360 * omega;
      real64 const A2 = -0.15742 - 1.0988 * omega;

      real64 const a0 = 1.0 + 0.017407*salinity;
      real64 const a1 = 1.0 + 0.033516*salinity;
      real64 const a2 = 1.0 + 0.011478*salinity;
      kij = A0*a0 + A1*a1*Tr + A2*a2*Tr*Tr;
      real64 const dkij_dTr = A1*a1 + 2.0*A2*a2*Tr;
      dkij_dT = dkij_dTr * dTr_dT;
    }
    else
    {
      // Tuned parameters
      real64 const A0 = 1.616705 - 1.884534 * power( omega, -0.1 );
      real64 const A1 = 0.8157014 + 0.8723632 * omega;
      real64 const A2 = -0.0887821 - 0.8767864 * omega;

      real64 const a0 = 1.0 + 0.09988448*salinity;
      real64 const a1 = 1.0 + 0.1485516*salinity;
      real64 const a2 = 1.0 + 0.2111324*salinity;

      kij = A0*a0 + A1*a1*Tr + A2*a2*Tr*Tr;
      real64 const dkij_dTr = A1*a1 + 2.0*A2*a2*Tr;
      dkij_dT = dkij_dTr * dTr_dT;
    }
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_IMPL_HPP_
