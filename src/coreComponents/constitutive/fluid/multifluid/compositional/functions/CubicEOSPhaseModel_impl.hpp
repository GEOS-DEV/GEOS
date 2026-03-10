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
 * @file CubicEOSPhaseModel_impl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_IMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_IMPL_HPP_

#include "CubicEOSPhaseModel.hpp"
#include "common/logger/Logger.hpp"

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
CubicEOSPhaseModel< EOS_TYPE >::
initialiseStack( integer const numComps,
                 real64 const & pressure,
                 real64 const & temperature,
                 ComponentProperties::KernelWrapper const & componentProperties,
                 StackVariables< DERIVATIVES > & stack )
{
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficients< DERIVATIVES >( ic,
                                            pressure,
                                            temperature,
                                            componentProperties,
                                            stack );
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                arraySlice1d< real64 > const & logFugacityCoefficients )
{
  real64 compressibilityFactor = 0.0;

  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSlice();

  // Step 1: Allocate the stack memory needed for the update
  StackVariables< false > stack( numComps, binaryInteractionCoefficients );
  initialiseStack( numComps,
                   pressure,
                   temperature,
                   componentProperties,
                   stack );

  // Step 2: Compute the mixture coefficients
  computeMixtureCoefficients( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  computeCompressibilityFactor( numComps,
                                composition,
                                stack,
                                compressibilityFactor,
                                nullptr /* No derivatives */ );

  // Step 4: Use mixture coefficients and compressibility factor to update fugacity coefficients
  computeLogFugacityCoefficients( numComps,
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
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficientsAndDerivs( integer const numComps,
                                         real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD1 > const & composition,
                                         ComponentProperties::KernelWrapper const & componentProperties,
                                         arraySlice1d< real64 > const & logFugacityCoefficients,
                                         arraySlice2d< real64, USD2 > const & logFugacityCoefficientDerivs )
{
  integer constexpr numMaxDofs = StackVariables< true >::maxNumDof;
  integer const numDofs = 2 + numComps;
  // Allocate space for the compressibility derivatives
  real64 compressibilityFactor = 0.0;
  StackArray< real64, 1, numMaxDofs > compressibilityFactorDerivs( numDofs );

  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSlice();

  // Step 1: Allocate the stack memory needed for the update
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  StackVariables< true > stack( numComps, binaryInteractionCoefficients, derivs.toSliceConst() );
  initialiseStack< true >( numComps,
                           pressure,
                           temperature,
                           componentProperties,
                           stack );

  // Step 2: Compute the mixture coefficients
  computeMixtureCoefficients< USD1, true >( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  computeCompressibilityFactor< USD1, true >( numComps,
                                              composition,
                                              stack,
                                              compressibilityFactor,
                                              compressibilityFactorDerivs.toSlice() );

  // Step 4: Use mixture coefficients and compressibility factor to update fugacity coefficients
  computeLogFugacityCoefficients< USD1, true >( numComps,
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
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 & compressibilityFactor )
{
  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSlice();

  // Step 1: Allocate the stack memory needed for the update
  StackVariables< false > stack( numComps, binaryInteractionCoefficients );
  initialiseStack( numComps,
                   pressure,
                   temperature,
                   componentProperties,
                   stack );

  // Step 2: Compute the mixture coefficients
  computeMixtureCoefficients( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  computeCompressibilityFactor( numComps,
                                composition,
                                stack,
                                compressibilityFactor,
                                nullptr /* No derivatives */ );
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactorAndDerivs( integer const numComps,
                                       real64 const & pressure,
                                       real64 const & temperature,
                                       arraySlice1d< real64 const, USD > const & composition,
                                       ComponentProperties::KernelWrapper const & componentProperties,
                                       real64 & compressibilityFactor,
                                       arraySlice1d< real64 > const & compressibilityFactorDerivs )
{
  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSlice();

  // Step 1: Allocate the stack memory needed for the update
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  StackVariables< true > stack( numComps, binaryInteractionCoefficients, derivs );
  initialiseStack< true >( numComps,
                           pressure,
                           temperature,
                           componentProperties,
                           stack );

  // Step 2: Compute the mixture coefficients
  computeMixtureCoefficients< USD, true >( numComps, composition, stack );

  // Step 3: Compute the compressibility factor (Z)
  computeCompressibilityFactor< USD, true >( numComps,
                                             composition,
                                             stack,
                                             compressibilityFactor,
                                             compressibilityFactorDerivs );
}

template< typename EOS_TYPE >
void
CubicEOSPhaseModel< EOS_TYPE >::
calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                 arraySlice1d< real64 > const & dimensionalVolumeShift )
{
  integer const numComps = componentProperties.getNumberOfComponents();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const Vs = componentProperties.getComponentVolumeShift()[ic];
    real64 const Pc = componentProperties.getComponentCriticalPressure()[ic];
    real64 const Tc = componentProperties.getComponentCriticalTemperature()[ic];
    real64 constexpr omegaB = EOS_TYPE::omegaB;
    dimensionalVolumeShift[ic] = constants::gasConstant * Vs * omegaB * Tc / Pc;
  }
}

template< typename EOS_TYPE >
template< bool DERIVATIVES >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         StackVariables< DERIVATIVES > & stack )
{
  real64 const Pc = componentProperties.m_componentCriticalPressure[ic];
  real64 const Tc = componentProperties.m_componentCriticalTemperature[ic];
  real64 const omega = componentProperties.m_componentAcentricFactor[ic];

  // Reduced properties
  real64 const Pr = pressure / Pc;
  real64 const Tr = temperature / Tc;
  real64 const sqrtTr = LvArray::math::sqrt( Tr );

  // Kappa and alpha
  real64 const kappa = EOS_TYPE::kappa( omega );
  real64 const sqrtAlpha = 1.0 + kappa * (1.0 - sqrtTr);
  real64 const alpha = sqrtAlpha*sqrtAlpha;

  // Values
  stack.aic[ic] = EOS_TYPE::omegaA * Pr * alpha / (Tr*Tr);
  stack.bic[ic] = EOS_TYPE::omegaB * Pr  / Tr;

  if constexpr (DERIVATIVES)
  {
    // Derivatives w.r.t pressure
    stack.daic_dp[ic] = EOS_TYPE::omegaA * alpha / (Pc * Tr * Tr);
    stack.dbic_dp[ic] = EOS_TYPE::omegaB / (Pc * Tr);

    // Derivatives w.r.t temperature
    // Derivatives of sqrtAlpha
    real64 const dsqrtAlpha_dT = -0.5 * kappa / (Tc * sqrtTr);
    real64 const d2sqrtAlpha_dT2 = -0.5 * dsqrtAlpha_dT / (Tr*Tc);

    // Derivatives of alpha
    real64 const dalpha_dT = 2.0 * dsqrtAlpha_dT * sqrtAlpha;
    real64 const d2alpha_dT2 = 2.0 * (dsqrtAlpha_dT * dsqrtAlpha_dT + sqrtAlpha * d2sqrtAlpha_dT2);

    real64 const da_dT = EOS_TYPE::omegaA * Pr * (dalpha_dT / (Tr * Tr) - 2.0 * alpha / (Tc * Tr * Tr * Tr));
    real64 const d2a_dT2 = EOS_TYPE::omegaA * Pr / (Tr * Tr) * (d2alpha_dT2 - 4.0 * dalpha_dT / (Tc * Tr ) + 6.0 * alpha / (Tc * Tc * Tr * Tr));

    stack.daic_dt[ic] = da_dT;
    stack.dbic_dt[ic] = -EOS_TYPE::omegaB * Pr / (Tc * Tr * Tr);

    // Second derivatives w.r.t temperature
    stack.d2aic_dt2[ic] = d2a_dT2;
    stack.d2bic_dt2[ic] =  2.0 * EOS_TYPE::omegaB * Pr / (Tc * Tc * Tr * Tr * Tr);
  }
}

template< typename EOS_TYPE >
template< integer USD, bool DERIVATIVES >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            arraySlice1d< real64 const, USD > const & composition,
                            StackVariables< DERIVATIVES > & stack )
{
  // Binary interaction coefficients
  arraySlice2d< real64 const > const & kij = stack.kij;
  stack.aMixture = 0.0;
  stack.bMixture = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      stack.aMixture += composition[ic] * composition[jc] * ( 1.0 - kij( ic, jc ) ) * LvArray::math::sqrt( stack.aic[ic] * stack.aic[jc] );
    }
    stack.bMixture += composition[ic] * stack.bic[ic];
  }
  if constexpr (DERIVATIVES)
  {
    LvArray::forValuesInSlice( stack.daMixture, setZero );
    LvArray::forValuesInSlice( stack.dbMixture, setZero );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      for( integer jc = 0; jc < numComps; ++jc )
      {
        real64 const sqrt_aiaj = LvArray::math::sqrt( stack.aic[ic] * stack.aic[jc] );
        real64 const kij_term = 1.0 - kij( ic, jc );
        real64 const aij = kij_term * sqrt_aiaj;
        real64 const coeff = 0.5 * kij_term / sqrt_aiaj;

        real64 const daij_dT = coeff * (stack.aic[jc] * stack.daic_dt[ic] + stack.aic[ic] * stack.daic_dt[jc]);

        real64 const daij_dP = coeff * (stack.aic[jc] * stack.daic_dp[ic] + stack.aic[ic] * stack.daic_dp[jc]);

        stack.daMixture[Deriv::dP] += composition[ic] * composition[jc] * daij_dP;
        stack.daMixture[Deriv::dT] += composition[ic] * composition[jc] * daij_dT;
        stack.daMixture[Deriv::dC+ic] += composition[jc] * aij;
        stack.daMixture[Deriv::dC+jc] += composition[ic] * aij;
      }
      stack.dbMixture[Deriv::dP] += composition[ic] * stack.dbic_dp[ic];
      stack.dbMixture[Deriv::dT] += composition[ic] * stack.dbic_dt[ic];
      stack.dbMixture[Deriv::dC+ic] = stack.bic[ic];
    }
    if( 0 < stack.dkij_dT.size( 0 ))
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        for( integer jc = 0; jc < numComps; ++jc )
        {
          real64 const sqrt_aiaj = LvArray::math::sqrt( stack.aic[ic] * stack.aic[jc] );
          real64 const dkij_term_dT = -stack.dkij_dT( ic, jc );
          stack.daMixture[Deriv::dT] += composition[ic] * composition[jc] * dkij_term_dT * sqrt_aiaj;
        }
      }
    }
  }
}

template< typename EOS_TYPE >
template< integer USD, bool DERIVATIVES >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              arraySlice1d< real64 const, USD > const & composition,
                              StackVariables< DERIVATIVES > const & stack,
                              real64 & compressibilityFactor,
                              StackDerivativeType< 1, DERIVATIVES > const & compressibilityFactorDerivs,
                              SelectedRoot const selectedRoot )
{
  // a Z^3 + b Z^2 + c Z + d = 0
  real64 const A = stack.aMixture;
  real64 const B = stack.bMixture;
  real64 constexpr d1pd2 = EOS_TYPE::delta1 + EOS_TYPE::delta2;
  real64 constexpr d1xd2 = EOS_TYPE::delta1 * EOS_TYPE::delta2;

  real64 const a = 1.0;
  real64 const b = (d1pd2 - 1.0) * B - 1.0;
  real64 const c = A + d1xd2 * B * B - d1pd2 * B * (B + 1.0);
  real64 const d = -(A * B + d1xd2 * B * B * (B + 1.0));

  real64 roots[3]{};
  integer numRoots = 0;
  solveCubicPolynomial( a, b, c, d, roots, numRoots );

  if( numRoots == 1 )
  {
    compressibilityFactor = roots[0];
  }
  else
  {
    // Choose
    real64 zMin = LvArray::NumericLimits< real64 >::max;
    real64 zMax = -LvArray::NumericLimits< real64 >::max;
    for( integer r = 0; r < numRoots; r++ )
    {
      if( roots[r] > B )
      {
        zMin = LvArray::math::min( zMin, roots[r] );
        zMax = LvArray::math::max( zMax, roots[r] );
      }
    }
    if( selectedRoot == SelectedRoot::MINIMUM )
    {
      compressibilityFactor = zMin;
    }
    else if( selectedRoot == SelectedRoot::MAXIMUM )
    {
      compressibilityFactor = zMax;
    }
    else
    {
      integer constexpr maxNumComp = StackVariables< true >::maxNumComp;
      real64 minEnergy = LvArray::NumericLimits< real64 >::max;
      StackArray< real64, 1, maxNumComp > logFugacityCoefficients( numComps );
      for( real64 const z : {zMin, zMax} )
      {
        computeLogFugacityCoefficients< USD, false, 0 >( numComps,
                                                         composition,
                                                         stack,
                                                         z,
                                                         nullptr,
                                                         logFugacityCoefficients.toSlice(),
                                                         nullptr );
        real64 dG = 0.0;
        for( integer ic = 0; ic < numComps; ++ic )
        {
          dG += composition[ic] * logFugacityCoefficients[ic];
        }
        if( dG < minEnergy )
        {
          minEnergy = dG;
          compressibilityFactor = z;
        }
      }
    }
  }
  if constexpr (DERIVATIVES)
  {
    real64 const Z = compressibilityFactor;
    // Implicit differentiation scale
    real64 const denominator = (3.0*a*Z + 2.0*b)*Z + c;
    real64 const scalingFactor = LvArray::math::abs( denominator ) < MultiFluidConstants::epsilon ? 0.0 : -1.0 / denominator;

    integer const numDofs = numComps + 2;

    for( integer idof = 0; idof < numDofs; ++idof )
    {
      // Given derivative of the mixture parameters a and b w.r.t. variable X, calculate the derivative of the
      // compressibility factor (z-factor) w.r.t. X
      real64 const da_dX = stack.daMixture[idof];
      real64 const db_dX = stack.dbMixture[idof];
      // a Z3 + b Z2 + cZ + d = 0
      // Derivatives for a,b,c,d
      real64 const dbdx = ( d1pd2 - 1.0 ) * db_dX;
      real64 const dcdx = da_dX + ( 2.0*(d1xd2-d1pd2) * B - d1pd2 )*db_dX;
      real64 const dddx = -(A*db_dX + da_dX*B + d1xd2*((3.0*B+2.0)*B*db_dX));
      compressibilityFactorDerivs[idof] = ((dbdx*Z + dcdx)*Z + dddx) * scalingFactor;
    }
  }
  else
  {
    GEOS_UNUSED_VAR( compressibilityFactorDerivs );
  }
}

template< typename EOS_TYPE >
template< integer USD1, bool DERIVATIVES, integer USD2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                arraySlice1d< real64 const, USD1 > const & composition,
                                StackVariables< DERIVATIVES > const & stack,
                                real64 const & compressibilityFactor,
                                StackConstDerivativeType< 1, DERIVATIVES > const & compressibilityFactorDerivs,
                                arraySlice1d< real64 > const & logFugacityCoefficients,
                                StackDerivativeType< 2, DERIVATIVES, USD2 > const & logFugacityCoefficientDerivs )
{
  constexpr integer maxNumComp = StackVariables< DERIVATIVES >::maxNumComp;
  StackArray< real64, 1, maxNumComp > ki( numComps );

  // ki
  for( integer ic = 0; ic < numComps; ++ic )
  {
    ki[ic] = 0.0;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      ki[ic] += composition[jc] * ( 1.0 - stack.kij( ic, jc ) ) * sqrt( stack.aic[ic] * stack.aic[jc] );
    }
  }

  real64 const A = stack.aMixture;
  real64 const B = stack.bMixture;
  real64 const Z = compressibilityFactor;
  // E
  real64 const expE = ( Z + EOS_TYPE::delta1 * B ) / ( Z + EOS_TYPE::delta2 * B );
  real64 const expF = Z - B;
  GEOS_ERROR_IF( expE < MultiFluidConstants::epsilon || expF < MultiFluidConstants::epsilon,
                 GEOS_FMT( "Cubic EOS failed with exp(E)={} and exp(F)={}", expE, expF ));
  real64 const E = log( expE );
  real64 const F = log( expF );
  real64 const G = 1.0 / ( ( EOS_TYPE::delta1 - EOS_TYPE::delta2 ) * B );

  // ln phi
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const Bi = stack.bic[ic] / B;
    logFugacityCoefficients[ic] = ( Z - 1 ) * Bi - F - G * ( 2 * ki[ic] - A * Bi ) * E;
  }

  if constexpr (DERIVATIVES)
  {
    constexpr integer maxNumDof = StackVariables< DERIVATIVES >::maxNumDof;
    integer const numDofs = numComps + 2;
    StackArray< real64, 2, maxNumComp * maxNumDof > dki( numComps, numDofs );

    // dki
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dki( ic, Deriv::dP ) = 0.0;
      dki( ic, Deriv::dT ) = 0.0;

      real64 const sqrtAic_ic = LvArray::math::sqrt( stack.aic[ic] );
      for( integer jc = 0; jc < numComps; ++jc )
      {
        real64 const sqrtAic_jc = LvArray::math::sqrt( stack.aic[jc] );
        real64 const kij = stack.kij( ic, jc );

        // Derivative with respect to pressure
        real64 const dSqrt_dP = 0.5 * (sqrtAic_jc / sqrtAic_ic * stack.daic_dp[ic] + sqrtAic_ic / sqrtAic_jc * stack.daic_dp[jc]);
        dki( ic, Deriv::dP ) += composition[jc] * (1.0 - kij) * dSqrt_dP;

        // Derivative with respect to temperature
        real64 const dSqrt_dT = 0.5 * (sqrtAic_jc / sqrtAic_ic * stack.daic_dt[ic] + sqrtAic_ic / sqrtAic_jc * stack.daic_dt[jc]);
        dki( ic, Deriv::dT ) += composition[jc] * (1.0 - kij) * dSqrt_dT;

        // Derivative with respect to composition
        dki( ic, Deriv::dC + jc ) = (1.0 - kij) * sqrtAic_ic * sqrtAic_jc;
      }
    }
    if( 0 < stack.dkij_dT.size() )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        for( integer jc = 0; jc < numComps; ++jc )
        {
          real64 const sqrt_aiaj = LvArray::math::sqrt( stack.aic[ic] * stack.aic[jc] );
          real64 const dkij_term_dT = -stack.dkij_dT( ic, jc );
          dki( ic, Deriv::dT ) += composition[jc] * dkij_term_dT * sqrt_aiaj;
        }
      }
    }

    auto const & dZ = compressibilityFactorDerivs;
    auto const & dA = stack.daMixture;
    auto const & dB = stack.dbMixture;

    auto const calculateDerivatives = [&]( integer const idof ){
      real64 const dE_dX =  (dZ[idof] + EOS_TYPE::delta1*dB[idof])/( Z + EOS_TYPE::delta1 * B )
                           -(dZ[idof] + EOS_TYPE::delta2*dB[idof])/( Z + EOS_TYPE::delta2 * B );

      //real64 const F = log( compressibilityFactor - bMixtureCoefficient );
      real64 const dF_dX = (dZ[idof] - dB[idof])/(Z - B);

      // real64 const G = 1.0 / ( ( EOS_TYPE::delta1 - EOS_TYPE::delta2 ) * B );
      real64 const dG_dX = -G * dB[idof] / B;

      real64 const dA_dX = dA[idof];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        real64 const Bi = stack.bic[ic] / B;
        real64 dBi_dX = -Bi*dB[idof] / B;
        if( idof == Deriv::dP )
        {
          dBi_dX += stack.dbic_dp[ic] / B;
        }
        else if( idof == Deriv::dT )
        {
          dBi_dX += stack.dbic_dt[ic] / B;
        }

        // lnPhi = ( Z - 1 ) * Bi - F - G * ( 2 * ki[ic] - A * Bi ) * E;
        logFugacityCoefficientDerivs( ic, idof ) =
          dZ[idof]*Bi + ( Z - 1 ) * dBi_dX
          - dF_dX
          - dG_dX * ( 2 * ki[ic] - A * Bi ) * E
          - G * ( 2 * dki( ic, idof ) - dA_dX * Bi - A * dBi_dX ) * E
          - G * ( 2 * ki[ic] - A * Bi ) * dE_dX;
      }
    };

    calculateDerivatives( Deriv::dP );
    calculateDerivatives( Deriv::dT );
    for( integer jc = 0; jc < numComps; ++jc )
    {
      calculateDerivatives( Deriv::dC+jc );
    }
  }
  else
  {
    GEOS_UNUSED_VAR( compressibilityFactorDerivs );
    GEOS_UNUSED_VAR( logFugacityCoefficientDerivs );
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
solveCubicPolynomial( real64 const & m3,
                      real64 const & m2,
                      real64 const & m1,
                      real64 const & m0,
                      real64 (& roots)[3],
                      integer & numRoots )
{
  // cubic equation : m3 * x^3 + m2 * x^2 + m1 *x + m0  = 0
  real64 const a1 = m2 / m3;
  real64 const a2 = m1 / m3;
  real64 const a3 = m0 / m3;

  real64 const q = ( a1 * a1 - 3 * a2 ) / 9;
  real64 const r = ( 2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3 ) / 54;
  real64 const qCubed = q * q * q;
  real64 const d = qCubed - r * r;

  // three real roots
  if( d >= 0 )
  {
    real64 const theta = acos( r / sqrt( qCubed ) );
    real64 const qSqrt = sqrt( q );
    roots[0] = -2 * qSqrt * cos( theta / 3 ) - a1 / 3;
    roots[1] = -2 * qSqrt * cos( ( theta + 2 * constants::pi ) / 3 ) - a1 / 3;
    roots[2] = -2 * qSqrt * cos( ( theta + 4 * constants::pi ) / 3 ) - a1 / 3;
    numRoots = 3;
  }
  // one real root
  else
  {
    real64 e = pow( sqrt( -d ) + LvArray::math::abs( r ), 1. / 3. );
    if( r > 0 )
    {
      e = -e;
    }
    roots[0] = ( e + q / e ) - a1 / 3.;
    numRoots = 1;
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_IMPL_HPP_
