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
 * @file CubicEOSPhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

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
  evaluate( real64 const & omega )
  {
    return ( omega < 0.49 )
      ? 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
      : 0.3796 + 1.485 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
  }

  static constexpr char const * catalogName(){ return "PengRobinson"; }
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
  evaluate( real64 const & omega )
  {
    return 0.480 + 1.574 * omega - 0.176 * omega * omega;
  }

  static constexpr char const * catalogName(){ return "SoaveRedlichKwong"; }
};

template< typename EOS_TYPE >
struct CubicEOSPhaseModel
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Generate a catalog name
   */
  static constexpr char const * catalogName(){ return EOS_TYPE::catalogName(); }

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
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arrayView1d< real64 const > const composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  arraySlice1d< real64 > const logFugacityCoefficients );

  /**
   * @brief Compute the mixture coefficients using pressure, temperature, composition and input
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] aPureCoefficient pure coefficient (A)
   * @param[out] bPureCoefficient pure coefficient (B)
   * @param[out] aMixtureCoefficient mixture coefficient (A)
   * @param[out] bMixtureCoefficient mixture coefficient (B)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arrayView1d< real64 const > const composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              arraySlice1d< real64 > const aPureCoefficient,
                              arraySlice1d< real64 > const bPureCoefficient,
                              real64 & aMixtureCoefficient,
                              real64 & bMixtureCoefficient );

  /**
   * @brief Compute the mixture coefficients derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] aMixtureCoefficientDerivs derivatives of mixture coefficient (A)
   * @param[out] bMixtureCoefficientDerivs derivatives of mixture coefficient (B)
   * @note Assumes that pressure and temperature are strictly positive
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arrayView1d< real64 const > const composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              arraySlice1d< real64 const > const aPureCoefficient,
                              arraySlice1d< real64 const > const bPureCoefficient,
                              real64 const aMixtureCoefficient,
                              real64 const bMixtureCoefficient,
                              arraySlice1d< real64 > const aMixtureCoefficientDerivs,
                              arraySlice1d< real64 > const bMixtureCoefficientDerivs );

  /**
   * @brief Compute the compressibility factor using compositions, BICs, and mixture coefficients
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] compressibilityFactor compressibility factor
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeCompressibilityFactor( integer const numComps,
                                arrayView1d< real64 const > const composition,
                                arrayView2d< real64 const > const & binaryInteractionCoefficients,
                                arraySlice1d< real64 const > const aPureCoefficient,
                                arraySlice1d< real64 const > const bPureCoefficient,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                real64 & compressibilityFactor );

  /**
   * @brief Compute the compressibility factor derivatives using mixture coefficients
   * @param[in] numComps number of components
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[in] compressibilityFactor the current compressibility factor
   * @param[in] aMixtureCoefficientDerivs derivatives of mixture coefficient (A)
   * @param[in] bMixtureCoefficientDerivs derivatives of mixture coefficient (B)
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   * @note Assumes that pressure and temperature are strictly positive
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                real64 const & compressibilityFactor,
                                arraySlice1d< real64 const > const aMixtureCoefficientDerivs,
                                arraySlice1d< real64 const > const bMixtureCoefficientDerivs,
                                arraySlice1d< real64 > const compressibilityFactorDerivs );

  /**
   * @brief Compute the log of the fugacity coefficients using compositions, BICs, compressibility factor and mixture coefficients
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[in] compressibilityFactor compressibility factor
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  arrayView1d< real64 const > const composition,
                                  arrayView2d< real64 const > const & binaryInteractionCoefficients,
                                  real64 const & compressibilityFactor,
                                  arraySlice1d< real64 const > const aPureCoefficient,
                                  arraySlice1d< real64 const > const bPureCoefficient,
                                  real64 const & aMixtureCoefficient,
                                  real64 const & bMixtureCoefficient,
                                  arraySlice1d< real64 > const logFugacityCoefficients );

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
  GEOS_FORCE_INLINE
  static void
    solveCubicPolynomial( real64 const & m3,
                          real64 const & m2,
                          real64 const & m1,
                          real64 const & m0,
                          real64 ( &roots )[3],
                          integer & numRoots );

};

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arrayView1d< real64 const > const composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                arraySlice1d< real64 > const logFugacityCoefficients )
{
  // step 0: allocate the stack memory needed for the update
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > aPureCoefficient( numComps );
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > bPureCoefficient( numComps );
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  real64 compressibilityFactor = 0.0;

  arrayView2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // step 1: compute the mixture coefficients aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // step 2: use mixture coefficients to update the compressibility factor
  computeCompressibilityFactor( numComps, // number of components
                                composition, // cell input
                                binaryInteractionCoefficients, // user input
                                aPureCoefficient, // computed by computeMixtureCoefficients
                                bPureCoefficient,
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                compressibilityFactor ); // output

  // step 3: use mixture coefficients and compressibility factor to update fugacity coefficients
  computeLogFugacityCoefficients( numComps, // number of components
                                  composition, // cell input
                                  binaryInteractionCoefficients, // user input
                                  compressibilityFactor, // computed by computeCompressibilityFactor
                                  aPureCoefficient, // computed by computeMixtureCoefficients
                                  bPureCoefficient,
                                  aMixtureCoefficient,
                                  bMixtureCoefficient,
                                  logFugacityCoefficients ); // output
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arrayView1d< real64 const > const composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            arraySlice1d< real64 > const aPureCoefficient,
                            arraySlice1d< real64 > const bPureCoefficient,
                            real64 & aMixtureCoefficient,
                            real64 & bMixtureCoefficient )
{
  arrayView1d< real64 const > const & criticalPressure = componentProperties.m_componentCriticalPressure;
  arrayView1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
  arrayView1d< real64 const > const & acentricFactor = componentProperties.m_componentAcentricFactor;
  arrayView2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // mixture coefficients
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const m = EOS_TYPE::evaluate( acentricFactor[ic] );
    real64 const pr = pressure / criticalPressure[ic];
    real64 const tr = temperature / criticalTemperature[ic];
    aPureCoefficient[ic] = EOS_TYPE::omegaA * pr / (tr*tr) * pow( 1.0 + m * ( 1.0 - sqrt( tr ) ), 2.0 );
    bPureCoefficient[ic] = EOS_TYPE::omegaB * pr / tr;
  }

  aMixtureCoefficient = 0.0;
  bMixtureCoefficient = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      aMixtureCoefficient += ( composition[ic] * composition[jc] * ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] ) );
    }
    bMixtureCoefficient += composition[ic] * bPureCoefficient[ic];
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arrayView1d< real64 const > const composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            arraySlice1d< real64 const > const aPureCoefficient,
                            arraySlice1d< real64 const > const bPureCoefficient,
                            real64 const aMixtureCoefficient,
                            real64 const bMixtureCoefficient,
                            arraySlice1d< real64 > const aMixtureCoefficientDerivs,
                            arraySlice1d< real64 > const bMixtureCoefficientDerivs )
{
  arrayView1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
  arrayView1d< real64 const > const & acentricFactor = componentProperties.m_componentAcentricFactor;
  arrayView2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > daPureCoefficient_dx( numComps );

  // Calculate pressure derivatives
  aMixtureCoefficientDerivs[Deriv::dP] = aMixtureCoefficient / pressure;
  bMixtureCoefficientDerivs[Deriv::dP] = bMixtureCoefficient / pressure;

  // Calculate temperature derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const m = EOS_TYPE::evaluate( acentricFactor[ic] );
    real64 const sqrtTr = sqrt( temperature / criticalTemperature[ic] );
    real64 const mt = 1.0 + m * (1.0 - sqrtTr);
    daPureCoefficient_dx[ic] = -aPureCoefficient[ic] * (2.0/temperature + m/(mt*sqrtTr*criticalTemperature[ic]));
  }
  aMixtureCoefficientDerivs[Deriv::dT] = 0.0;
  bMixtureCoefficientDerivs[Deriv::dT] = -bMixtureCoefficient / temperature;

  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const coeff = composition[ic] * composition[jc] * ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) / sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      aMixtureCoefficientDerivs[Deriv::dT] += 0.5 * coeff * (daPureCoefficient_dx[ic]*aPureCoefficient[jc] + daPureCoefficient_dx[jc]*aPureCoefficient[ic]);
    }
  }

  // Calculate composition derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    aMixtureCoefficientDerivs[Deriv::dC+ic] = 0.0;
    bMixtureCoefficientDerivs[Deriv::dC+ic] = 0.0;
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const coeff = ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      aMixtureCoefficientDerivs[Deriv::dC+ic] += coeff * composition[jc];
      aMixtureCoefficientDerivs[Deriv::dC+jc] += coeff * composition[ic];
    }
    bMixtureCoefficientDerivs[Deriv::dC+ic] = bPureCoefficient[ic];
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              arrayView1d< real64 const > const composition,
                              arrayView2d< real64 const > const & binaryInteractionCoefficients,
                              arraySlice1d< real64 const > const aPureCoefficient,
                              arraySlice1d< real64 const > const bPureCoefficient,
                              real64 const & aMixtureCoefficient,
                              real64 const & bMixtureCoefficient,
                              real64 & compressibilityFactor )
{
  // a Z3 + b Z2 + cZ + d = 0
  real64 const a = 1.0;
  real64 const b = ( EOS_TYPE::delta1 + EOS_TYPE::delta2 - 1.0 ) * bMixtureCoefficient - 1.0;
  real64 const c = aMixtureCoefficient + EOS_TYPE::delta1 * EOS_TYPE::delta2 * bMixtureCoefficient * bMixtureCoefficient
                   - ( EOS_TYPE::delta1 + EOS_TYPE::delta2 ) * bMixtureCoefficient * ( bMixtureCoefficient + 1.0 );
  real64 const d = -( aMixtureCoefficient * bMixtureCoefficient
                      + EOS_TYPE::delta1 * EOS_TYPE::delta2 * bMixtureCoefficient * bMixtureCoefficient * ( bMixtureCoefficient + 1.0 ) );

  real64 roots[3]{};
  integer numRoots = 0;
  solveCubicPolynomial( a, b, c, d, roots, numRoots );

  if( numRoots == 1 )
  {
    compressibilityFactor = roots[0];
  }
  else
  {
    real64 zMin =  LvArray::NumericLimits< real64 >::max;
    real64 zMax = -LvArray::NumericLimits< real64 >::max;

    for( integer i = 0; i < 3; ++i )
    {
      // skip unphysical roots
      if( roots[i] > bMixtureCoefficient )
      {
        // choose the root according to Gibbs' free energy minimization
        if( zMin > roots[i] )
        {
          zMin = roots[i];
        }
        if( zMax < roots[i] )
        {
          zMax = roots[i];
        }
      }
    }

    stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > logFugacityCoefficientsMax( numComps );
    stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > logFugacityCoefficientsMin( numComps );
    computeLogFugacityCoefficients( numComps, composition, binaryInteractionCoefficients, zMin,
                                    aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient,
                                    logFugacityCoefficientsMin.toSlice() );
    computeLogFugacityCoefficients( numComps, composition, binaryInteractionCoefficients, zMax,
                                    aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient,
                                    logFugacityCoefficientsMax.toSlice() );

    real64 dG = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dG += composition[ic] * ( logFugacityCoefficientsMin[ic] - logFugacityCoefficientsMax[ic] );
    }
    compressibilityFactor = ( dG < 0 ) ? zMin : zMax;
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & aMixtureCoefficient,
                              real64 const & bMixtureCoefficient,
                              real64 const & compressibilityFactor,
                              arraySlice1d< real64 const > const aMixtureCoefficientDerivs,
                              arraySlice1d< real64 const > const bMixtureCoefficientDerivs,
                              arraySlice1d< real64 > const compressibilityFactorDerivs )
{
  constexpr real64 d1pd2 = EOS_TYPE::delta1 + EOS_TYPE::delta2;
  constexpr real64 d1xd2 = EOS_TYPE::delta1 * EOS_TYPE::delta2;

  constexpr real64 a = 1.0;
  real64 const b = ( d1pd2 - 1.0 ) * bMixtureCoefficient - 1.0;
  real64 const c = aMixtureCoefficient + d1xd2 * bMixtureCoefficient * bMixtureCoefficient
                   - d1pd2 * bMixtureCoefficient * ( bMixtureCoefficient + 1.0 );

  // Implicit differentiation scale
  real64 const denominator = (3.0*a*compressibilityFactor + 2.0*b)*compressibilityFactor + c;
  constexpr real64 epsilon = MultiFluidConstants::epsilon;
  real64 const scalingFactor = fabs( denominator ) < epsilon ? 0.0 : -1.0 / denominator;

  // Given derivative of the mixture parameters a and b w.r.t. variable X, calculate the derivative of the
  // compressibility factor (z-factor) w.r.t. X
  auto const calculateDerivative = [&]( real64 const da_dX, real64 const db_dX ) -> real64 {
    // a Z3 + b Z2 + cZ + d = 0
    // Derivatives for a,b,c,d
    real64 const dbdx = ( d1pd2 - 1.0 ) * db_dX;
    real64 const dcdx = da_dX + ( 2.0*(d1xd2-d1pd2) * bMixtureCoefficient - d1pd2 )*db_dX;
    real64 const dddx = -(aMixtureCoefficient*db_dX + da_dX*bMixtureCoefficient
                          + d1xd2*((3.0*bMixtureCoefficient+2.0)*bMixtureCoefficient*db_dX));
    return (((dbdx*compressibilityFactor) + dcdx)*compressibilityFactor + dddx) * scalingFactor;
  };

  // Pressure and temperature derivatives
  for( integer const kc : {Deriv::dP, Deriv::dT} )
  {
    compressibilityFactorDerivs[kc] = calculateDerivative( aMixtureCoefficientDerivs[kc], bMixtureCoefficientDerivs[kc] );
  }

  // Composition derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const kc = Deriv::dC+ic;
    compressibilityFactorDerivs[kc] = calculateDerivative( aMixtureCoefficientDerivs[kc], bMixtureCoefficientDerivs[kc] );
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                arrayView1d< real64 const > const composition,
                                arrayView2d< real64 const > const & binaryInteractionCoefficients,
                                real64 const & compressibilityFactor,
                                arraySlice1d< real64 const > const aPureCoefficient,
                                arraySlice1d< real64 const > const bPureCoefficient,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                arraySlice1d< real64 > const logFugacityCoefficients )
{
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > ki( numComps );

  // ki
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      ki[ic] += composition[jc] * ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
    }
  }

  // E
  real64 const E = log( ( compressibilityFactor + EOS_TYPE::delta1 * bMixtureCoefficient )
                        / ( compressibilityFactor + EOS_TYPE::delta2 * bMixtureCoefficient ) );
  real64 const F = log( compressibilityFactor - bMixtureCoefficient );
  real64 const G = 1.0 / ( ( EOS_TYPE::delta1 - EOS_TYPE::delta2 ) * bMixtureCoefficient );
  real64 const A = aMixtureCoefficient;

  // ln phi
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const B = bPureCoefficient[ic] / bMixtureCoefficient;
    logFugacityCoefficients[ic] = ( compressibilityFactor - 1 ) * B - F - G * ( 2 * ki[ic] - A * B ) * E;
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
    roots[1] = -2 * qSqrt * cos( ( theta + 2 * MultiFluidConstants::pi ) / 3 ) - a1 / 3;
    roots[2] = -2 * qSqrt * cos( ( theta + 4 * MultiFluidConstants::pi ) / 3 ) - a1 / 3;
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

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
