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

namespace geos
{

namespace constitutive
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
};

template< typename EOS_TYPE >
struct CubicEOSPhaseModel
{
public:
  /**
   * @brief Main entry point of the cubic EOS model
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] criticalPressure critical pressures
   * @param[in] criticalTemperature critical temperatures
   * @param[in] acentricFactor acentric factors
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  compute( integer const numComps,
           real64 const & pressure,
           real64 const & temperature,
           arrayView1d< real64 const > const composition,
           arrayView1d< real64 const > const criticalPressure,
           arrayView1d< real64 const > const criticalTemperature,
           arrayView1d< real64 const > const acentricFactor,
           real64 const & binaryInteractionCoefficients,
           arraySlice1d< real64 > const logFugacityCoefficients );

  /**
   * @brief Compute the mixture coefficients using pressure, temperature, composition and input
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] criticalPressure critical pressures
   * @param[in] criticalTemperature critical temperatures
   * @param[in] acentricFactor acentric factors
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
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
                              arrayView1d< real64 const > const criticalPressure,
                              arrayView1d< real64 const > const criticalTemperature,
                              arrayView1d< real64 const > const acentricFactor,
                              real64 const & binaryInteractionCoefficients,
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
   * @param[in] criticalPressure critical pressures
   * @param[in] criticalTemperature critical temperatures
   * @param[in] acentricFactor acentric factors
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] daMixtureCoefficient_dp derivative of mixture coefficient (A) wrt pressure
   * @param[out] dbMixtureCoefficient_dp derivative of mixture coefficient (B) wrt pressure
   * @param[out] daMixtureCoefficient_dt derivative of mixture coefficient (A) wrt temperature
   * @param[out] dbMixtureCoefficient_dt derivative of mixture coefficient (B) wrt temperature
   * @param[out] daMixtureCoefficient_dz derivative of mixture coefficient (A) wrt composition
   * @param[out] dbMixtureCoefficient_dz derivative of mixture coefficient (B) wrt composition
   * @note Assumes that pressure and temperature are strictly positive
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arrayView1d< real64 const > const composition,
                              arrayView1d< real64 const > const criticalPressure,
                              arrayView1d< real64 const > const criticalTemperature,
                              arrayView1d< real64 const > const acentricFactor,
                              real64 const & binaryInteractionCoefficients,
                              arraySlice1d< real64 const > const aPureCoefficient,
                              arraySlice1d< real64 const > const bPureCoefficient,
                              real64 const aMixtureCoefficient,
                              real64 const bMixtureCoefficient,
                              real64 & daMixtureCoefficient_dp,
                              real64 & dbMixtureCoefficient_dp,
                              real64 & daMixtureCoefficient_dt,
                              real64 & dbMixtureCoefficient_dt,
                              arraySlice1d< real64 > const daMixtureCoefficient_dz,
                              arraySlice1d< real64 > const dbMixtureCoefficient_dz );

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
                                real64 const & binaryInteractionCoefficients,
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
   * @param[in] daMixtureCoefficient_dp derivative of mixture coefficient (A) wrt pressure
   * @param[in] dbMixtureCoefficient_dp derivative of mixture coefficient (B) wrt pressure
   * @param[in] daMixtureCoefficient_dt derivative of mixture coefficient (A) wrt temperature
   * @param[in] dbMixtureCoefficient_dt derivative of mixture coefficient (B) wrt temperature
   * @param[in] daMixtureCoefficient_dz derivative of mixture coefficient (A) wrt composition
   * @param[in] dbMixtureCoefficient_dz derivative of mixture coefficient (B) wrt composition
   * @param[out] dCompressibilityFactor_dp derivative of the compressibility factor wrt pressure
   * @param[out] dCompressibilityFactor_dt derivative of the compressibility factor wrt temperature
   * @param[out] dCompressibilityFactor_dz derivative of the compressibility factor wrt composition
   * @note Assumes that pressure and temperature are strictly positive
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                real64 const & compressibilityFactor,
                                real64 const & daMixtureCoefficient_dp,
                                real64 const & dbMixtureCoefficient_dp,
                                real64 const & daMixtureCoefficient_dt,
                                real64 const & dbMixtureCoefficient_dt,
                                arraySlice1d< real64 const > const daMixtureCoefficient_dz,
                                arraySlice1d< real64 const > const dbMixtureCoefficient_dz,
                                real64 & dCompressibilityFactor_dp,
                                real64 & dCompressibilityFactor_dt,
                                arraySlice1d< real64 > const dCompressibilityFactor_dz );

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
                                  real64 const & binaryInteractionCoefficients,
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
compute( integer const numComps,
         real64 const & pressure,
         real64 const & temperature,
         arrayView1d< real64 const > const composition,
         arrayView1d< real64 const > const criticalPressure,
         arrayView1d< real64 const > const criticalTemperature,
         arrayView1d< real64 const > const acentricFactor,
         real64 const & binaryInteractionCoefficients,
         arraySlice1d< real64 > const logFugacityCoefficients )
{
  // step 0: allocate the stack memory needed for the update
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > aPureCoefficient( numComps );
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > bPureCoefficient( numComps );
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  real64 compressibilityFactor = 0.0;

  // step 1: compute the mixture coefficients aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              criticalPressure, // user input
                              criticalTemperature,
                              acentricFactor,
                              binaryInteractionCoefficients,
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
                            arrayView1d< real64 const > const criticalPressure,
                            arrayView1d< real64 const > const criticalTemperature,
                            arrayView1d< real64 const > const acentricFactor,
                            real64 const & binaryInteractionCoefficients,
                            arraySlice1d< real64 > const aPureCoefficient,
                            arraySlice1d< real64 > const bPureCoefficient,
                            real64 & aMixtureCoefficient,
                            real64 & bMixtureCoefficient )
{
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
      aMixtureCoefficient += ( composition[ic] * composition[jc] * ( 1.0 - binaryInteractionCoefficients ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] ) );
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
                            arrayView1d< real64 const > const criticalPressure,
                            arrayView1d< real64 const > const criticalTemperature,
                            arrayView1d< real64 const > const acentricFactor,
                            real64 const & binaryInteractionCoefficients,
                            arraySlice1d< real64 const > const aPureCoefficient,
                            arraySlice1d< real64 const > const bPureCoefficient,
                            real64 const aMixtureCoefficient,
                            real64 const bMixtureCoefficient,
                            real64 & daMixtureCoefficient_dp,
                            real64 & dbMixtureCoefficient_dp,
                            real64 & daMixtureCoefficient_dt,
                            real64 & dbMixtureCoefficient_dt,
                            arraySlice1d< real64 > const daMixtureCoefficient_dz,
                            arraySlice1d< real64 > const dbMixtureCoefficient_dz )
{
  GEOS_UNUSED_VAR( criticalPressure );
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > daPureCoefficient_dx( numComps );
  // Calculate pressure derivatives
  daMixtureCoefficient_dp = aMixtureCoefficient / pressure;
  dbMixtureCoefficient_dp = bMixtureCoefficient / pressure;
  // Calculate temperature derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const m = EOS_TYPE::evaluate( acentricFactor[ic] );
    real64 const sqrtTr = sqrt( temperature / criticalTemperature[ic] );
    real64 const mt = 1.0 + m * (1.0 - sqrtTr);
    daPureCoefficient_dx[ic] = -aPureCoefficient[ic] * (2.0/temperature + m/(mt*sqrtTr*criticalTemperature[ic]));
  }
  daMixtureCoefficient_dt = 0.0;
  dbMixtureCoefficient_dt = -bMixtureCoefficient / temperature;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const coeff = composition[ic] * composition[jc] * ( 1.0 - binaryInteractionCoefficients ) / sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      daMixtureCoefficient_dt += 0.5 * coeff * (daPureCoefficient_dx[ic]*aPureCoefficient[jc] + daPureCoefficient_dx[jc]*aPureCoefficient[ic]);
    }
  }
  // Calculate composition derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    daMixtureCoefficient_dz[ic] = 0.0;
    dbMixtureCoefficient_dz[ic] = 0.0;
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const coeff = ( 1.0 - binaryInteractionCoefficients ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      daMixtureCoefficient_dz[ic] += coeff * composition[jc];
      daMixtureCoefficient_dz[jc] += coeff * composition[ic];
    }
    dbMixtureCoefficient_dz[ic] = bPureCoefficient[ic];
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              arrayView1d< real64 const > const composition,
                              real64 const & binaryInteractionCoefficients,
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
                              real64 const & daMixtureCoefficient_dp,
                              real64 const & dbMixtureCoefficient_dp,
                              real64 const & daMixtureCoefficient_dt,
                              real64 const & dbMixtureCoefficient_dt,
                              arraySlice1d< real64 const > const daMixtureCoefficient_dz,
                              arraySlice1d< real64 const > const dbMixtureCoefficient_dz,
                              real64 & dcompressibilityFactor_dp,
                              real64 & dcompressibilityFactor_dt,
                              arraySlice1d< real64 > const dcompressibilityFactor_dz )
{
  // a Z3 + b Z2 + cZ + d = 0
  // Derivatives for a,b,c,d
  // dadx is zero;
  real64 dbdx = 0.0;
  real64 dcdx = 0.0;
  real64 dddx = 0.0;

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

  // Pressure derivatives
  dbdx = ( d1pd2 - 1.0 ) * dbMixtureCoefficient_dp;
  dcdx = daMixtureCoefficient_dp + (2.0*(d1xd2-d1pd2)*bMixtureCoefficient-d1pd2)*dbMixtureCoefficient_dp;
  dddx = -(aMixtureCoefficient*dbMixtureCoefficient_dp + daMixtureCoefficient_dp*bMixtureCoefficient
           + d1xd2*((3.0*bMixtureCoefficient+2.0)*bMixtureCoefficient*dbMixtureCoefficient_dp));
  dcompressibilityFactor_dp = (((dbdx*compressibilityFactor) + dcdx)*compressibilityFactor + dddx) * scalingFactor;

  // Temperature derivatives
  dbdx = ( d1pd2 - 1.0 ) * dbMixtureCoefficient_dt;
  dcdx = daMixtureCoefficient_dt + (2.0*(d1xd2-d1pd2)*bMixtureCoefficient-d1pd2)*dbMixtureCoefficient_dt;
  dddx = -(aMixtureCoefficient*dbMixtureCoefficient_dt + daMixtureCoefficient_dt*bMixtureCoefficient
           + d1xd2*((3.0*bMixtureCoefficient+2.0)*bMixtureCoefficient*dbMixtureCoefficient_dt));
  dcompressibilityFactor_dt = (((dbdx*compressibilityFactor) + dcdx)*compressibilityFactor + dddx) * scalingFactor;

  // Composition derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dbdx = ( d1pd2 - 1.0 ) * dbMixtureCoefficient_dz[ic];
    dcdx = daMixtureCoefficient_dz[ic] + (2.0*(d1xd2-d1pd2)*bMixtureCoefficient-d1pd2)*dbMixtureCoefficient_dz[ic];
    dddx = -(aMixtureCoefficient*dbMixtureCoefficient_dz[ic] + daMixtureCoefficient_dz[ic]*bMixtureCoefficient
             + d1xd2*((3.0*bMixtureCoefficient+2.0)*bMixtureCoefficient*dbMixtureCoefficient_dz[ic]));
    dcompressibilityFactor_dz[ic] = (((dbdx*compressibilityFactor) + dcdx)*compressibilityFactor + dddx) * scalingFactor;
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                arrayView1d< real64 const > const composition,
                                real64 const & binaryInteractionCoefficients,
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
      ki[ic] += composition[jc] * ( 1.0 - binaryInteractionCoefficients ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
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


} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
