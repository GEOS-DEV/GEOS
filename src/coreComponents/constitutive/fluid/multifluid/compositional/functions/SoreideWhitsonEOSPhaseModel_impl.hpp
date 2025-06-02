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
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 const & salinity,
                         real64 & aCoefficient,
                         real64 & bCoefficient )
{
  real64 daCoefficient_dp = 0.0;
  real64 dbCoefficient_dp = 0.0;
  real64 daCoefficient_dt = 0.0;
  real64 dbCoefficient_dt = 0.0;
  computePureCoefficientsAndDerivs( ic,
                                    pressure,
                                    temperature,
                                    componentProperties,
                                    salinity,
                                    aCoefficient,
                                    bCoefficient,
                                    daCoefficient_dp,
                                    dbCoefficient_dp,
                                    daCoefficient_dt,
                                    dbCoefficient_dt );
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
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
                                  real64 & dbCoefficient_dt )
{
  arraySlice1d< integer const > const & componentType = componentProperties.m_componentType;
  if( !isComponentType( componentType[ic], ComponentType::Water ))
  {
    CubicModel::computePureCoefficients( ic,
                                         pressure,
                                         temperature,
                                         componentProperties,
                                         aCoefficient,
                                         bCoefficient,
                                         daCoefficient_dp,
                                         dbCoefficient_dp,
                                         daCoefficient_dt,
                                         dbCoefficient_dt );
    return;
  }
  arraySlice1d< real64 const > const & criticalPressure = componentProperties.m_componentCriticalPressure;
  arraySlice1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;

  real64 const pr = pressure / criticalPressure[ic];
  real64 const tr = temperature / criticalTemperature[ic];

  real64 const trToMinus2 = 1.0/(tr*tr);
  real64 const trToMinus3 = trToMinus2/tr;

  real64 const csw = power( salinity, 1.1 );

  real64 const sqrtAlpha = 1.0 + 0.4530*(1.0 - tr*(1.0 - 0.0103*csw)) + 0.0034*(trToMinus3 - 1.0);
  real64 const dSqrtAlpha_dtr = -0.4530*(1.0 - 0.0103*csw) - 0.0102*trToMinus3/tr;
  real64 const alpha = sqrtAlpha * sqrtAlpha;
  real64 const dAlpha_dt = 2.0 * sqrtAlpha * dSqrtAlpha_dtr / criticalTemperature[ic];

  aCoefficient = EOS_TYPE::omegaA * pr * trToMinus2 * alpha;
  bCoefficient = EOS_TYPE::omegaB * pr / tr;

  daCoefficient_dp = aCoefficient / pressure;
  dbCoefficient_dp = bCoefficient / pressure;

  daCoefficient_dt = EOS_TYPE::omegaA * pr * (-2.0*trToMinus3 * alpha / criticalTemperature[ic] + trToMinus2 * dAlpha_dt);
  dbCoefficient_dt = -bCoefficient / temperature;
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

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d< real64 const, USD > const & composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            real64 const & salinity,
                            arraySlice1d< real64 > const & aPureCoefficient,
                            arraySlice1d< real64 > const & bPureCoefficient,
                            real64 & aMixtureCoefficient,
                            real64 & bMixtureCoefficient )
{
  // Pure component coefficients
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficients( ic, pressure, temperature, componentProperties, salinity, aPureCoefficient[ic], bPureCoefficient[ic] );
  }

  // Mixture coefficients
  aMixtureCoefficient = 0.0;
  for( integer jc = 0; jc < numComps; ++jc )
  {
    real64 const aj = sqrt( aPureCoefficient[jc] );
    aMixtureCoefficient += composition[jc] * composition[jc] * aj * aj;
    for( integer ic = jc+1; ic < numComps; ++ic )
    {
      real64 const ai = sqrt( aPureCoefficient[ic] );
      real64 kij = 0.0;
      real64 dkij_dT = 0.0;
      getBinaryInteractionCoefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij, dkij_dT );
      aMixtureCoefficient += 2.0 * composition[ic] * composition[jc] * ( 1.0 - kij ) * ai * aj;
    }
  }

  bMixtureCoefficient = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    bMixtureCoefficient += composition[ic] * bPureCoefficient[ic];
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
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
                                 arraySlice1d< real64 > const & bMixtureCoefficientDerivs )
{
  // a parameter derivatives
  aMixtureCoefficientDerivs[Deriv::dP] = aMixtureCoefficient / pressure;

  real64 aCoefficient = 0.0;
  real64 bCoefficient = 0.0;
  real64 dummy = 0.0;
  stackArray1d< real64, maxNumComps > daPureCoefficient_dt( numComps );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficientsAndDerivs( ic, pressure, temperature, componentProperties, salinity,
                                      aCoefficient, bCoefficient, dummy, dummy, daPureCoefficient_dt[ic], dummy );
  }

  aMixtureCoefficientDerivs[Deriv::dT] = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    aMixtureCoefficientDerivs[Deriv::dC+ic] = 0.0;
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 kij = 0.0;
      real64 dkij_dT = 0.0;
      getBinaryInteractionCoefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij, dkij_dT );

      real64 const aij = sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      real64 const coeff = composition[ic] * composition[jc] * ( 1.0 - kij ) / aij;
      aMixtureCoefficientDerivs[Deriv::dT] += 0.5 * coeff * (daPureCoefficient_dt[ic]*aPureCoefficient[jc] + daPureCoefficient_dt[jc]*aPureCoefficient[ic]);
      aMixtureCoefficientDerivs[Deriv::dT] -= composition[ic] * composition[jc] * dkij_dT * aij;

      aMixtureCoefficientDerivs[Deriv::dC+ic] += composition[jc] * ( 1.0 - kij ) * aij;
      aMixtureCoefficientDerivs[Deriv::dC+jc] += composition[ic] * ( 1.0 - kij ) * aij;
    }
  }

  // b parameter derivatives
  bMixtureCoefficientDerivs[Deriv::dP] = bMixtureCoefficient / pressure;
  bMixtureCoefficientDerivs[Deriv::dT] = -bMixtureCoefficient / temperature;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    bMixtureCoefficientDerivs[Deriv::dC+ic] = bPureCoefficient[ic];
  }
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
  stackArray2d< real64, 2*maxNumComps > tempData( 2, numComps );
  arraySlice1d< real64 > aPureCoefficient = tempData[0];
  arraySlice1d< real64 > bPureCoefficient = tempData[1];
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;

  // Step 1: compute the mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              salinity,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // Step 2: Extract binary interaction coefficients
  stackArray2d< real64, maxNumComps *maxNumComps > kij( numComps, numComps );
  real64 dkij_dT = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    kij( ic, ic ) = 0.0;
    for( integer jc = ic+1; jc < numComps; ++jc )
    {
      getBinaryInteractionCoefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij( ic, jc ), dkij_dT );
      kij( jc, ic ) = kij( ic, jc );
    }
  }

  // Step 3: Calculate the compressibility factor from the underlying cubic EOS
  CubicModel::computeCompressibilityFactor( numComps,
                                            composition,
                                            kij.toSliceConst(),
                                            aPureCoefficient,
                                            bPureCoefficient,
                                            aMixtureCoefficient,
                                            bMixtureCoefficient,
                                            compressibilityFactor );
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
  integer constexpr numMaxDofs = maxNumComps + 2;

  stackArray2d< real64, 4*numMaxDofs > tempData( 4, numMaxDofs );
  arraySlice1d< real64 > aPureCoefficient = tempData[0];
  arraySlice1d< real64 > bPureCoefficient = tempData[1];
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  arraySlice1d< real64 > aMixtureCoefficientDerivs = tempData[2];
  arraySlice1d< real64 > bMixtureCoefficientDerivs = tempData[3];

  // step 1: compute the mixture coefficients
  // 1.1: Compute the pure and mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              salinity,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // 1.2: Compute the pure and mixture coefficient derivatives
  computeMixtureCoefficientDerivs( numComps,
                                   pressure,
                                   temperature,
                                   composition,
                                   componentProperties,
                                   salinity,
                                   aPureCoefficient.toSliceConst(),
                                   bPureCoefficient.toSliceConst(),
                                   aMixtureCoefficient,
                                   bMixtureCoefficient,
                                   aMixtureCoefficientDerivs, // output
                                   bMixtureCoefficientDerivs );

  // 2.1: Update the compressibility factor
  computeCompressibilityFactor( numComps,
                                pressure,
                                temperature,
                                composition,
                                componentProperties,
                                salinity,
                                compressibilityFactor );

  // 2.2: Calculate derivatives
  CubicModel::computeCompressibilityFactor( numComps,
                                            aMixtureCoefficient,
                                            bMixtureCoefficient,
                                            compressibilityFactor,
                                            aMixtureCoefficientDerivs.toSliceConst(),
                                            bMixtureCoefficientDerivs.toSliceConst(),
                                            compressibilityFactorDerivs );
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
  // step 0: allocate the stack memory needed for the update
  stackArray2d< real64, 2*maxNumComps > pureCoefficients( 2, numComps );
  arraySlice1d< real64 > aPureCoefficient = pureCoefficients[0];
  arraySlice1d< real64 > bPureCoefficient = pureCoefficients[1];
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  real64 compressibilityFactor = 0.0;

  // step 1: calculate the dynamic binary interaction coefficients
  stackArray2d< real64, maxNumComps *maxNumComps > kij( numComps, numComps );
  real64 dkij_dT = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    kij( ic, ic ) = 0.0;
    for( integer jc = ic+1; jc < numComps; ++jc )
    {
      getBinaryInteractionCoefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij( ic, jc ), dkij_dT );
      kij( jc, ic ) = kij( ic, jc );
    }
  }

  // step 2: compute the mixture coefficients
  computeMixtureCoefficients( numComps,
                              pressure,
                              temperature,
                              composition,
                              componentProperties,
                              salinity,
                              aPureCoefficient,
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // step 3: compute the compressibility factor
  CubicModel::computeCompressibilityFactor( numComps,
                                            composition,
                                            kij.toSliceConst(),
                                            aPureCoefficient,
                                            bPureCoefficient,
                                            aMixtureCoefficient,
                                            bMixtureCoefficient,
                                            compressibilityFactor );

  // step 4: use mixture coefficients and compressibility factor to update fugacity coefficients
  CubicModel::computeLogFugacityCoefficients( numComps,
                                              composition,
                                              kij.toSliceConst(),
                                              compressibilityFactor,
                                              aPureCoefficient,
                                              bPureCoefficient,
                                              aMixtureCoefficient,
                                              bMixtureCoefficient,
                                              logFugacityCoefficients );
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficientDerivs( integer const numComps,
                                     real64 const & pressure,
                                     real64 const & temperature,
                                     arraySlice1d< real64 const, USD > const & composition,
                                     ComponentProperties::KernelWrapper const & componentProperties,
                                     real64 const & salinity,
                                     arraySlice1d< real64 const > const & logFugacityCoefficients,
                                     arraySlice2d< real64 > const & logFugacityCoefficientDerivs )
{
  integer constexpr numMaxComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  integer constexpr numMaxDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  integer const numDofs = 2 + numComps;

  GEOS_UNUSED_VAR( logFugacityCoefficients );

  stackArray3d< real64, 3*2*numMaxComps > pureCoefficients( 3, 2, numComps );
  arraySlice1d< real64 > aPureCoefficient = pureCoefficients[0][0];
  arraySlice1d< real64 > bPureCoefficient = pureCoefficients[0][1];
  arraySlice2d< real64 > aPureCoefficientDerivs = pureCoefficients[1];
  arraySlice2d< real64 > bPureCoefficientDerivs = pureCoefficients[2];
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  real64 compressibilityFactor = 0.0;
  stackArray2d< real64, 3*numMaxDofs > scalarDerivs( 3, numDofs );
  arraySlice1d< real64 > aMixtureCoefficientDerivs = scalarDerivs[0];
  arraySlice1d< real64 > bMixtureCoefficientDerivs = scalarDerivs[1];
  arraySlice1d< real64 > compressibilityFactorDerivs = scalarDerivs[2];

  // 1: Compute pure coefficients and derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficientsAndDerivs( ic,
                                      pressure,
                                      temperature,
                                      componentProperties,
                                      salinity,
                                      aPureCoefficient[ic],
                                      bPureCoefficient[ic],
                                      aPureCoefficientDerivs( Deriv::dP, ic ),
                                      bPureCoefficientDerivs( Deriv::dP, ic ),
                                      aPureCoefficientDerivs( Deriv::dT, ic ),
                                      bPureCoefficientDerivs( Deriv::dT, ic ));
  }

  // 2.1: Compute mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              salinity,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // 2.2: Compute mixture coefficient derivatives
  computeMixtureCoefficientDerivs( numComps,
                                   pressure,
                                   temperature,
                                   composition,
                                   componentProperties,
                                   salinity,
                                   aPureCoefficient,
                                   bPureCoefficient,
                                   aMixtureCoefficient,
                                   bMixtureCoefficient,
                                   aMixtureCoefficientDerivs,
                                   bMixtureCoefficientDerivs );

  // 3: Extract binary interaction coefficients
  stackArray1d< real64, maxNumComps > ki( numComps );
  stackArray2d< real64, numMaxComps * numMaxDofs > dki( numComps, numDofs );
  real64 dkij_dT = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dki( ic, ic ) = 0.0;
    for( integer jc = ic+1; jc < numComps; ++jc )
    {
      getBinaryInteractionCoefficient( pressure, temperature, componentProperties, salinity, ic, jc, dki( ic, jc ), dkij_dT );
      dki( jc, ic ) = dki( ic, jc );
    }
  }

  // 4.1: Calculate the compressibility factor from the underlying cubic EOS
  CubicModel::computeCompressibilityFactor( numComps,
                                            composition,
                                            dki.toSliceConst(),
                                            aPureCoefficient,
                                            bPureCoefficient,
                                            aMixtureCoefficient,
                                            bMixtureCoefficient,
                                            compressibilityFactor );

  // 4.2: Compute the compressibility factor derivatives
  CubicModel::computeCompressibilityFactor( numComps,
                                            aMixtureCoefficient,
                                            bMixtureCoefficient,
                                            compressibilityFactor,
                                            aMixtureCoefficientDerivs.toSliceConst(),
                                            bMixtureCoefficientDerivs.toSliceConst(),
                                            compressibilityFactorDerivs );


  // 5. Calculate interaction parameters
  // ki
  for( integer ic = 0; ic < numComps; ++ic )
  {
    ki[ic] = 0.0;
    dki( ic, Deriv::dP ) = 0.0;
    dki( ic, Deriv::dT ) = 0.0;
  }

  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const Ai = aPureCoefficient[ic];
    real64 const ai = sqrt( Ai );
    real64 const dai_dP = aPureCoefficientDerivs( Deriv::dP, ic );
    real64 const dai_dT = aPureCoefficientDerivs( Deriv::dT, ic );

    ki[ic] += composition[ic] * Ai;
    dki( ic, Deriv::dC + ic ) = Ai;

    dki( ic, Deriv::dP ) += composition[ic] * dai_dP;
    dki( ic, Deriv::dT ) += composition[ic] * dai_dT;

    for( integer jc = ic+1; jc < numComps; ++jc )
    {
      real64 const aj = sqrt( aPureCoefficient[jc] );
      real64 const daj_dP = aPureCoefficientDerivs( Deriv::dP, jc );
      real64 const daj_dT = aPureCoefficientDerivs( Deriv::dT, jc );

      real64 kij = 0.0;
      getBinaryInteractionCoefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij, dkij_dT );

      real64 const bicValue = ( 1.0 - kij ) * ai * aj;

      ki[ic] += composition[jc] * bicValue;
      ki[jc] += composition[ic] * bicValue;

      dki( ic, Deriv::dC + jc ) = bicValue;
      dki( jc, Deriv::dC + ic ) = bicValue;

      dki( ic, Deriv::dP ) += 0.5 * composition[jc] * bicValue * (dai_dP/aPureCoefficient[ic] + daj_dP/aPureCoefficient[jc]);
      dki( jc, Deriv::dP ) += 0.5 * composition[ic] * bicValue * (dai_dP/aPureCoefficient[ic] + daj_dP/aPureCoefficient[jc]);

      dki( ic, Deriv::dT ) += 0.5 * composition[jc] * bicValue * (dai_dT/aPureCoefficient[ic] + daj_dT/aPureCoefficient[jc])
                              - composition[jc] * dkij_dT * ai * aj;
      dki( jc, Deriv::dT ) += 0.5 * composition[ic] * bicValue * (dai_dT/aPureCoefficient[ic] + daj_dT/aPureCoefficient[jc])
                              - composition[ic] * dkij_dT * ai * aj;
    }
  }

  auto const calculateDerivatives = [&]( integer const idof ){
    real64 constexpr d1 = EOS_TYPE::delta1;
    real64 constexpr d2 = EOS_TYPE::delta2;
    real64 const E = log( compressibilityFactor + d1 * bMixtureCoefficient )
                     - log( compressibilityFactor + d2 * bMixtureCoefficient );

    real64 const dE_dX = (compressibilityFactorDerivs[idof] + d1*bMixtureCoefficientDerivs[idof]) /
                         ( compressibilityFactor + d1 * bMixtureCoefficient )
                         -(compressibilityFactorDerivs[idof] + d2*bMixtureCoefficientDerivs[idof]) /
                         ( compressibilityFactor + d2 * bMixtureCoefficient );

    //real64 const F = log( compressibilityFactor - bMixtureCoefficient );
    real64 const dF_dX = (compressibilityFactorDerivs[idof] - bMixtureCoefficientDerivs[idof]) / (compressibilityFactor - bMixtureCoefficient);

    real64 const G = 1.0 / ( ( d1 - d2 ) * bMixtureCoefficient );
    real64 const dG_dX = -G * bMixtureCoefficientDerivs[idof] / bMixtureCoefficient;

    real64 const A = aMixtureCoefficient;
    real64 const dA_dX = aMixtureCoefficientDerivs[idof];

    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const B = bPureCoefficient[ic] / bMixtureCoefficient;
      real64 dB_dX = -B*bMixtureCoefficientDerivs[idof] / bMixtureCoefficient;
      if( idof < Deriv::dC )
      {
        dB_dX += bPureCoefficientDerivs( idof, ic ) / bMixtureCoefficient;
      }

      // lnPhi = ( compressibilityFactor - 1 ) * B - F - G * ( 2 * ki[ic] - A * B ) * E;
      logFugacityCoefficientDerivs( ic, idof ) =
        compressibilityFactorDerivs[idof]*B + ( compressibilityFactor - 1 ) * dB_dX
        - dF_dX
        - dG_dX * ( 2 * ki[ic] - A * B ) * E
        - G * ( 2 * dki( ic, idof ) - dA_dX * B - A * dB_dX ) * E
        - G * ( 2 * ki[ic] - A * B ) * dE_dX;
    }
  };

  calculateDerivatives( Deriv::dP );
  calculateDerivatives( Deriv::dT );
  for( integer jc = 0; jc < numComps; ++jc )
  {
    calculateDerivatives( Deriv::dC+jc );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSPHASEMODEL_IMPL_HPP_
