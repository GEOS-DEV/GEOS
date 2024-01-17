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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "LohrenzBrayClarkViscosity.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

namespace geos
{

namespace constitutive
{

using Deriv = multifluid::DerivativeOffset;

namespace compositional
{

LohrenzBrayClarkViscosityUpdate::LohrenzBrayClarkViscosityUpdate( MixingType const mixing_type )
  : m_mixing_type( mixing_type )
{}

LohrenzBrayClarkViscosity::LohrenzBrayClarkViscosity( string const & name,
                                                      ComponentProperties const & componentProperties ):
  FunctionBase( name, componentProperties )
{}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_mixing_type );
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                               real64 const & pressure,
                                               real64 const & temperature,
                                               arraySlice1d< real64 const > const & phaseComposition,
                                               arraySlice2d< real64 const > const & dPhaseComposition,
                                               real64 const & density,
                                               arraySlice1d< real64 const > const & dDensity,
                                               real64 & viscosity,
                                               arraySlice1d< real64 > const & dViscosity,
                                               bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure );   // No-direct pressure dependence (instead through density)
  GEOS_UNUSED_VAR( useMass );

  integer constexpr maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  integer constexpr maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  integer const numComponents = componentProperties.m_componentMolarWeight.size();
  integer const numDofs = numComponents + 2;

  // Space for temporary derivatives
  stackArray1d< real64, maxNumDofs > tempDerivs( numDofs );

  // Estimate pure component properties at dilute-gas conditions (pressure near atmospheric) using
  // Stiel and Thodos [1961] correlation: https://doi.org/10.1002/aic.690070416
  // Dilute viscosity is solely temperature dependent
  // Units are converted so componentViscosity is in centipoise to match original reference

  stackArray1d< real64, maxNumComps > componentDiluteViscosity( numComponents );

  computeComponentDiluteViscosity_StielThodos( numComponents,
                                               componentProperties,
                                               temperature,
                                               componentDiluteViscosity,
                                               tempDerivs );

  // Estimate phase viscosity (in cp) at dilute gas conditions using either
  // the Herning and Zipperer [1936], Wilke [1950], or Brokaw[1968] mixture rule.
  // The classic LBC model uses Herning-Zipperer, but the other two may be more accurate.

  if( m_mixing_type == MixingType::HERNING_ZIPPERER )
  {
    computePhaseDiluteViscosity_HerningZipperer( numComponents,
                                                 componentProperties,
                                                 temperature,
                                                 phaseComposition,
                                                 componentDiluteViscosity,
                                                 tempDerivs,
                                                 viscosity,
                                                 dViscosity );
  }
  else if( m_mixing_type == MixingType::WILKE )
  {
    computePhaseDiluteViscosity_Wilke( numComponents,
                                       componentProperties,
                                       temperature,
                                       phaseComposition,
                                       componentDiluteViscosity,
                                       tempDerivs,
                                       viscosity,
                                       dViscosity );
  }
  else if( m_mixing_type == MixingType::BROKAW )
  {
    computePhaseDiluteViscosity_Brokaw( numComponents,
                                        componentProperties,
                                        temperature,
                                        phaseComposition,
                                        componentDiluteViscosity,
                                        tempDerivs,
                                        viscosity,
                                        dViscosity );
  }

  // Estimate phase viscosity at given (P,T) conditions using LBC [1964] correlation.
  // This is an additional term added to the dilute gas estimate above.

  computePhaseViscosity_LohrenzBrayClark( numComponents,
                                          componentProperties,
                                          phaseComposition,
                                          density,
                                          dDensity,
                                          viscosity,
                                          dViscosity );

  // Convert derivatives from phase to total composition
  convertDerivativesToTotalMoleFraction( numComponents, dPhaseComposition, dViscosity, tempDerivs );

  // Scale centipoise to pascal.seconds
  viscosity *= CP_TO_PAS;
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    dViscosity[kc] *= CP_TO_PAS;
  }
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::computeComponentDiluteViscosity_StielThodos( integer const numComponents,
                                                                                   ComponentProperties::KernelWrapper const & componentProperties,
                                                                                   real64 const temperature,
                                                                                   arraySlice1d< real64 > const componentDiluteViscosity,
                                                                                   arraySlice1d< real64 > const dComponentDiluteViscosity_dTemperature ) const
{
  for( integer ic = 0; ic < numComponents; ++ic )
  {
    real64 const criticalPressure = componentProperties.m_componentCriticalPressure[ic];
    real64 const criticalTemperature = componentProperties.m_componentCriticalTemperature[ic];
    real64 const molarWeight = componentProperties.m_componentMolarWeight[ic];

    real64 reducedTemperature = temperature / criticalTemperature;
    real64 inverseComponentChi = 0.0;
    real64 discardDerivative = 0.0;
    inverseChiParameter( criticalPressure, criticalTemperature, molarWeight, inverseComponentChi,
                         discardDerivative, discardDerivative, discardDerivative );

    if( molarWeight < 2.1e-3 )  // hydrogen correlation, Stiel & Thodos, 1961, Eq. 12
    {
      componentDiluteViscosity[ic] = 90.71e-5 * pow( 0.1375*temperature - 1.67, 0.625 );
      dComponentDiluteViscosity_dTemperature[ic] = 90.71e-5 * 0.625 * 0.1375 * pow( 0.1375*temperature - 1.67, -0.375 );
    }
    else if( reducedTemperature <= 1.5 )   // nonpolar gas correlation at low temp, Eq. 9
    {
      componentDiluteViscosity[ic] = 34e-5 * pow( reducedTemperature, 0.94 ) * inverseComponentChi;
      dComponentDiluteViscosity_dTemperature[ic] = 34e-5 * 0.94 * pow( reducedTemperature, -0.06 ) * inverseComponentChi / criticalTemperature;
    }
    else   // nonpolar gas correlation at high temp, Eq. 10
    {
      componentDiluteViscosity[ic] = 17.78e-5 * pow( 4.58*reducedTemperature-1.67, 0.625 ) * inverseComponentChi;
      dComponentDiluteViscosity_dTemperature[ic] = 17.78e-5 * 4.58 * 0.625 * pow( 4.58*reducedTemperature-1.67, -0.375 ) * inverseComponentChi / criticalTemperature;
    }
  }
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::computePhaseDiluteViscosity_HerningZipperer( integer const numComponents,
                                                                                   ComponentProperties::KernelWrapper const & componentProperties,
                                                                                   real64 const temperature,
                                                                                   arraySlice1d< real64 const > const & phaseComposition,
                                                                                   arraySlice1d< real64 const > const & componentDiluteViscosity,
                                                                                   arraySlice1d< real64 const > const & dComponentDiluteViscosity_dTemperature,
                                                                                   real64 & phaseViscosity,
                                                                                   arraySlice1d< real64 > const & dPhaseViscosity ) const
{
  GEOS_UNUSED_VAR( temperature );

  real64 A = 0.0;
  real64 dA_dT = 0.0;
  real64 B = 0.0;

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    real64 const sqrtMolarWeight = sqrt( componentProperties.m_componentMolarWeight[ic] );
    A += phaseComposition[ic] * sqrtMolarWeight * componentDiluteViscosity[ic];
    B += phaseComposition[ic] * sqrtMolarWeight;
    dA_dT += phaseComposition[ic] * sqrtMolarWeight * dComponentDiluteViscosity_dTemperature[ic];
  }

  phaseViscosity = A/B;
  dPhaseViscosity[Deriv::dP] = 0.0;
  dPhaseViscosity[Deriv::dT] = dA_dT/B;

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    real64 const sqrtMolarWeight = sqrt( componentProperties.m_componentMolarWeight[ic] );
    dPhaseViscosity[Deriv::dC+ic] = sqrtMolarWeight * componentDiluteViscosity[ic] / B - A * sqrtMolarWeight / (B * B);
  }
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::computePhaseDiluteViscosity_Wilke( integer const numComponents,
                                                                         ComponentProperties::KernelWrapper const & componentProperties,
                                                                         real64 const temperature,
                                                                         arraySlice1d< real64 const > const & phaseComposition,
                                                                         arraySlice1d< real64 const > const & componentDiluteViscosity,
                                                                         arraySlice1d< real64 const > const & dComponentDiluteViscosity_dTemperature,
                                                                         real64 & phaseViscosity,
                                                                         arraySlice1d< real64 > const & dPhaseViscosity ) const
{
  GEOS_UNUSED_VAR( temperature );

  // compute the "phi" interaction matrix (and its temperature derivatives)
  integer constexpr maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  stackArray2d< real64, maxNumComps *maxNumComps > phi( numComponents, numComponents );
  stackArray2d< real64, maxNumComps *maxNumComps > dPhi_dT( numComponents, numComponents );

  phi.zero();
  dPhi_dT.zero();

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    for( integer jc = 0; jc < numComponents; ++jc )
    {
      real64 const mw_i = componentProperties.m_componentMolarWeight[ic];
      real64 const mw_j = componentProperties.m_componentMolarWeight[jc];

      real64 const weightRatio = mw_i / mw_j;

      real64 const invVisc_j = 1.0 / componentDiluteViscosity[jc];
      real64 const viscosityRatio = componentDiluteViscosity[ic] * invVisc_j;

      real64 const dViscosityRatio_dT = dComponentDiluteViscosity_dTemperature[ic] * invVisc_j
                                        - componentDiluteViscosity[ic] * dComponentDiluteViscosity_dTemperature[jc] * invVisc_j * invVisc_j;

      real64 const A = 1.0 + sqrt( viscosityRatio )*pow( weightRatio, -0.25 );
      real64 const dA_dT = 0.5*pow( weightRatio, -0.25 )/sqrt( viscosityRatio )*dViscosityRatio_dT;

      real64 const B = A*A;
      real64 const dB_dT = 2.0*A*dA_dT;
      real64 const C = sqrt( 8.0 )*sqrt( 1.0 + weightRatio );

      phi( ic, jc ) = B/C;
      dPhi_dT( ic, jc ) = dB_dT/C;
    }
  }

  // compute phase viscosity via Wilke mixing rule

  phaseViscosity = 0;
  LvArray::forValuesInSlice( dPhaseViscosity, setZero );

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    real64 A = 0;
    real64 dA_dT = 0;

    for( integer jc = 0; jc < numComponents; ++jc )
    {
      A += phi( ic, jc )*phaseComposition[jc];
      dA_dT += dPhi_dT( ic, jc )*phaseComposition[jc];
    }

    phaseViscosity += phaseComposition[ic] * componentDiluteViscosity[ic] / A;
    dPhaseViscosity[Deriv::dT] += phaseComposition[ic] * dComponentDiluteViscosity_dTemperature[ic] / A
                                  - phaseComposition[ic] * componentDiluteViscosity[ic] / (A*A) * dA_dT;

    dPhaseViscosity[Deriv::dC+ic] += componentDiluteViscosity[ic] / A;

    // the following is some tricky loop merging.  the derivatives for other components "jc" will depend on the value of A
    // computed above for component "ic", so we add these entries for other derivatives immediately.

    for( integer jc = 0; jc < numComponents; ++jc )
    {
      dPhaseViscosity[Deriv::dC+jc] -= phaseComposition[ic] * componentDiluteViscosity[ic] / (A*A) * phi( ic, jc );
    }
  }
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::computePhaseDiluteViscosity_Brokaw( integer const numComponents,
                                                                          ComponentProperties::KernelWrapper const & componentProperties,
                                                                          real64 const temperature,
                                                                          arraySlice1d< real64 const > const & phaseComposition,
                                                                          arraySlice1d< real64 const > const & componentDiluteViscosity,
                                                                          arraySlice1d< real64 const > const & dComponentDiluteViscosity_dTemperature,
                                                                          real64 & phaseViscosity,
                                                                          arraySlice1d< real64 > const & dPhaseViscosity ) const
{
  GEOS_UNUSED_VAR( temperature );

  // Compute the "phi" interaction matrix (constant, as only function of molecular weights)
  integer constexpr maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  stackArray2d< real64, maxNumComps *maxNumComps > phi( numComponents, numComponents );
  phi.zero();

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    for( integer jc = 0; jc < ic; ++jc )
    {
      real64 const mw_i = componentProperties.m_componentMolarWeight[ic];
      real64 const mw_j = componentProperties.m_componentMolarWeight[jc];
      real64 const A = mw_i / mw_j;
      real64 const B = pow( 4.0*mw_i*mw_j/((mw_i+mw_j)*(mw_i+mw_j)), 0.25 );
      real64 const A45 = pow( A, 0.45 );

      phi( ic, jc ) = B/sqrt( A ) * (1.0 + (A-A45)/(2.0 + 2.0*A + B*(1+A45) / (1+B)));
      phi( jc, ic ) = phi( ic, jc );
    }
    phi( ic, ic ) = 1.0;
  }

  phaseViscosity = 0.0;
  LvArray::forValuesInSlice( dPhaseViscosity, setZero );

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    real64 A = 0;
    real64 dA_dT = 0;

    for( integer jc = 0; jc < numComponents; ++jc )
    {
      real64 const aij = phi( ic, jc )*phaseComposition[jc] / sqrt( componentDiluteViscosity[jc] );
      A += aij;
      dA_dT -= 0.5*aij*dComponentDiluteViscosity_dTemperature[jc] / componentDiluteViscosity[jc];
    }

    real64 const sqrtComponentDiluteViscosity = sqrt( componentDiluteViscosity[ic] );
    phaseViscosity += phaseComposition[ic] * sqrtComponentDiluteViscosity / A;

    dPhaseViscosity[Deriv::dT] += 0.5*phaseComposition[ic] / sqrtComponentDiluteViscosity * dComponentDiluteViscosity_dTemperature[ic] / A
                                  - phaseComposition[ic] * sqrtComponentDiluteViscosity / (A*A) * dA_dT;

    dPhaseViscosity[Deriv::dC+ic] += sqrtComponentDiluteViscosity / A;

    // the following is some tricky loop merging.  the derivatives for other components "jc" will depend on the value of A
    // computed above for component "ic", so we add these entries for other derivatives immediately.

    for( integer jc = 0; jc < numComponents; ++jc )
    {
      dPhaseViscosity[Deriv::dC+jc] -= phaseComposition[ic] * phi( ic, jc ) * sqrt( componentDiluteViscosity[ic]/componentDiluteViscosity[jc] ) / (A*A);
    }
  }
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::computePhaseViscosity_LohrenzBrayClark( integer const numComponents,
                                                                              ComponentProperties::KernelWrapper const & componentProperties,
                                                                              arraySlice1d< real64 const > const & phaseComposition,
                                                                              real64 const phaseDensity,
                                                                              arraySlice1d< real64 const > const & dPhaseDensity,
                                                                              real64 & phaseViscosity,
                                                                              arraySlice1d< real64 > const & dPhaseViscosity ) const
{
  // Compute phase pseudo properties via Kay's mixing rule
  real64 phaseCriticalPressure = 0.0;
  real64 phaseCriticalTemperature = 0.0;
  real64 phaseCriticalVolume = 0.0;
  real64 phaseMolarWeight = 0.0;

  auto const & criticalPressure    = componentProperties.m_componentCriticalPressure;
  auto const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
  auto const & criticalVolume      = componentProperties.m_componentCriticalVolume;
  auto const & molarWeight         = componentProperties.m_componentMolarWeight;

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    phaseCriticalPressure    += phaseComposition[ic] * criticalPressure[ic];
    phaseCriticalTemperature += phaseComposition[ic] * criticalTemperature[ic];
    phaseCriticalVolume      += phaseComposition[ic] * criticalVolume[ic];
    phaseMolarWeight         += phaseComposition[ic] * molarWeight[ic];
  }

  // Compute LBC polynomial

  real64 reducedDensity = phaseDensity * phaseCriticalVolume / phaseMolarWeight;
  real64 inversePhaseChi, dInversePhaseChi_dPc, dInversePhaseChi_dTc, dInversePhaseChi_dMw;

  inverseChiParameter( phaseCriticalPressure,
                       phaseCriticalTemperature,
                       phaseMolarWeight,
                       inversePhaseChi,
                       dInversePhaseChi_dPc,
                       dInversePhaseChi_dTc,
                       dInversePhaseChi_dMw );

  real64 polynomialOne =   0.1023000
                         + 0.0233640*reducedDensity
                         + 0.0585330*pow( reducedDensity, 2 )
                         - 0.0407580*pow( reducedDensity, 3 )
                         + 0.0093324*pow( reducedDensity, 4 );

  real64 polynomialTwo = pow( polynomialOne, 4.0 ) - 1.0e-4;

  // add polynomial contribution to dilute term to get final phase viscosity

  phaseViscosity += polynomialTwo * inversePhaseChi;

  // get derivatives, noting phaseDensity is a function of (pressure, temperature, composition)
  // and inversePhaseChi is a function of composition.
  // these derivatives are *added* to the ones already present from the dilute terms above.

  real64 dPolynomialOne_dReducedDensity = 0.023364 + 0.117066*reducedDensity - 0.122274*pow( reducedDensity, 2 ) + 0.0373296*pow( reducedDensity, 3 );
  real64 dPolynomialTwo_dReducedDensity = 4.0 * pow( polynomialOne, 3.0 ) * dPolynomialOne_dReducedDensity;

  real64 dViscosity_dDensity       = dPolynomialTwo_dReducedDensity * inversePhaseChi * phaseCriticalVolume / phaseMolarWeight;
  real64 dViscosity_dCriticalRatio = dPolynomialTwo_dReducedDensity * inversePhaseChi * phaseDensity;

  real64 dViscosity_dPc = polynomialTwo * dInversePhaseChi_dPc;
  real64 dViscosity_dTc = polynomialTwo * dInversePhaseChi_dTc;
  real64 dViscosity_dMw = polynomialTwo * dInversePhaseChi_dMw - dViscosity_dCriticalRatio * phaseCriticalVolume / pow( phaseMolarWeight, 2 );
  real64 dViscosity_dVc = dViscosity_dCriticalRatio / phaseMolarWeight;

  dPhaseViscosity[Deriv::dP] += dViscosity_dDensity * dPhaseDensity[Deriv::dP];
  dPhaseViscosity[Deriv::dT] += dViscosity_dDensity * dPhaseDensity[Deriv::dT];

  for( integer ic = 0; ic < numComponents; ++ic )
  {
    dPhaseViscosity[Deriv::dC+ic] += dViscosity_dDensity * dPhaseDensity[Deriv::dC+ic]
                                     + dViscosity_dPc * criticalPressure[ic]
                                     + dViscosity_dTc * criticalTemperature[ic]
                                     + dViscosity_dMw * molarWeight[ic]
                                     + dViscosity_dVc * criticalVolume[ic];
  }
}

GEOS_HOST_DEVICE
void LohrenzBrayClarkViscosityUpdate::inverseChiParameter( real64 const criticalPressure,
                                                           real64 const criticalTemperature,
                                                           real64 const molarWeight,
                                                           real64 & value,
                                                           real64 & derivP,
                                                           real64 & derivT,
                                                           real64 & derivM ) const
{
  real64 T  = pow( criticalTemperature, 1.0/6.0 );
  real64 dT = (1.0/6.0) * pow( criticalTemperature, -5.0/6.0 );

  real64 constexpr sqrt1000 = sqrt( 1000.0 ); // note: kg/mol to atomic mass units
  real64 M  = sqrt1000 * sqrt( molarWeight );
  real64 dM = 0.5 * sqrt1000 / sqrt( molarWeight );

  real64 P  = pow( criticalPressure / PA_TO_ATM, 2.0/3.0 );   // note: pascal to atm conversion
  real64 dP = pow( PA_TO_ATM, -2.0/3.0 ) * pow( criticalPressure, -1.0/3.0 ) * 2.0/3.0;

  value  = M*P/T;
  derivP = M*dP/T;
  derivT = -M*P*dT/(T*T);
  derivM = dM*P/T;
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
