/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SurfaceInformedPolymerHelpers.hpp
 * @brief Shared scalar functions for the surface-informed polymer continuum and cohesive-zone models.
 *
 * The functions in this file are intentionally small, header-only, and device callable.  They are used
 * by the continuum and cohesive-zone versions of the model so that the same reference-temperature,
 * crystallinity, hardening, and pressure-asymmetry functions are evaluated in both places.  Keeping
 * these functions in one file reduces the chance that the thin-film cohesive projection drifts away
 * from the bulk material model.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SURFACEINFORMEDPOLYMERHELPERS_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SURFACEINFORMEDPOLYMERHELPERS_HPP_

#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{
namespace surfaceInformedPolymerHelpers
{

/**
 * @brief Return a value bounded to a range that is safe for exp on host and device.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 boundedExpArgument( real64 const x )
{
  return LvArray::math::max( -50.0, LvArray::math::min( 50.0, x ) );
}

/**
 * @brief Logistic switch that approaches zero for negative arguments and one for positive arguments.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 logistic( real64 const x )
{
  return 1.0 / ( 1.0 + LvArray::math::exp( -boundedExpArgument( x ) ) );
}

/**
 * @brief Temperature scale normalized to one at the model transition/reference temperature.
 *
 * The scalar material parameters are interpreted as reference values at @p glassTransitionTemperature.
 * The scale is defined in logarithmic space to ensure positivity:
 *
 *   log S_T = a_cold max(Tg - T, 0) - a_hot max(T - Tg, 0)
 *             - J ( sigmoid((T-Tg)/w) - 1/2 ).
 *
 * The two slopes control the approximately linear behavior of log strength below and above the
 * reference temperature.  The last term allows a smooth drop centered on the reference temperature.
 * Setting @p temperatureTransitionMagnitude to zero reduces the expression to a two-sided Arrhenius-
 * like log-linear scale.  The expression is exactly normalized so that S_T(Tg)=1.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 temperatureScale( real64 const temperature,
                         real64 const glassTransitionTemperature,
                         real64 const temperatureColdSlope,
                         real64 const temperatureHotSlope,
                         real64 const temperatureTransitionMagnitude,
                         real64 const temperatureTransitionWidth )
{
  real64 const below = LvArray::math::max( glassTransitionTemperature - temperature, 0.0 );
  real64 const above = LvArray::math::max( temperature - glassTransitionTemperature, 0.0 );

  real64 logScale = temperatureColdSlope * below - temperatureHotSlope * above;
  if( temperatureTransitionWidth > 1.0e-12 && LvArray::math::abs( temperatureTransitionMagnitude ) > 1.0e-16 )
  {
    real64 const switchValue = logistic( ( temperature - glassTransitionTemperature ) / temperatureTransitionWidth );
    logScale -= temperatureTransitionMagnitude * ( switchValue - 0.5 );
  }

  return LvArray::math::exp( boundedExpArgument( logScale ) );
}

/**
 * @brief Smoothly turn on crystallinity sensitivity near the reference/transition temperature.
 *
 * Crystallinity has a weak effect far below the reference transition in many semicrystalline polymers,
 * where the amorphous phase is already glassy.  This activation keeps the correction optional and
 * bounded.  A non-positive width disables the switch and returns one so that a user may request a
 * temperature-independent crystallinity correction.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 crystallinityActivation( real64 const temperature,
                                real64 const glassTransitionTemperature,
                                real64 const crystallinityTransitionWidth )
{
  if( crystallinityTransitionWidth <= 1.0e-12 )
  {
    return 1.0;
  }
  return logistic( ( temperature - glassTransitionTemperature ) / crystallinityTransitionWidth );
}

/**
 * @brief Linear crystallinity multiplier bounded away from zero.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 crystallinityScale( real64 const temperature,
                           real64 const glassTransitionTemperature,
                           real64 const crystallinity,
                           real64 const referenceCrystallinity,
                           real64 const crystallinityCoeff,
                           real64 const crystallinityTransitionWidth )
{
  real64 const activation = crystallinityActivation( temperature, glassTransitionTemperature, crystallinityTransitionWidth );
  real64 const scale = 1.0 + crystallinityCoeff * ( crystallinity - referenceCrystallinity ) * activation;
  return LvArray::math::max( 1.0e-8, scale );
}

/**
 * @brief Temperature-localized pressure/tension-compression asymmetry coefficient.
 *
 * The pressure term is intended to modify the scalar yield envelope, not to introduce volumetric
 * plastic flow by itself.  The continuum update uses a non-associated radial return that keeps the
 * mean stress from the elastic trial state.  The cohesive-zone projection uses the same scalar term
 * with the film normal stress as the local one-dimensional mean-stress proxy.  The raw pressure
 * sensitivity is a low-pressure asymmetry and should be bounded before it is extrapolated into
 * strongly confined compression.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 pressureAsymmetry( real64 const temperature,
                          real64 const glassTransitionTemperature,
                          real64 const pressureAsymmetryAmplitude,
                          real64 const pressureAsymmetryWidth )
{
  if( pressureAsymmetryWidth <= 1.0e-12 || LvArray::math::abs( pressureAsymmetryAmplitude ) <= 1.0e-16 )
  {
    return 0.0;
  }
  real64 const x = ( temperature - glassTransitionTemperature ) / pressureAsymmetryWidth;
  return pressureAsymmetryAmplitude * LvArray::math::exp( boundedExpArgument( -0.5 * x * x ) );
}

/**
 * @brief Mean stress used by the pressure-asymmetry term after compressive limiting.
 *
 * The polymer pressure-asymmetry coefficient is calibrated from low-pressure tension/compression
 * differences.  If it is applied linearly at very large compressive pressure, it can create an
 * unbounded shear strength that is not supported by the calibration data.  This helper keeps the
 * tensile-positive mean stress unchanged in tension and clips only the compressive side.  A
 * negative @p compressivePressureStrengtheningCap recovers the uncapped model.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 pressureStrengtheningMeanStress( real64 const meanStress,
                                        real64 const compressivePressureStrengtheningCap )
{
  if( compressivePressureStrengtheningCap < 0.0 )
  {
    return meanStress;
  }
  return LvArray::math::max( meanStress, -compressivePressureStrengtheningCap );
}

/**
 * @brief Shear-softening contribution as a function of equivalent plastic strain.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 softeningContribution( real64 const shearSofteningMagnitude,
                              real64 const equivalentPlasticStrain,
                              real64 const shapeParameter1,
                              real64 const shapeParameter2 )
{
  real64 const positiveR1 = LvArray::math::max( shapeParameter1, 1.0e-16 );
  real64 const exponent = LvArray::math::pow( LvArray::math::max( equivalentPlasticStrain, 0.0 ) / positiveR1,
                                             LvArray::math::max( shapeParameter2, 1.0e-16 ) );
  return shearSofteningMagnitude * LvArray::math::exp( boundedExpArgument( -exponent ) );
}


/**
 * @brief Temperature-dependent fracture stretch for stretch-triggered damage.
 *
 * If @p fractureStretchLambda0 is non-positive, the legacy constant @p maximumStretch is returned.
 * Otherwise the damage trigger is
 *
 *   lambda_f = lambda_min + lambda_0 exp((T - T_0)/a).
 *
 * The exponent is numerically bounded for host/device robustness, but the material law itself does
 * not impose a high-temperature cap.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 fractureStretch( real64 const temperature,
                        real64 const maximumStretch,
                        real64 const fractureStretchLambdaMin,
                        real64 const fractureStretchLambda0,
                        real64 const fractureStretchT0,
                        real64 const fractureStretchTemperatureScale )
{
  if( fractureStretchLambda0 <= 0.0 || fractureStretchTemperatureScale <= 0.0 )
  {
    return maximumStretch;
  }
  real64 const exponent = ( temperature - fractureStretchT0 ) / fractureStretchTemperatureScale;
  return fractureStretchLambdaMin + fractureStretchLambda0 * LvArray::math::exp( boundedExpArgument( exponent ) );
}

/**
 * @brief Stretch-hardening driver used by both the continuum and cohesive-zone models.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 stretchHardeningMeasure( real64 const lambda )
{
  real64 const lambdaEff = LvArray::math::max( lambda, 1.0 );
  return lambdaEff * lambdaEff - 1.0 / lambdaEff;
}

/**
 * @brief Magnitude of a symmetric strain-like tensor stored in GEOS engineering Voigt form.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 strainMagnitude( real64 const ( &strain )[6] )
{
  real64 mag2 = 0.0;
  for( localIndex i = 0; i < 6; ++i )
  {
    // Normal components have one tensor entry.  Engineering shear components represent two symmetric
    // tensor entries after dividing by two, so their contribution is gamma^2/2.
    mag2 += 0.5 * ( 1 + ( i < 3 ) ) * strain[i] * strain[i];
  }
  return LvArray::math::sqrt( LvArray::math::max( mag2, 0.0 ) );
}

/**
 * @brief Maximum principal stretch from an unrotated deformation gradient.
 *
 * The right stretch tensor is computed from C = F^T F.  This avoids treating the non-symmetric
 * deformation gradient itself as a stretch tensor and makes the stretch calculation independent of
 * superposed rigid-body rotation.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 maximumPrincipalStretch( real64 const ( &unrotatedDeformationGradient )[3][3] )
{
  real64 C[3][3] = { { 0.0 } };
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      for( localIndex a = 0; a < 3; ++a )
      {
        C[i][j] += unrotatedDeformationGradient[a][i] * unrotatedDeformationGradient[a][j];
      }
    }
  }

  real64 Csym[6] = { 0.0 };
  LvArray::tensorOps::denseToSymmetric< 3 >( Csym, C );

  real64 principalStretchSquared[3] = { 0.0 };
  real64 eigenVectors[3][3] = { { 0.0 } };
  LvArray::tensorOps::symEigenvectors< 3 >( principalStretchSquared, eigenVectors, Csym );

  real64 maximumStretch = 0.0;
  for( localIndex i = 0; i < 3; ++i )
  {
    maximumStretch = LvArray::math::max( maximumStretch, LvArray::math::sqrt( LvArray::math::max( principalStretchSquared[i], 0.0 ) ) );
  }
  return maximumStretch;
}

/**
 * @brief Chain-stretch measure for hardening and finite-extensibility failure.
 *
 * The raw maximum principal stretch preserves the previous tensile calibration and captures ordinary
 * tensile extension.  The isochoric branch J^{-1/3} lambda_max suppresses hydrostatic compression
 * while still allowing constrained compression and shear under pressure to stretch chains.  The max
 * of the two measures keeps pure hydrostatic compression at unity, keeps tensile behavior unchanged
 * relative to the original maximum-principal-stretch driver, and makes finite-extensibility failure
 * pressure independent.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 chainStretch( real64 const ( &unrotatedDeformationGradient )[3][3] )
{
  real64 const lambdaMax = maximumPrincipalStretch( unrotatedDeformationGradient );
  real64 const J = LvArray::tensorOps::determinant< 3 >( unrotatedDeformationGradient );
  real64 const Jsafe = LvArray::math::max( J, 1.0e-16 );
  real64 const lambdaIsochoric = lambdaMax / LvArray::math::pow( Jsafe, 1.0 / 3.0 );
  return LvArray::math::max( lambdaMax, lambdaIsochoric );
}

} // namespace surfaceInformedPolymerHelpers
} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_SOLID_SURFACEINFORMEDPOLYMERHELPERS_HPP_
