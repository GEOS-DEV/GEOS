/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "constitutive/fluid/multifluid/compositional/functions/SoreideWhitsonEOSPhaseModel.hpp"
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"

namespace geos
{

using namespace constitutive;
using namespace constitutive::compositional;

namespace testing
{

template< integer NC >
struct TestData : public TestDataBase< NC >
{
  real64 salinity;        ///< Salinity
  real64 expectedZFactor; ///< Expected compressibility factor
};

template< TestFluidName FLUID, typename EOS_TYPE >
class SoreideWhitsonEOSPhaseModelTestFixture : public TestFluidFixture< FLUID, TestData >
{
public:
  using Base = TestFluidFixture< FLUID, TestData >;
  using Base::numComps;
  using Base::numDofs;
  using Base::relTol;
  static constexpr real64 absTol = 1.0e-4;
  using ParamType = TestData< numComps >;
  using EOS = SoreideWhitsonEOSPhaseModel< EOS_TYPE >;
  using CubicModel = typename EOS::CubicModel;
  using Deriv = typename EOS::Deriv;
public:
  ~SoreideWhitsonEOSPhaseModelTestFixture() = default;

  void testPureCoefficients( ParamType const & testData );
  void testBinaryInteractionCoefficients( ParamType const & testData );
  void testMixtureCoefficients( ParamType const & testData );
  void testATemperatureCoefficients( ParamType const & testData );
  void testCompressibilityFactor( ParamType const & testData );
  void testCompressibilityFactorValue( ParamType const & testData );
  void testLogFugacityCoefficients( ParamType const & testData );
  void testEnthalpy( ParamType const & testData );
};

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testPureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;

  typename EOS::template StackVariables< true > stack( numComps );
  EOS::initialiseStack( numComps,
                        pressure,
                        temperature,
                        componentProperties,
                        salinity,
                        stack );

  integer constexpr numValues = 2*numComps;
  stackArray1d< real64, numValues > derivatives( numValues );

  auto concatValues = []( auto const & a, auto const & b, auto & v, real64 const scale = 1.0 ){
    for( integer ic = 0; ic < numComps; ic++ )
    {
      v[2*ic+0] = scale*a[ic];
      v[2*ic+1] = scale*b[ic];
    }
  };

  // Pressure derivatives
  real64 constexpr pressureScale = 1.0e6;
  real64 const dp = 1.0e-4 * pressure;
  concatValues( stack.daic_dp, stack.dbic_dp, derivatives, pressureScale );
  internal::testNumericalDerivative< numValues >( pressure, dp, derivatives.toSliceConst(), [&]( real64 const p, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps );
    EOS::initialiseStack( numComps, p, temperature, componentProperties, salinity, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values, pressureScale );
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  concatValues( stack.daic_dt, stack.dbic_dt, derivatives );
  internal::testNumericalDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, salinity, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values );
  }, absTol, relTol );

  // Second order temperature derivatives
  concatValues( stack.d2aic_dt2, stack.d2bic_dt2, derivatives );
  internal::testNumericalSecondDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto const & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, salinity, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values );
  }, absTol, relTol );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testBinaryInteractionCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;

  // Test symmetry
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = ic; jc < numComps; ++jc )
    {
      real64 kij = 0.0;
      real64 kji = 0.0;
      real64 dk_dT = 0.0;
      EOS::getBinaryInteractionCoefficient( pressure,
                                            temperature,
                                            componentProperties,
                                            salinity,
                                            ic,
                                            jc,
                                            kij,
                                            dk_dT );
      EOS::getBinaryInteractionCoefficient( pressure,
                                            temperature,
                                            componentProperties,
                                            salinity,
                                            ic,
                                            jc,
                                            kji,
                                            dk_dT );
      EXPECT_NEAR( kij, kji, absTol );
      if( ic == jc )
      {
        EXPECT_NEAR( kij, 0.0, absTol );
      }
    }
  }

  // Test numerical derivatives
  real64 const dT = 1.0e-6 * temperature;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 kij = 0.0;
      real64 dkij_dT = 0.0;
      EOS::getBinaryInteractionCoefficient( pressure,
                                            temperature,
                                            componentProperties,
                                            salinity,
                                            ic,
                                            jc,
                                            kij,
                                            dkij_dT );

      internal::testNumericalDerivative( temperature, dT, dkij_dT,
                                         [&]( real64 t ) -> real64 {
        real64 l_kij, l_dkij_dT;
        EOS::getBinaryInteractionCoefficient( pressure,
                                              t,
                                              componentProperties,
                                              salinity,
                                              ic,
                                              jc,
                                              l_kij,
                                              l_dkij_dT );
        return l_kij;
      }, absTol, relTol );
    }
  }
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testMixtureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  typename EOS::template StackVariables< true > stack( numComps );
  EOS::initialiseStack( numComps,
                        pressure,
                        temperature,
                        componentProperties,
                        salinity,
                        stack );
  CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

  integer constexpr numValues = 2;
  stackArray1d< real64, numValues > derivatives( numValues );

  auto concatValues = []( auto const & s, auto & v, real64 const scale = 1.0 ){
    v[0] = scale*s.aMixture;
    v[1] = scale*s.bMixture;
  };
  auto concatDerivatives = []( auto const & s, integer const dof, auto & v, real64 const scale = 1.0 ){
    v[0] = scale*s.daMixture[dof];
    v[1] = scale*s.dbMixture[dof];
  };

  // Pressure derivatives
  real64 constexpr pressureScale = 1.0e6;
  real64 const dp = 1.0e-4 * pressure;
  concatDerivatives( stack, Deriv::dP, derivatives, pressureScale );
  internal::testNumericalDerivative< numValues >( pressure, dp, derivatives.toSliceConst(), [&]( real64 const p, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps );
    EOS::initialiseStack( numComps, p, temperature, componentProperties, salinity, valueStack );
    CubicModel::computeMixtureCoefficients( numComps, composition.toSliceConst(), valueStack );
    concatValues( valueStack, values, pressureScale );
  }, absTol, relTol );

  // Temperature derivatives
  real64 constexpr temperatureScale = 1.0e4;
  real64 const dT = 1.0e-4 * temperature;
  concatDerivatives( stack, Deriv::dT, derivatives, temperatureScale );
  internal::testNumericalDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, salinity, valueStack );
    CubicModel::computeMixtureCoefficients( numComps, composition.toSliceConst(), valueStack );
    concatValues( valueStack, values, temperatureScale );
  }, absTol, relTol );

  // Composition derivatives
  real64 const dz = 1.0e-6;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const idof = Deriv::dC + ic;
    concatDerivatives( stack, idof, derivatives );
    internal::testNumericalDerivative< numValues >( 0, dz, derivatives.toSliceConst(), [&]( real64 const z, auto & values )
    {
      real64 const z_orig = composition[ic];
      composition[ic] += z;
      typename EOS::template StackVariables< false > valueStack( numComps );
      EOS::initialiseStack( numComps, pressure, temperature, componentProperties, salinity, valueStack );
      CubicModel::computeMixtureCoefficients( numComps, composition.toSliceConst(), valueStack );
      concatValues( valueStack, values );
      composition[ic] = z_orig;
    }, absTol, relTol );
  }
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testATemperatureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  typename EOS::template StackVariables< true > stack( numComps );

  EOS::initialiseStack( numComps,
                        pressure,
                        temperature,
                        componentProperties,
                        salinity,
                        stack );
  CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

  integer constexpr numValues = 1;
  stackArray1d< real64, numValues > derivatives( numValues );

  auto concatValues = []( auto const & s, auto & v, real64 const scale = 1.0 ){
    v[0] = scale*s.daMixture[Deriv::dT];
  };
  auto concatDerivatives = []( auto const & s, integer const dof, auto & v, real64 const scale = 1.0 ){
    v[0] = scale*s.ddaMixture[dof];
  };

  // Pressure derivatives
  real64 const dp = 1.0e-4 * pressure;
  concatDerivatives( stack, Deriv::dP, derivatives );
  internal::testNumericalDerivative< numValues >( pressure, dp, derivatives.toSliceConst(), [&]( real64 const p, auto & values )
  {
    typename EOS::template StackVariables< true > valueStack( numComps );
    EOS::initialiseStack( numComps, p, temperature, componentProperties, salinity, valueStack );
    CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), valueStack );
    concatValues( valueStack, values );
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  concatDerivatives( stack, Deriv::dT, derivatives );
  internal::testNumericalDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto & values )
  {
    typename EOS::template StackVariables< true > valueStack( numComps );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, salinity, valueStack );
    CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), valueStack );
    concatValues( valueStack, values );
  }, absTol, relTol );

  // Composition derivatives
  real64 const dz = 1.0e-6;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const idof = Deriv::dC + ic;
    concatDerivatives( stack, idof, derivatives );
    internal::testNumericalDerivative< numValues >( 0, dz, derivatives.toSliceConst(), [&]( real64 const z, auto & values )
    {
      real64 const z_orig = composition[ic];
      composition[ic] += z;
      typename EOS::template StackVariables< true > valueStack( numComps );
      EOS::initialiseStack( numComps, pressure, temperature, componentProperties, salinity, valueStack );
      CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), valueStack );
      concatValues( valueStack, values );
      composition[ic] = z_orig;
    }, absTol, relTol );
  }
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testCompressibilityFactor( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, numDofs > compressibilityFactorDerivs( numDofs );

  EOS::computeCompressibilityFactorAndDerivs( numComps,
                                              pressure,
                                              temperature,
                                              composition.toSliceConst(),
                                              componentProperties,
                                              salinity,
                                              compressibilityFactor,
                                              compressibilityFactorDerivs.toSlice() );
  // Pressure derivative
  real64 constexpr pressureScale = 1.0e6;
  real64 const dp = 1.0e-4 * pressure;
  internal::testNumericalDerivative( pressure, dp, pressureScale*compressibilityFactorDerivs[Deriv::dP],
                                     [&]( real64 const p ) -> real64
  {
    real64 zfactor = 0.0;
    EOS::computeCompressibilityFactor( numComps, p, temperature, composition.toSliceConst(), componentProperties, salinity, zfactor );
    return pressureScale*zfactor;
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  internal::testNumericalDerivative( temperature, dT, compressibilityFactorDerivs[Deriv::dT],
                                     [&]( real64 const t ) -> real64
  {
    real64 zfactor = 0.0;
    EOS::computeCompressibilityFactor( numComps, pressure, t, composition.toSliceConst(), componentProperties, salinity, zfactor );
    return zfactor;
  }, absTol, relTol );

  real64 const dz = 1.0e-7;
  for( integer kc = 0; kc < numComps; kc++ )
  {
    internal::testNumericalDerivative( 0.0, dz, compressibilityFactorDerivs[Deriv::dC+kc],
                                       [&]( real64 const z ) -> real64
    {
      real64 const z_old = composition[kc];
      composition[kc] += z;
      real64 zfactor = 0.0;
      EOS::computeCompressibilityFactor( numComps, pressure, temperature, composition.toSliceConst(), componentProperties, salinity, zfactor );
      composition[kc] = z_old;
      return zfactor;
    }, absTol, relTol );
  }
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testCompressibilityFactorValue( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );
  real64 const expectedZFactor = testData.expectedZFactor;

  real64 zFactor = 0.0;
  EOS::computeCompressibilityFactor( numComps,
                                     pressure,
                                     temperature,
                                     composition.toSliceConst(),
                                     componentProperties,
                                     salinity,
                                     zFactor );
  checkRelativeError( zFactor, expectedZFactor, relTol, absTol );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testLogFugacityCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;
  StackArray< real64, 1, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  StackArray< real64, 1, numComps > logFugacityCoefficients( numComps );
  StackArray< real64, 2, numComps *numDofs > logFugacityCoefficientDerivs( numComps, numDofs );

  // Inflate the values of the log fugacity coefficients to catch errors
  real64 constexpr logFugacityScale = 1.0e3;
  StackArray< real64, 1, numComps > derivatives( numComps );

  EOS::computeLogFugacityCoefficientsAndDerivs( numComps,
                                                pressure,
                                                temperature,
                                                composition.toSliceConst(),
                                                componentProperties,
                                                salinity,
                                                logFugacityCoefficients.toSlice(),
                                                logFugacityCoefficientDerivs.toSlice() );

  auto const concatDerivatives = [&]( integer idof, real64 const scale )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      derivatives[ic] = scale * logFugacityCoefficientDerivs( ic, idof );
    }
  };

  // Pressure derivatives
  real64 constexpr pressureScale = 1.0e6;
  concatDerivatives( Deriv::dP, pressureScale*logFugacityScale );
  real64 const dp = 1.0e-4 * pressure;
  internal::testNumericalDerivative< numComps >( pressure, dp, derivatives.toSliceConst(),
                                                 [&]( real64 const p, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps,
                                         p,
                                         temperature,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         salinity,
                                         values );
    LvArray::forValuesInSlice( values.toSlice(), [=]( real64 & a ){ a *= pressureScale*logFugacityScale; } );
  }, absTol, relTol );

  // Temperature derivatives
  concatDerivatives( Deriv::dT, logFugacityScale );
  real64 const dT = 1.0e-6 * temperature;
  internal::testNumericalDerivative< numComps >( temperature, dT, derivatives.toSliceConst(),
                                                 [&]( real64 const t, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps,
                                         pressure,
                                         t,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         salinity,
                                         values );
    LvArray::forValuesInSlice( values.toSlice(), [=]( real64 & a ){ a *= logFugacityScale; } );
  }, absTol, relTol );

  // Composition derivatives
  real64 constexpr dz = 1.0e-6;
  for( integer jc = 0; jc < numComps; ++jc )
  {
    concatDerivatives( Deriv::dC + jc, logFugacityScale );
    internal::testNumericalDerivative< numComps >( 0.0, dz, derivatives.toSliceConst(),
                                                   [&]( real64 const z, auto & values )
    {
      real64 const zj_old = composition[jc];
      composition[jc] += z;
      EOS::computeLogFugacityCoefficients( numComps,
                                           pressure,
                                           temperature,
                                           composition.toSliceConst(),
                                           componentProperties,
                                           salinity,
                                           values );
      composition[jc] = zj_old;
      LvArray::forValuesInSlice( values.toSlice(), [=]( real64 & a ){ a *= logFugacityScale; } );
    }, absTol, relTol );
  }
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testEnthalpy( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  real64 const salinity = testData.salinity;
  StackArray< real64, 1, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  constexpr real64 enthalpyScale = 5.0e-4;

  real64 enthalpy = 0.0;
  stackArray1d< real64, numDofs > enthalpyDerivs( numDofs );

  typename EOS::template StackVariables< true > stack( numComps );

  EOS::initialiseStack( numComps,
                        pressure,
                        temperature,
                        componentProperties,
                        salinity,
                        stack );
  CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

  CubicModel::template computeEnthalpy< 0, true >( numComps,
                                                   pressure,
                                                   temperature,
                                                   composition.toSliceConst(),
                                                   componentProperties,
                                                   stack,
                                                   enthalpy,
                                                   enthalpyDerivs.toSlice() );

  // Pressure derivative
  real64 const dp = 1.0e-4 * pressure;
  internal::testNumericalDerivative( pressure, dp, enthalpyDerivs[Deriv::dP],
                                     [&]( real64 const p ) -> real64
  {
    typename EOS::template StackVariables< true > valueStack( numComps );
    EOS::initialiseStack( numComps, p, temperature, componentProperties, salinity, valueStack );
    CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), valueStack );

    real64 h = 0.0;
    CubicModel::template computeEnthalpy< 0, false >( numComps,
                                                      p,
                                                      temperature,
                                                      composition.toSliceConst(),
                                                      componentProperties,
                                                      valueStack,
                                                      h,
                                                      nullptr );
    return h;
  }, absTol, relTol );

  // Temperature derivative
  real64 const dT = 1.0e-4 * temperature;
  internal::testNumericalDerivative( temperature, dT, enthalpyScale * enthalpyDerivs[Deriv::dT],
                                     [&]( real64 const t ) -> real64
  {
    typename EOS::template StackVariables< true > valueStack( numComps );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, salinity, valueStack );
    CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), valueStack );

    real64 h = 0.0;
    CubicModel::template computeEnthalpy< 0, false >( numComps,
                                                      pressure,
                                                      t,
                                                      composition.toSliceConst(),
                                                      componentProperties,
                                                      valueStack,
                                                      h,
                                                      nullptr );
    return enthalpyScale * h;
  }, 10*absTol, 10*relTol );

  // Composition derivatives
  real64 const dz = 1.0e-7;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const idof = Deriv::dC + ic;
    internal::testNumericalDerivative( 0.0, dz, enthalpyScale * enthalpyDerivs[idof],
                                       [&]( real64 const z ) -> real64
    {
      real64 const z_orig = composition[ic];
      composition[ic] += z;

      typename EOS::template StackVariables< true > valueStack( numComps );
      EOS::initialiseStack( numComps, pressure, temperature, componentProperties, salinity, valueStack );
      CubicModel::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), valueStack );

      real64 h = 0.0;
      CubicModel::template computeEnthalpy< 0, false >( numComps,
                                                        pressure,
                                                        temperature,
                                                        composition.toSliceConst(),
                                                        componentProperties,
                                                        valueStack,
                                                        h,
                                                        nullptr );
      composition[ic] = z_orig;
      return enthalpyScale * h;
    }, 10*absTol, relTol );
  }
}

using PengRobinson2 = SoreideWhitsonEOSPhaseModelTestFixture< TestFluidName::OilWater, PengRobinsonEOS >;
using PengRobinson4 = SoreideWhitsonEOSPhaseModelTestFixture< TestFluidName::AquiferGasStorage, PengRobinsonEOS >;
using SoaveRedlichKwong3 = SoreideWhitsonEOSPhaseModelTestFixture< TestFluidName::HydrogenStorage, SoaveRedlichKwongEOS >;

TEST_P( PengRobinson2, testSWModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testBinaryInteractionCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testATemperatureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
  testEnthalpy( testParam );
}

TEST_P( PengRobinson4, testSWModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testBinaryInteractionCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
}

TEST_P( SoaveRedlichKwong3, testSWModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testBinaryInteractionCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
}

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSPhaseModelTest, PengRobinson2,
  ::testing::ValuesIn<TestData< 2 >>({
    {{1.00e+05, 297.15, {0.995000, 0.005000}}, 0.0,  0.00090026},
    {{1.00e+05, 297.15, {0.995000, 0.005000}}, 1.7,  0.00089960},
    {{1.00e+05, 363.00, {0.995000, 0.005000}}, 0.0,  0.00077645},
    {{1.00e+05, 363.00, {0.995000, 0.005000}}, 1.7,  0.00077531},
    {{1.84e+06, 297.15, {0.995000, 0.005000}}, 0.0,  0.01656101},
    {{1.84e+06, 297.15, {0.995000, 0.005000}}, 1.7,  0.01654886},
    {{1.84e+06, 363.00, {0.995000, 0.005000}}, 0.0,  0.01428025},
    {{1.84e+06, 363.00, {0.995000, 0.005000}}, 1.7,  0.01425936},
    {{1.84e+08, 297.15, {0.995000, 0.005000}}, 0.0,  1.62425580},
    {{1.84e+08, 297.15, {0.995000, 0.005000}}, 1.7,  1.62345870},
    {{1.84e+08, 363.00, {0.995000, 0.005000}}, 0.0,  1.37947706},
    {{1.84e+08, 363.00, {0.995000, 0.005000}}, 1.7,  1.37827630},
    {{1.00e+05, 297.15, {1.000000, 0.000000}}, 0.0,  0.00085955},
    {{1.00e+05, 297.15, {1.000000, 0.000000}}, 1.7,  0.00085886},
    {{1.00e+05, 363.00, {1.000000, 0.000000}}, 0.0,  0.00073994},
    {{1.00e+05, 363.00, {1.000000, 0.000000}}, 1.7,  0.00073876},
    {{1.84e+06, 297.15, {1.000000, 0.000000}}, 0.0,  0.01581228},
    {{1.84e+06, 297.15, {1.000000, 0.000000}}, 1.7,  0.01579967},
    {{1.84e+06, 363.00, {1.000000, 0.000000}}, 0.0,  0.01360943},
    {{1.84e+06, 363.00, {1.000000, 0.000000}}, 1.7,  0.01358772},
    {{1.84e+08, 297.15, {1.000000, 0.000000}}, 0.0,  1.55295404},
    {{1.84e+08, 297.15, {1.000000, 0.000000}}, 1.7,  1.55210598},
    {{1.84e+08, 363.00, {1.000000, 0.000000}}, 0.0,  1.31823992},
    {{1.84e+08, 363.00, {1.000000, 0.000000}}, 1.7,  1.31694721},
    {{1.00e+05, 297.15, {0.002000, 0.998000}}, 0.0,  0.00856008},
    {{1.00e+05, 297.15, {0.002000, 0.998000}}, 1.7,  0.00856009},
    {{1.00e+05, 363.00, {0.002000, 0.998000}}, 0.0,  0.00738480},
    {{1.00e+05, 363.00, {0.002000, 0.998000}}, 1.7,  0.00738484},
    {{1.84e+06, 297.15, {0.002000, 0.998000}}, 0.0,  0.15720092},
    {{1.84e+06, 297.15, {0.002000, 0.998000}}, 1.7,  0.15720122},
    {{1.84e+06, 363.00, {0.002000, 0.998000}}, 0.0,  0.13535973},
    {{1.84e+06, 363.00, {0.002000, 0.998000}}, 1.7,  0.13536038},
    {{1.84e+08, 297.15, {0.002000, 0.998000}}, 0.0, 14.71187824},
    {{1.84e+08, 297.15, {0.002000, 0.998000}}, 1.7, 14.71188210},
    {{1.84e+08, 363.00, {0.002000, 0.998000}}, 0.0, 12.18401104},
    {{1.84e+08, 363.00, {0.002000, 0.998000}}, 1.7, 12.18401699}
  })
);

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSPhaseModelTest, PengRobinson4,
  ::testing::ValuesIn<TestData< 4 >>({
    {{1.00e+05, 297.15, {0.030933, 0.319683, 0.637861, 0.011523}}, 0.0, 0.99572092},
    {{1.00e+05, 297.15, {0.030933, 0.319683, 0.637861, 0.011523}}, 1.7, 0.99572434},
    {{1.00e+05, 363.00, {0.030933, 0.319683, 0.637861, 0.011523}}, 0.0, 0.99778476},
    {{1.00e+05, 363.00, {0.030933, 0.319683, 0.637861, 0.011523}}, 1.7, 0.99778741},
    {{1.84e+06, 297.15, {0.030933, 0.319683, 0.637861, 0.011523}}, 0.0, 0.91930630},
    {{1.84e+06, 297.15, {0.030933, 0.319683, 0.637861, 0.011523}}, 1.7, 0.91937745},
    {{1.84e+06, 363.00, {0.030933, 0.319683, 0.637861, 0.011523}}, 0.0, 0.95971259},
    {{1.84e+06, 363.00, {0.030933, 0.319683, 0.637861, 0.011523}}, 1.7, 0.95976376},
    {{1.84e+08, 297.15, {0.030933, 0.319683, 0.637861, 0.011523}}, 0.0, 2.53720258},
    {{1.84e+08, 297.15, {0.030933, 0.319683, 0.637861, 0.011523}}, 1.7, 2.53739289},
    {{1.84e+08, 363.00, {0.030933, 0.319683, 0.637861, 0.011523}}, 0.0, 2.24823017},
    {{1.84e+08, 363.00, {0.030933, 0.319683, 0.637861, 0.011523}}, 1.7, 2.24847504},
    {{1.00e+05, 297.15, {0.000000, 0.349686, 0.637891, 0.012423}}, 0.0, 0.99562285},
    {{1.00e+05, 297.15, {0.000000, 0.349686, 0.637891, 0.012423}}, 1.7, 0.99562638},
    {{1.00e+05, 363.00, {0.000000, 0.349686, 0.637891, 0.012423}}, 0.0, 0.99772526},
    {{1.00e+05, 363.00, {0.000000, 0.349686, 0.637891, 0.012423}}, 1.7, 0.99772807},
    {{1.84e+06, 297.15, {0.000000, 0.349686, 0.637891, 0.012423}}, 0.0, 0.91728093},
    {{1.84e+06, 297.15, {0.000000, 0.349686, 0.637891, 0.012423}}, 1.7, 0.91735461},
    {{1.84e+06, 363.00, {0.000000, 0.349686, 0.637891, 0.012423}}, 0.0, 0.95857373},
    {{1.84e+06, 363.00, {0.000000, 0.349686, 0.637891, 0.012423}}, 1.7, 0.95862796},
    {{1.84e+08, 297.15, {0.000000, 0.349686, 0.637891, 0.012423}}, 0.0, 2.53887060},
    {{1.84e+08, 297.15, {0.000000, 0.349686, 0.637891, 0.012423}}, 1.7, 2.53906366},
    {{1.84e+08, 363.00, {0.000000, 0.349686, 0.637891, 0.012423}}, 0.0, 2.24859657},
    {{1.84e+08, 363.00, {0.000000, 0.349686, 0.637891, 0.012423}}, 1.7, 2.24885192},
    {{1.00e+05, 297.15, {0.000000, 0.349686, 0.650314, 0.000000}}, 0.0, 0.99574609},
    {{1.00e+05, 297.15, {0.000000, 0.349686, 0.650314, 0.000000}}, 1.7, 0.99574609},
    {{1.00e+05, 363.00, {0.000000, 0.349686, 0.650314, 0.000000}}, 0.0, 0.99778787},
    {{1.00e+05, 363.00, {0.000000, 0.349686, 0.650314, 0.000000}}, 1.7, 0.99778787},
    {{1.84e+06, 297.15, {0.000000, 0.349686, 0.650314, 0.000000}}, 0.0, 0.91986981},
    {{1.84e+06, 297.15, {0.000000, 0.349686, 0.650314, 0.000000}}, 1.7, 0.91986981},
    {{1.84e+06, 363.00, {0.000000, 0.349686, 0.650314, 0.000000}}, 0.0, 0.95979173},
    {{1.84e+06, 363.00, {0.000000, 0.349686, 0.650314, 0.000000}}, 1.7, 0.95979173},
    {{1.84e+08, 297.15, {0.000000, 0.349686, 0.650314, 0.000000}}, 0.0, 2.55428138},
    {{1.84e+08, 297.15, {0.000000, 0.349686, 0.650314, 0.000000}}, 1.7, 2.55428138},
    {{1.84e+08, 363.00, {0.000000, 0.349686, 0.650314, 0.000000}}, 0.0, 2.26138556},
    {{1.84e+08, 363.00, {0.000000, 0.349686, 0.650314, 0.000000}}, 1.7, 2.26138556},
    {{1.00e+05, 297.15, {0.000000, 0.000000, 0.000000, 1.000000}}, 0.0, 0.00085955},
    {{1.00e+05, 297.15, {0.000000, 0.000000, 0.000000, 1.000000}}, 1.7, 0.00085886},
    {{1.00e+05, 363.00, {0.000000, 0.000000, 0.000000, 1.000000}}, 0.0, 0.00073994},
    {{1.00e+05, 363.00, {0.000000, 0.000000, 0.000000, 1.000000}}, 1.7, 0.00073876},
    {{1.84e+06, 297.15, {0.000000, 0.000000, 0.000000, 1.000000}}, 0.0, 0.01581228},
    {{1.84e+06, 297.15, {0.000000, 0.000000, 0.000000, 1.000000}}, 1.7, 0.01579967},
    {{1.84e+06, 363.00, {0.000000, 0.000000, 0.000000, 1.000000}}, 0.0, 0.01360943},
    {{1.84e+06, 363.00, {0.000000, 0.000000, 0.000000, 1.000000}}, 1.7, 0.01358772},
    {{1.84e+08, 297.15, {0.000000, 0.000000, 0.000000, 1.000000}}, 0.0, 1.55295404},
    {{1.84e+08, 297.15, {0.000000, 0.000000, 0.000000, 1.000000}}, 1.7, 1.55210598},
    {{1.84e+08, 363.00, {0.000000, 0.000000, 0.000000, 1.000000}}, 0.0, 1.31823992},
    {{1.84e+08, 363.00, {0.000000, 0.000000, 0.000000, 1.000000}}, 1.7, 1.31694721}
  })
);

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSPhaseModelTest, SoaveRedlichKwong3,
  ::testing::ValuesIn<TestData< 3 >>({
    {{1.00e+05, 297.15, {0.995000, 0.000000, 0.005000}}, 0.0, 0.00097447},
    {{1.00e+05, 297.15, {0.995000, 0.000000, 0.005000}}, 1.7, 0.00097360},
    {{1.00e+05, 363.00, {0.995000, 0.000000, 0.005000}}, 0.0, 0.99162416},
    {{1.00e+05, 363.00, {0.995000, 0.000000, 0.005000}}, 1.7, 0.99155411},
    {{1.84e+06, 297.15, {0.995000, 0.000000, 0.005000}}, 0.0, 0.01792455},
    {{1.84e+06, 297.15, {0.995000, 0.000000, 0.005000}}, 1.7, 0.01790864},
    {{1.84e+06, 363.00, {0.995000, 0.000000, 0.005000}}, 0.0, 0.01552264},
    {{1.84e+06, 363.00, {0.995000, 0.000000, 0.005000}}, 1.7, 0.01549547},
    {{1.84e+08, 297.15, {0.995000, 0.000000, 0.005000}}, 0.0, 1.74740244},
    {{1.84e+08, 297.15, {0.995000, 0.000000, 0.005000}}, 1.7, 1.74642816},
    {{1.84e+08, 363.00, {0.995000, 0.000000, 0.005000}}, 0.0, 1.48626001},
    {{1.84e+08, 363.00, {0.995000, 0.000000, 0.005000}}, 1.7, 1.48481033},
    {{1.00e+05, 297.15, {0.990000, 0.005000, 0.005000}}, 0.0, 0.00097714},
    {{1.00e+05, 297.15, {0.990000, 0.005000, 0.005000}}, 1.7, 0.00097626},
    {{1.00e+05, 363.00, {0.990000, 0.005000, 0.005000}}, 0.0, 0.99165029},
    {{1.00e+05, 363.00, {0.990000, 0.005000, 0.005000}}, 1.7, 0.99158070},
    {{1.84e+06, 297.15, {0.990000, 0.005000, 0.005000}}, 0.0, 0.01797356},
    {{1.84e+06, 297.15, {0.990000, 0.005000, 0.005000}}, 1.7, 0.01795757},
    {{1.84e+06, 363.00, {0.990000, 0.005000, 0.005000}}, 0.0, 0.01557193},
    {{1.84e+06, 363.00, {0.990000, 0.005000, 0.005000}}, 1.7, 0.01554458},
    {{1.84e+08, 297.15, {0.990000, 0.005000, 0.005000}}, 0.0, 1.75167280},
    {{1.84e+08, 297.15, {0.990000, 0.005000, 0.005000}}, 1.7, 1.75069753},
    {{1.84e+08, 363.00, {0.990000, 0.005000, 0.005000}}, 0.0, 1.49014828},
    {{1.84e+08, 363.00, {0.990000, 0.005000, 0.005000}}, 1.7, 1.48869666},
    {{1.00e+05, 297.15, {0.970000, 0.025000, 0.005000}}, 0.0, 0.00098788},
    {{1.00e+05, 297.15, {0.970000, 0.025000, 0.005000}}, 1.7, 0.00098699},
    {{1.00e+05, 363.00, {0.970000, 0.025000, 0.005000}}, 0.0, 0.99175435},
    {{1.00e+05, 363.00, {0.970000, 0.025000, 0.005000}}, 1.7, 0.99168656},
    {{1.84e+06, 297.15, {0.970000, 0.025000, 0.005000}}, 0.0, 0.01817086},
    {{1.84e+06, 297.15, {0.970000, 0.025000, 0.005000}}, 1.7, 0.01815457},
    {{1.84e+06, 363.00, {0.970000, 0.025000, 0.005000}}, 0.0, 0.01577122},
    {{1.84e+06, 363.00, {0.970000, 0.025000, 0.005000}}, 1.7, 0.01574312},
    {{1.84e+08, 297.15, {0.970000, 0.025000, 0.005000}}, 0.0, 1.76879530},
    {{1.84e+08, 297.15, {0.970000, 0.025000, 0.005000}}, 1.7, 1.76781650},
    {{1.84e+08, 363.00, {0.970000, 0.025000, 0.005000}}, 0.0, 1.50574959},
    {{1.84e+08, 363.00, {0.970000, 0.025000, 0.005000}}, 1.7, 1.50429088}
  })
);

/* UNCRUSTIFY-ON */

} // namespace testing
} // namespace geos
