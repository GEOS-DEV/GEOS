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

// Source includes
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
struct TestData : public TestDataBase< NC >
{
  real64 expectedZFactor; ///< Expected compressibility factor
};

template< TestFluidName FLUID, typename EOS_TYPE >
class CubicEOSPhaseModelTestFixture : public TestFluidFixture< FLUID, TestData >
{
public:
  using Base = TestFluidFixture< FLUID, TestData >;
  using Base::numComps;
  using Base::numDofs;
  using Base::relTol;
  static constexpr real64 absTol = 1.0e-4;
  using ParamType = TestData< numComps >;
  using EOS = CubicEOSPhaseModel< EOS_TYPE >;
  using Deriv = typename EOS::Deriv;

public:
  ~CubicEOSPhaseModelTestFixture() = default;

  void testPureCoefficients( ParamType const & testData );
  void testMixtureCoefficients( ParamType const & testData );
  void testATemperatureCoefficients( ParamType const & testData );
  void testCompressibilityFactor( ParamType const & testData );
  void testCompressibilityFactorValue( ParamType const & testData );
  void testLogFugacityCoefficients( ParamType const & testData );
  void testEnthalpy( ParamType const & testData );

  /// @return the step sizes and tolerances used for the derivative checks
  static DerivativeCheckSettings derivativeSettings()
  {
    DerivativeCheckSettings settings;
    settings.temperatureStep = 1.0e-4;
    settings.compositionStep = 1.0e-6;
    settings.absTolerance = absTol;
    settings.relTolerance = relTol;
    return settings;
  }
};

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testPureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;

  auto const binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSliceConst();
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  typename EOS::template StackVariables< true > stack( numComps, binaryInteractionCoefficients, derivs );
  EOS::template initialiseStack< true >(
    numComps,
    pressure,
    temperature,
    componentProperties,
    stack );

  // Values: the pure component coefficients a_i and b_i, interleaved
  integer constexpr numValues = 2*numComps;
  auto concatValues = []( auto const & a, auto const & b, auto & v ){
    for( integer ic = 0; ic < numComps; ic++ )
    {
      v[2*ic+0] = a[ic];
      v[2*ic+1] = b[ic];
    }
  };

  // The pure coefficients depend on pressure and temperature only
  stackArray2d< real64, numValues *numDofs > derivatives( numValues, numDofs );
  LvArray::forValuesInSlice( derivatives.toSlice(), []( real64 & a ){ a = 0.0; } );
  for( integer ic = 0; ic < numComps; ic++ )
  {
    derivatives( 2*ic+0, Deriv::dP ) = stack.daic_dp[ic];
    derivatives( 2*ic+1, Deriv::dP ) = stack.dbic_dp[ic];
    derivatives( 2*ic+0, Deriv::dT ) = stack.daic_dt[ic];
    derivatives( 2*ic+1, Deriv::dT ) = stack.dbic_dt[ic];
  }
  DerivativeCheckSettings settings = derivativeSettings();
  settings.pressureScale = 1.0e6;
  settings.checkComposition = false;
  Base::template checkNumericalDerivatives< numValues >( testData, derivatives.toSliceConst(),
                                                         [&]( real64 const p, real64 const t, auto const &, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, p, t, componentProperties, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values );
  }, settings );

  // Second order temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  stackArray1d< real64, numValues > secondDerivatives( numValues );
  concatValues( stack.d2aic_dt2, stack.d2bic_dt2, secondDerivatives );
  internal::testNumericalSecondDerivative< numValues >( temperature, dT, secondDerivatives.toSliceConst(), [&]( real64 const t, auto const & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values );
  }, absTol, relTol );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testMixtureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  auto const binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSliceConst();
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  typename EOS::template StackVariables< true > stack( numComps, binaryInteractionCoefficients, derivs );
  EOS::template initialiseStack< true >(
    numComps,
    testData.pressure,
    testData.temperature,
    componentProperties,
    stack );
  EOS::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

  // Values: the mixture coefficients a and b
  stackArray2d< real64, 2*numDofs > derivatives( 2, numDofs );
  for( integer idof = 0; idof < numDofs; ++idof )
  {
    derivatives( 0, idof ) = stack.daMixture[idof];
    derivatives( 1, idof ) = stack.dbMixture[idof];
  }
  DerivativeCheckSettings settings = derivativeSettings();
  settings.pressureScale = 1.0e6;
  Base::template checkNumericalDerivatives< 2 >( testData, derivatives.toSliceConst(),
                                                 [&]( real64 const p, real64 const t, auto const & z, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, p, t, componentProperties, valueStack );
    EOS::computeMixtureCoefficients( numComps, z, valueStack );
    values[0] = valueStack.aMixture;
    values[1] = valueStack.bMixture;
  }, settings );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testATemperatureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  auto const binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSliceConst();
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  typename EOS::template StackVariables< true > stack( numComps, binaryInteractionCoefficients, derivs );

  EOS::template initialiseStack< true >(
    numComps,
    testData.pressure,
    testData.temperature,
    componentProperties,
    stack );
  EOS::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

  // Value: the temperature derivative of the mixture coefficient a
  stackArray2d< real64, numDofs > derivatives( 1, numDofs );
  for( integer idof = 0; idof < numDofs; ++idof )
  {
    derivatives( 0, idof ) = stack.ddaMixture[idof];
  }
  DerivativeCheckSettings settings = derivativeSettings();
  settings.pressureScale = 1.0e6;
  Base::template checkNumericalDerivatives< 1 >( testData, derivatives.toSliceConst(),
                                                 [&]( real64 const p, real64 const t, auto const & z, auto & values )
  {
    typename EOS::template StackVariables< true > valueStack( numComps, binaryInteractionCoefficients, derivs );
    EOS::template initialiseStack< true >( numComps, p, t, componentProperties, valueStack );
    EOS::template computeMixtureCoefficients< 0, true >( numComps, z, valueStack );
    values[0] = valueStack.daMixture[Deriv::dT];
  }, settings );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testCompressibilityFactor( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  real64 compressibilityFactor = 0.0;
  stackArray2d< real64, numDofs > derivatives( 1, numDofs );
  EOS::computeCompressibilityFactorAndDerivs( numComps,
                                              testData.pressure,
                                              testData.temperature,
                                              composition.toSliceConst(),
                                              componentProperties,
                                              compressibilityFactor,
                                              derivatives[0] );

  DerivativeCheckSettings settings = derivativeSettings();
  settings.pressureScale = 1.0e6;
  Base::template checkNumericalDerivatives< 1 >( testData, derivatives.toSliceConst(),
                                                 [&]( real64 const p, real64 const t, auto const & z, auto & values )
  {
    EOS::computeCompressibilityFactor( numComps, p, t, z, componentProperties, values[0] );
  }, settings );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testLogFugacityCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  stackArray1d< real64, numComps > logFugacityCoefficients( numComps );
  stackArray2d< real64, numComps *numDofs > logFugacityCoefficientDerivs( numComps, numDofs );

  EOS::computeLogFugacityCoefficientsAndDerivs( numComps,
                                                testData.pressure,
                                                testData.temperature,
                                                composition.toSliceConst(),
                                                componentProperties,
                                                logFugacityCoefficients.toSlice(),
                                                logFugacityCoefficientDerivs.toSlice() );

  Base::template checkNumericalDerivatives< numComps >( testData, logFugacityCoefficientDerivs.toSliceConst(),
                                                        [&]( real64 const p, real64 const t, auto const & z, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps, p, t, z, componentProperties, values.toSlice() );
  }, derivativeSettings() );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testCompressibilityFactorValue( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );
  real64 const expectedZFactor = testData.expectedZFactor;

  real64 zFactor = 0.0;
  EOS::computeCompressibilityFactor( numComps,
                                     pressure,
                                     temperature,
                                     composition.toSliceConst(),
                                     componentProperties,
                                     zFactor );
  checkRelativeError( zFactor, expectedZFactor, relTol, absTol );
}

template< TestFluidName FLUID, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< FLUID, EOS_TYPE >::testEnthalpy( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, testData.composition );

  auto const binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff.toSliceConst();
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  typename EOS::template StackVariables< true > stack( numComps, binaryInteractionCoefficients, derivs );

  EOS::template initialiseStack< true >( numComps,
                                         testData.pressure,
                                         testData.temperature,
                                         componentProperties,
                                         stack );
  EOS::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

  real64 enthalpy = 0.0;
  stackArray2d< real64, numDofs > derivatives( 1, numDofs );
  EOS::template computeEnthalpy< 0, true >( numComps,
                                            testData.pressure,
                                            testData.temperature,
                                            composition.toSliceConst(),
                                            componentProperties,
                                            stack,
                                            enthalpy,
                                            derivatives[0] );

  DerivativeCheckSettings settings = derivativeSettings();
  settings.temperatureScale = 1.0e-2;
  Base::template checkNumericalDerivatives< 1 >( testData, derivatives.toSliceConst(),
                                                 [&]( real64 const p, real64 const t, auto const & z, auto & values )
  {
    typename EOS::template StackVariables< true > valueStack( numComps, binaryInteractionCoefficients, derivs );
    EOS::template initialiseStack< true >( numComps, p, t, componentProperties, valueStack );
    EOS::template computeMixtureCoefficients< 0, true >( numComps, z, valueStack );
    EOS::template computeEnthalpy< 0, false >( numComps, p, t, z, componentProperties, valueStack, values[0], nullptr );
  }, settings );
}

using PengRobinson4 = CubicEOSPhaseModelTestFixture< TestFluidName::AquiferGasStorage, PengRobinsonEOS >;
using SoaveRedlichKwong2 = CubicEOSPhaseModelTestFixture< TestFluidName::CarbonStorage, SoaveRedlichKwongEOS >;

TEST_P( PengRobinson4, testCubicModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testATemperatureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
  testEnthalpy( testParam );
}

TEST_P( SoaveRedlichKwong2, testCubicModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testATemperatureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
}

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P( CubicEOSPhaseModelTestFixture, PengRobinson4,
  ::testing::ValuesIn<TestData< 4 >>({
    {{1.000e+05, 297.15, {0.0309330, 0.3196830, 0.6378610, 0.0115230}}, 0.9958115},
    {{1.000e+05, 363.00, {0.0309330, 0.3196830, 0.6378610, 0.0115230}}, 0.9978260},
    {{1.840e+06, 297.15, {0.0309330, 0.3196830, 0.6378610, 0.0115230}}, 0.9211853},
    {{1.840e+06, 363.00, {0.0309330, 0.3196830, 0.6378610, 0.0115230}}, 0.9605058},
    {{1.840e+08, 297.15, {0.0309330, 0.3196830, 0.6378610, 0.0115230}}, 2.5422714},
    {{1.840e+08, 363.00, {0.0309330, 0.3196830, 0.6378610, 0.0115230}}, 2.2520420},
    {{1.000e+05, 297.15, {0.0000000, 0.3496860, 0.6378910, 0.0124230}}, 0.9957205},
    {{1.000e+05, 363.00, {0.0000000, 0.3496860, 0.6378910, 0.0124230}}, 0.9977697},
    {{1.840e+06, 297.15, {0.0000000, 0.3496860, 0.6378910, 0.0124230}}, 0.9193153},
    {{1.840e+06, 363.00, {0.0000000, 0.3496860, 0.6378910, 0.0124230}}, 0.9594314},
    {{1.840e+08, 297.15, {0.0000000, 0.3496860, 0.6378910, 0.0124230}}, 2.5442495},
    {{1.840e+08, 363.00, {0.0000000, 0.3496860, 0.6378910, 0.0124230}}, 2.2526541},
    {{1.000e+05, 297.15, {0.0000000, 0.3496860, 0.6503140, 0.0000000}}, 0.9957461},
    {{1.000e+05, 363.00, {0.0000000, 0.3496860, 0.6503140, 0.0000000}}, 0.9977879},
    {{1.840e+06, 297.15, {0.0000000, 0.3496860, 0.6503140, 0.0000000}}, 0.9198698},
    {{1.840e+06, 363.00, {0.0000000, 0.3496860, 0.6503140, 0.0000000}}, 0.9597917},
    {{1.840e+08, 297.15, {0.0000000, 0.3496860, 0.6503140, 0.0000000}}, 2.5542814},
    {{1.840e+08, 363.00, {0.0000000, 0.3496860, 0.6503140, 0.0000000}}, 2.2613856},
    {{1.000e+05, 297.15, {0.0000000, 0.0000000, 0.0000000, 1.0000000}}, 0.0008587},
    {{1.000e+05, 363.00, {0.0000000, 0.0000000, 0.0000000, 1.0000000}}, 0.0007388},
    {{1.840e+06, 297.15, {0.0000000, 0.0000000, 0.0000000, 1.0000000}}, 0.0157961},
    {{1.840e+06, 363.00, {0.0000000, 0.0000000, 0.0000000, 1.0000000}}, 0.0135883},
    {{1.840e+08, 297.15, {0.0000000, 0.0000000, 0.0000000, 1.0000000}}, 1.5518670},
    {{1.840e+08, 363.00, {0.0000000, 0.0000000, 0.0000000, 1.0000000}}, 1.3169829}
  })
);

INSTANTIATE_TEST_SUITE_P( CubicEOSPhaseModelTestFixture, SoaveRedlichKwong2,
  ::testing::ValuesIn<TestData< 2 >>({
    {{1.000e+05, 297.15, {0.9950000, 0.0050000}}, 0.0009671},
    {{1.000e+05, 363.00, {0.9950000, 0.0050000}}, 0.0008352},
    {{1.840e+06, 297.15, {0.9950000, 0.0050000}}, 0.0177905},
    {{1.840e+06, 363.00, {0.9950000, 0.0050000}}, 0.0153605},
    {{1.840e+08, 297.15, {0.9950000, 0.0050000}}, 1.7406285},
    {{1.840e+08, 363.00, {0.9950000, 0.0050000}}, 1.4790046},
    {{1.000e+05, 297.15, {0.9900000, 0.0100000}}, 0.0009700},
    {{1.000e+05, 363.00, {0.9900000, 0.0100000}}, 0.0008382},
    {{1.840e+06, 297.15, {0.9900000, 0.0100000}}, 0.0178437},
    {{1.840e+06, 363.00, {0.9900000, 0.0100000}}, 0.0154151},
    {{1.840e+08, 297.15, {0.9900000, 0.0100000}}, 1.7451857},
    {{1.840e+08, 363.00, {0.9900000, 0.0100000}}, 1.4832211},
    {{1.000e+05, 297.15, {0.1000000, 0.9000000}}, 0.9945206},
    {{1.000e+05, 363.00, {0.1000000, 0.9000000}}, 0.9972003},
    {{1.840e+06, 297.15, {0.1000000, 0.9000000}}, 0.8911050},
    {{1.840e+06, 363.00, {0.1000000, 0.9000000}}, 0.9473679},
    {{1.840e+08, 297.15, {0.1000000, 0.9000000}}, 2.6597629},
    {{1.840e+08, 363.00, {0.1000000, 0.9000000}}, 2.3395824}
  })
);

/* UNCRUSTIFY-ON */

} // namespace testing
} // namespace geos
