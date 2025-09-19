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
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
struct FluidData {};

template<>
struct FluidData< 2 >
{
  static std::unique_ptr< TestFluid< 2 > > create()
  {
    auto fluid = TestFluid< 2 >::create( {Fluid::H2O, Fluid::CO2} );
    fluid->setBinaryCoefficients( Feed< 1 >{ 0.1896 } );
    return fluid;
  }

  static std::array< Feed< 2 >, 3 > constexpr feeds = {
    Feed< 2 >{0.995, 0.005},
    Feed< 2 >{0.990, 0.010},
    Feed< 2 >{0.100, 0.900}
  };
};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > create()
  {
    auto fluid = TestFluid< 4 >::create( {Fluid::N2, Fluid::CH4, Fluid::CO2, Fluid::H2O} );
    fluid->setBinaryCoefficients( Feed< 6 >{ 0.0, 0.0, 0.0, 0.4778, 0.4850, 0.1896 } );
    return fluid;
  }

  static std::array< Feed< 4 >, 4 > constexpr feeds = {
    Feed< 4 >{0.030933, 0.319683, 0.637861, 0.011523},
    Feed< 4 >{0.000000, 0.349686, 0.637891, 0.012423},
    Feed< 4 >{0.000000, 0.349686, 0.650314, 0.000000},
    Feed< 4 >{0.000000, 0.000000, 0.000000, 1.000000}
  };
};

template< int NC >
using TestData = std::tuple<
  real64 const,         // Pressure
  real64 const,         // Temperature
  Feed< NC > const,     // Input composition
  real64 const          // Expected Z-factor
  >;

template< integer NC, typename EOS_TYPE >
class CubicEOSPhaseModelTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDofs = NC + 2;
  static constexpr real64 absTol = 1.0e-4;
  static constexpr real64 relTol = 1.0e-5;
  using ParamType = TestData< NC >;
  using EOS = CubicEOSPhaseModel< EOS_TYPE >;
  using Deriv = typename EOS::Deriv;

public:
  CubicEOSPhaseModelTestFixture():
    m_fluid( FluidData< NC >::create() )
  {}
  ~CubicEOSPhaseModelTestFixture() = default;

  void testPureCoefficients( ParamType const & testData );
  void testMixtureCoefficients( ParamType const & testData );
  void testCompressibilityFactor( ParamType const & testData );
  void testCompressibilityFactorValue( ParamType const & testData );
  void testLogFugacityCoefficients( ParamType const & testData );

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

template< integer NC, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< NC, EOS_TYPE >::testPureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );

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
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, p, temperature, componentProperties, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values, pressureScale );
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  concatValues( stack.daic_dt, stack.dbic_dt, derivatives );
  internal::testNumericalDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values );
  }, absTol, relTol );

  // Second order temperature derivatives
  concatValues( stack.d2aic_dt2, stack.d2bic_dt2, derivatives );
  internal::testNumericalSecondDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto const & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, valueStack );
    concatValues( valueStack.aic, valueStack.bic, values );
  }, absTol, relTol );
}

template< integer NC, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< NC, EOS_TYPE >::testMixtureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( testData ));

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
  EOS::template computeMixtureCoefficients< 0, true >( numComps, composition.toSliceConst(), stack );

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
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, p, temperature, componentProperties, valueStack );
    EOS::computeMixtureCoefficients( numComps, composition.toSliceConst(), valueStack );
    concatValues( valueStack, values, pressureScale );
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  concatDerivatives( stack, Deriv::dT, derivatives );
  internal::testNumericalDerivative< numValues >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto & values )
  {
    typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
    EOS::initialiseStack( numComps, pressure, t, componentProperties, valueStack );
    EOS::computeMixtureCoefficients( numComps, composition.toSliceConst(), valueStack );
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
      typename EOS::template StackVariables< false > valueStack( numComps, binaryInteractionCoefficients );
      EOS::initialiseStack( numComps, pressure, temperature, componentProperties, valueStack );
      EOS::computeMixtureCoefficients( numComps, composition.toSliceConst(), valueStack );
      concatValues( valueStack, values );
      composition[ic] = z_orig;
    }, absTol, relTol );
  }
}

template< integer NC, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< NC, EOS_TYPE >::testCompressibilityFactor( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( testData ));

  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, numDofs > compressibilityFactorDerivs( numDofs );

  EOS::computeCompressibilityFactorAndDerivs( numComps,
                                              pressure,
                                              temperature,
                                              composition.toSliceConst(),
                                              componentProperties,
                                              compressibilityFactor,
                                              compressibilityFactorDerivs.toSlice() );

  // Pressure derivative
  real64 constexpr pressureScale = 1.0e6;
  real64 const dp = 1.0e-4 * pressure;
  internal::testNumericalDerivative( pressure, dp, pressureScale*compressibilityFactorDerivs[Deriv::dP],
                                     [&]( real64 const p ) -> real64
  {
    real64 z_factor = 0.0;
    EOS::computeCompressibilityFactor( numComps,
                                       p,
                                       temperature,
                                       composition.toSliceConst(),
                                       componentProperties,
                                       z_factor );
    return pressureScale*z_factor;
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  internal::testNumericalDerivative( temperature, dT, compressibilityFactorDerivs[Deriv::dT],
                                     [&]( real64 const t ) -> real64
  {
    real64 z_factor = 0.0;
    EOS::computeCompressibilityFactor( numComps,
                                       pressure,
                                       t,
                                       composition.toSliceConst(),
                                       componentProperties,
                                       z_factor );
    return z_factor;
  }, absTol, relTol );

  // Composition derivatives
  real64 const dz = 1.0e-6;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    integer const idof = Deriv::dC + ic;
    internal::testNumericalDerivative( 0.0, dz, compressibilityFactorDerivs[idof],
                                       [&]( real64 const z ) -> real64
    {
      real64 z_factor = 0.0;
      real64 const z_orig = composition[ic];
      composition[ic] += z;
      EOS::computeCompressibilityFactor( numComps,
                                         pressure,
                                         temperature,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         z_factor );
      composition[ic] = z_orig;
      return z_factor;
    }, absTol, relTol );
  }
}

template< integer NC, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< NC, EOS_TYPE >::testLogFugacityCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( testData ));

  stackArray1d< real64, numComps > logFugacityCoefficients( numComps );
  stackArray2d< real64, numComps *numDofs > logFugacityCoefficientDerivs( numComps, numDofs );

  EOS::computeLogFugacityCoefficientsAndDerivs( numComps,
                                                pressure,
                                                temperature,
                                                composition.toSliceConst(),
                                                componentProperties,
                                                logFugacityCoefficients.toSlice(),
                                                logFugacityCoefficientDerivs.toSlice() );

  stackArray1d< real64, numComps > derivatives( numComps );

  auto concatDerivatives = [&]( integer const dof, auto & v ){
    for( integer ic = 0; ic < numComps; ic++ )
    {
      v[ic] = logFugacityCoefficientDerivs( ic, dof );
    }
  };

  // Pressure derivatives
  real64 const dp = 1.0e-4 * pressure;
  concatDerivatives( Deriv::dP, derivatives );
  internal::testNumericalDerivative< numComps >( pressure, dp, derivatives.toSliceConst(), [&]( real64 const p, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps,
                                         p,
                                         temperature,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         values.toSlice() );
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-4 * temperature;
  concatDerivatives( Deriv::dT, derivatives );
  internal::testNumericalDerivative< numComps >( temperature, dT, derivatives.toSliceConst(), [&]( real64 const t, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps,
                                         pressure,
                                         t,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         values.toSlice() );
  }, absTol, relTol );

  // Composition derivatives
  real64 const dz = 1.0e-6;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    concatDerivatives( Deriv::dC + ic, derivatives );
    internal::testNumericalDerivative< numComps >( 0, dz, derivatives.toSliceConst(), [&]( real64 const z, auto & values )
    {
      real64 const z_orig = composition[ic];
      composition[ic] += z;
      EOS::computeLogFugacityCoefficients( numComps,
                                           pressure,
                                           temperature,
                                           composition.toSliceConst(),
                                           componentProperties,
                                           values.toSlice() );
      composition[ic] = z_orig;
    }, absTol, relTol );
  }
}

template< integer NC, typename EOS_TYPE >
void
CubicEOSPhaseModelTestFixture< NC, EOS_TYPE >::testCompressibilityFactorValue( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( testData ));
  real64 const expectedZFactor = std::get< 3 >( testData );

  real64 zFactor = 0.0;
  EOS::computeCompressibilityFactor( numComps,
                                     pressure,
                                     temperature,
                                     composition.toSliceConst(),
                                     componentProperties,
                                     zFactor );
  checkRelativeError( zFactor, expectedZFactor, relTol, absTol );
}

using PengRobinson4 = CubicEOSPhaseModelTestFixture< 4, PengRobinsonEOS >;
using SoaveRedlichKwong2 = CubicEOSPhaseModelTestFixture< 2, SoaveRedlichKwongEOS >;

TEST_P( PengRobinson4, testCubicModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
}

TEST_P( SoaveRedlichKwong2, testCubicModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testCompressibilityFactorValue( testParam );
  testLogFugacityCoefficients( testParam );
}

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P( CubicEOSPhaseModelTestFixture, PengRobinson4,
  ::testing::ValuesIn<TestData< 4 >>({
    {1.000e+05, 297.15, { 0.0309330, 0.3196830, 0.6378610, 0.0115230 }, 0.9958115},
    {1.000e+05, 363.00, { 0.0309330, 0.3196830, 0.6378610, 0.0115230 }, 0.9978260},
    {1.840e+06, 297.15, { 0.0309330, 0.3196830, 0.6378610, 0.0115230 }, 0.9211853},
    {1.840e+06, 363.00, { 0.0309330, 0.3196830, 0.6378610, 0.0115230 }, 0.9605058},
    {1.840e+08, 297.15, { 0.0309330, 0.3196830, 0.6378610, 0.0115230 }, 2.5422714},
    {1.840e+08, 363.00, { 0.0309330, 0.3196830, 0.6378610, 0.0115230 }, 2.2520420},
    {1.000e+05, 297.15, { 0.0000000, 0.3496860, 0.6378910, 0.0124230 }, 0.9957205},
    {1.000e+05, 363.00, { 0.0000000, 0.3496860, 0.6378910, 0.0124230 }, 0.9977697},
    {1.840e+06, 297.15, { 0.0000000, 0.3496860, 0.6378910, 0.0124230 }, 0.9193153},
    {1.840e+06, 363.00, { 0.0000000, 0.3496860, 0.6378910, 0.0124230 }, 0.9594314},
    {1.840e+08, 297.15, { 0.0000000, 0.3496860, 0.6378910, 0.0124230 }, 2.5442495},
    {1.840e+08, 363.00, { 0.0000000, 0.3496860, 0.6378910, 0.0124230 }, 2.2526541},
    {1.000e+05, 297.15, { 0.0000000, 0.3496860, 0.6503140, 0.0000000 }, 0.9957461},
    {1.000e+05, 363.00, { 0.0000000, 0.3496860, 0.6503140, 0.0000000 }, 0.9977879},
    {1.840e+06, 297.15, { 0.0000000, 0.3496860, 0.6503140, 0.0000000 }, 0.9198698},
    {1.840e+06, 363.00, { 0.0000000, 0.3496860, 0.6503140, 0.0000000 }, 0.9597917},
    {1.840e+08, 297.15, { 0.0000000, 0.3496860, 0.6503140, 0.0000000 }, 2.5542814},
    {1.840e+08, 363.00, { 0.0000000, 0.3496860, 0.6503140, 0.0000000 }, 2.2613856},
    {1.000e+05, 297.15, { 0.0000000, 0.0000000, 0.0000000, 1.0000000 }, 0.0008587},
    {1.000e+05, 363.00, { 0.0000000, 0.0000000, 0.0000000, 1.0000000 }, 0.0007388},
    {1.840e+06, 297.15, { 0.0000000, 0.0000000, 0.0000000, 1.0000000 }, 0.0157961},
    {1.840e+06, 363.00, { 0.0000000, 0.0000000, 0.0000000, 1.0000000 }, 0.0135883},
    {1.840e+08, 297.15, { 0.0000000, 0.0000000, 0.0000000, 1.0000000 }, 1.5518670},
    {1.840e+08, 363.00, { 0.0000000, 0.0000000, 0.0000000, 1.0000000 }, 1.3169829}
  })
);

INSTANTIATE_TEST_SUITE_P( CubicEOSPhaseModelTestFixture, SoaveRedlichKwong2,
  ::testing::ValuesIn<TestData< 2 >>({
    {1.000e+05, 297.15, { 0.9950000, 0.0050000 }, 0.0009671},
    {1.000e+05, 363.00, { 0.9950000, 0.0050000 }, 0.0008352},
    {1.840e+06, 297.15, { 0.9950000, 0.0050000 }, 0.0177905},
    {1.840e+06, 363.00, { 0.9950000, 0.0050000 }, 0.0153605},
    {1.840e+08, 297.15, { 0.9950000, 0.0050000 }, 1.7406285},
    {1.840e+08, 363.00, { 0.9950000, 0.0050000 }, 1.4790046},
    {1.000e+05, 297.15, { 0.9900000, 0.0100000 }, 0.0009700},
    {1.000e+05, 363.00, { 0.9900000, 0.0100000 }, 0.0008382},
    {1.840e+06, 297.15, { 0.9900000, 0.0100000 }, 0.0178437},
    {1.840e+06, 363.00, { 0.9900000, 0.0100000 }, 0.0154151},
    {1.840e+08, 297.15, { 0.9900000, 0.0100000 }, 1.7451857},
    {1.840e+08, 363.00, { 0.9900000, 0.0100000 }, 1.4832211},
    {1.000e+05, 297.15, { 0.1000000, 0.9000000 }, 0.9945206},
    {1.000e+05, 363.00, { 0.1000000, 0.9000000 }, 0.9972003},
    {1.840e+06, 297.15, { 0.1000000, 0.9000000 }, 0.8911050},
    {1.840e+06, 363.00, { 0.1000000, 0.9000000 }, 0.9473679},
    {1.840e+08, 297.15, { 0.1000000, 0.9000000 }, 2.6597629},
    {1.840e+08, 363.00, { 0.1000000, 0.9000000 }, 2.3395824}
  })
);

/* UNCRUSTIFY-ON */

} // namespace testing
} // namespace geos
