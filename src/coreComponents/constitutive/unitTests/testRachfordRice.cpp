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


// Source includes
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

static constexpr real64 relTol = 1e-4;

TEST( RachfordRiceTest, testRachfordRiceTwoComponents )
{
  constexpr integer numComps = 2;

  array1d< real64 > kValues( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  ////////////////////////////////////////

  kValues[0] = 1.12223;
  kValues[1] = 1.12223;

  feed[0] = 0.5;
  feed[1] = 0.5;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction1 = 1;
  real64 const vaporFraction1 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction1, expectedVaporFraction1, relTol );

  ////////////////////////////////////////

  kValues[0] = 56.1091;
  kValues[1] = 56.1091;

  feed[0] = 0.5;
  feed[1] = 0.5;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction2 = 1;
  real64 const vaporFraction2 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction2, expectedVaporFraction2, relTol );

  ////////////////////////////////////////

  kValues[0] = 54.866;
  kValues[1] = 54.866;

  feed[0] = 0.1;
  feed[1] = 0.9;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction3 = 1;
  real64 const vaporFraction3 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction3, expectedVaporFraction3, relTol );

  ////////////////////////////////////////

  kValues[0] = 1.09733;
  kValues[1] = 1.09733;

  feed[0] = 0.1;
  feed[1] = 0.9;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction4 = 1;
  real64 const vaporFraction4 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction4, expectedVaporFraction4, relTol );
}


TEST( RachfordRiceTest, testRachfordRiceFourComponents )
{
  constexpr integer numComps = 4;

  array1d< real64 > kValues( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  ////////////////////////////////////////

  kValues[0] = 537.526;
  kValues[1] = 0.00297518;
  kValues[2] = 6.19123e-08;
  kValues[3] = 1.44696;

  feed[0] = 0.099;
  feed[1] = 0.3;
  feed[2] = 0.6;
  feed[3] = 0.001;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction1 = 0.0975568;
  real64 const vaporFraction1 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction1, expectedVaporFraction1, relTol );

  ////////////////////////////////////////

  kValues[0] = 17.4329;
  kValues[1] = 0.000770753;
  kValues[2] = 1.18694e-06;
  kValues[3] = 0.00968376;

  feed[0] = 0.1;
  feed[1] = 0.1;
  feed[2] = 0.1;
  feed[3] = 0.7;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction2 = 0.045999;
  real64 const vaporFraction2 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction2, expectedVaporFraction2, relTol );

  ////////////////////////////////////////

  kValues[0] = 11.5506;
  kValues[1] = 0.000210682;
  kValues[2] = 1.57686e-08;
  kValues[3] = 0.0442361;

  feed[0] = 0.0984186;
  feed[1] = 0.297297;
  feed[2] = 0.593142;
  feed[3] = 0.0111427;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction3 = 0.0130259;
  real64 const vaporFraction3 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction3, expectedVaporFraction3, relTol );

}

// -----------------------------------------------------------------
// Derivative tests
// -----------------------------------------------------------------
template< int NC >
using TestData = std::tuple< real64 const, real64 const, Feed< NC > const >;

template< int NC >
struct TestFeed
{
  static std::array< real64 const, 5*NC > const crookston;
  static std::array< Feed< NC >, 4 > const feeds;
};

// Crookston correlation parameters obtained from fitting the PR EOS
template<>
std::array< real64 const, 10 > const TestFeed< 2 >::crookston = {
  -1.06433325e+01, 5.37086674e+07, 2.27201259e-06, 2.74613443e+00, 5.33213460e+02,
  4.41588002e-02, 6.20560890e+03, 3.68740831e-08, -5.54046510e+00, 5.98850203e+02
};

template<>
std::array< real64 const, 20 > const TestFeed< 4 >::crookston = {
  -1.06433325e+01, 5.37086674e+07, 2.27201259e-06, 2.74613443e+00, 5.33213460e+02,
  9.62608740e-01, 1.73365888e+05, -2.23745761e-07, -2.54802769e+01, 5.05617852e+02,
  4.41588002e-02, 6.20560890e+03, 3.68740831e-08, -5.54046510e+00, 5.98850203e+02,
  4.43507838e+00, 2.95487478e+06, -1.11104605e-06, 2.88802590e+01, 2.73149987e+02
};

template<>
std::array< Feed< 2 >, 4 > const TestFeed< 2 >::feeds = {
  Feed< 2 >{0.100000, 0.900000},
  Feed< 2 >{0.500000, 0.500000},
  Feed< 2 >{0.001000, 0.999000},
  Feed< 2 >{0.000010, 9.999990}
};

template<>
std::array< Feed< 4 >, 4 > const TestFeed< 4 >::feeds = {
  Feed< 4 >{0.030933, 0.319683, 0.637861, 0.011523},
  Feed< 4 >{0.999780, 0.000006, 0.000006, 0.000210},
  Feed< 4 >{0.000100, 0.348586, 0.637891, 0.012423},
  Feed< 4 >{0.100000, 0.347486, 0.550314, 0.000200}
};

template< int NC >
std::vector< TestData< NC > > generateTestData()
{
  std::vector< TestData< NC > > testData;
  std::array< real64 const, 3 > pressures( {1.83959e+06, 1.83959e+07, 1.83959e+08} );
  std::array< real64 const, 2 > temperatures( {2.97150e+02, 3.63000e+02} );
  for( const real64 pressure : pressures )
  {
    for( const real64 temperature : temperatures )
    {
      for( const auto & composition : TestFeed< NC >::feeds )
      {
        testData.emplace_back( pressure, temperature, composition );
      }
    }
  }
  return testData;
}

template< int NC >
class DerivativeTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDof = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  DerivativeTestFixture()
    : crookston( numComps, 5 )
  {
    auto const & data = TestFeed< NC >::crookston;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      for( integer kc = 0; kc < 5; ++kc )
      {
        crookston( ic, kc ) = data[5*ic + kc];
      }
    }
  }
  ~DerivativeTestFixture() = default;

  void evaluateCrookston( real64 const pressure, real64 const temperature,
                          arraySlice1d< real64 > const kValues ) const
  {
    real64 const p = pressure;
    real64 const t = temperature;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      auto [A, B, C, D, E] = get( ic );

      real64 const presTerm = A + B / p + C * p;
      real64 const tempTerm = LvArray::math::exp( -D / (t - E));

      kValues[ic] = presTerm * tempTerm;
    }
  }

  void evaluateCrookston( real64 const pressure, real64 const temperature,
                          arraySlice1d< real64 > const kValues, arraySlice2d< real64 > const kValueDerivs ) const
  {
    LvArray::forValuesInSlice( kValueDerivs, []( real64 & val ){ val = 0.0; } );
    real64 const p = pressure;
    real64 const t = temperature;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      auto [A, B, C, D, E] = get( ic );

      real64 const presTerm = A + B / p + C * p;
      real64 const tempTerm = LvArray::math::exp( -D / (t - E));

      kValues[ic] = presTerm * tempTerm;

      kValueDerivs( ic, Deriv::dP ) = (-B/(p*p) + C )*tempTerm;
      kValueDerivs( ic, Deriv::dT ) = presTerm * tempTerm * D /((t-E)*(t-E));
    }
  }

protected:
  std::tuple< real64 const, real64 const, real64 const, real64 const, real64 const > get( integer const ic ) const
  {
    return std::make_tuple( crookston( ic, 0 ), crookston( ic, 1 ), crookston( ic, 2 ), crookston( ic, 3 ), crookston( ic, 4 ));
  }
protected:
  stackArray2d< real64, 5*NC > crookston;
};

template< int NC >
class VapourFractionDerivativeTestFixture : public DerivativeTestFixture< NC >
{
public:
  using DerivativeTestFixture< NC >::numComps;
  using DerivativeTestFixture< NC >::numDof;
  using Deriv = typename DerivativeTestFixture< NC >::Deriv;
  using ParamType = typename DerivativeTestFixture< NC >::ParamType;
  using EOS = compositional::CubicEOSPhaseModel< compositional::PengRobinsonEOS >;

public:
  void testNumericalDerivatives( ParamType const & testData ) const
  {
    array1d< real64 > composition;
    real64 const pressure = std::get< 0 >( testData );
    real64 const temperature = std::get< 1 >( testData );
    TestFluid< NC >::createArray( composition, std::get< 2 >( testData ));

    real64 vapourFraction = -1.0;
    stackArray1d< integer, numComps > presentComponents;
    stackArray1d< real64, numComps > kValues( numComps );
    stackArray2d< real64, numComps *numDof > kValueDerivs( numComps, numDof );
    stackArray1d< real64, numDof > vapourFractionDerivs( numDof );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      presentComponents.emplace_back( ic );
    }

    this->evaluateCrookston( pressure, temperature, kValues, kValueDerivs );

    // Calculate vapour fraction
    vapourFraction = RachfordRice::solve( kValues, composition, presentComponents );

    // Calculate vapour fraction derivatives
    RachfordRice::computeDerivatives( kValues, kValueDerivs, composition, presentComponents, vapourFraction, vapourFractionDerivs );

    // Calculation function
    auto const calculateVapourFraction = [&]( real64 const p, real64 const t, auto const & zmf ) -> real64 {
      stackArray1d< real64, numComps > displacedKValues( numComps );
      this->evaluateCrookston( p, t, displacedKValues );
      return RachfordRice::solve( displacedKValues, zmf, presentComponents );
    };

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative(
      pressure, dp, vapourFractionDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return calculateVapourFraction( p, temperature, composition );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    geos::testing::internal::testNumericalDerivative(
      temperature, dT, vapourFractionDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return calculateVapourFraction( pressure, t, composition );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      geos::testing::internal::testNumericalDerivative(
        0.0, dz, vapourFractionDerivs[Deriv::dC+ic],
        [&]( real64 const z ) -> real64 {
        composition[ic] += z;
        real64 const fraction = calculateVapourFraction( pressure, temperature, composition );
        composition[ic] -= z;
        return fraction;
      } );
    }
  }
};

using VapourFractionDerivative2 = VapourFractionDerivativeTestFixture< 2 >;
using VapourFractionDerivative4 = VapourFractionDerivativeTestFixture< 4 >;

TEST_P( VapourFractionDerivative2, testNumericalDerivatives )
{
  testNumericalDerivatives( GetParam() );
}
TEST_P( VapourFractionDerivative4, testNumericalDerivatives )
{
  testNumericalDerivatives( GetParam() );
}

INSTANTIATE_TEST_SUITE_P(
  VapourFractionTest,
  VapourFractionDerivative2,
  ::testing::ValuesIn( generateTestData< 2 >() )
  );
INSTANTIATE_TEST_SUITE_P(
  VapourFractionTest,
  VapourFractionDerivative4,
  ::testing::ValuesIn( generateTestData< 4 >() )
  );
