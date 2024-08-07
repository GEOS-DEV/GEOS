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
#include "constitutive/fluid/multifluid/compositional/functions/SoreideWhitsonEOSModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

namespace geos
{

using namespace constitutive;
using namespace constitutive::compositional;

namespace testing
{

template< int NC >
using TestData = std::tuple<
  real64 const,         // Pressure
  real64 const,         // Temperature
  Feed< NC > const      // Input composition
  >;

template< integer NC, SoreideWhitsonPhaseType PHASE_TYPE, EquationOfStateType EOS_TYPE >
class SoreideWhitsonEOSModelTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDof = NC + 2;
  static constexpr real64 absTol = 1.0e-4;
  static constexpr real64 relTol = 1.0e-5;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using ParamType = TestData< NC >;
public:
  SoreideWhitsonEOSModelTestFixture();
  ~SoreideWhitsonEOSModelTestFixture() = default;

  void testPureCoefficients( ParamType const & testData )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();
    real64 const pressure = std::get< 0 >( testData );
    real64 const temperature = std::get< 1 >( testData );

    real64 aCoefficient = 0.0;
    real64 bCoefficient = 0.0;
    real64 daCoefficient_dp = 0.0;
    real64 dbCoefficient_dp = 0.0;
    real64 daCoefficient_dt = 0.0;
    real64 dbCoefficient_dt = 0.0;

    for( integer ic = 0; ic < numComps; ++ic )
    {
      SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::computePureCoefficients(
        ic,
        pressure,
        temperature,
        componentProperties,
        aCoefficient,
        bCoefficient,
        daCoefficient_dp,
        dbCoefficient_dp,
        daCoefficient_dt,
        dbCoefficient_dt );

      real64 const dp = 1.0e-4 * pressure;
      internal::testNumericalDerivative( pressure, dp, daCoefficient_dp,
                                         [&]( real64 p ) -> real64 {
        real64 a = 0.0, b = 0.0;
        SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::computePureCoefficients( ic, p, temperature, componentProperties, a, b );
        return a;
      }, absTol, relTol );
      internal::testNumericalDerivative( pressure, dp, dbCoefficient_dp,
                                         [&]( real64 p ) -> real64 {
        real64 a = 0.0, b = 0.0;
        SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::computePureCoefficients( ic, p, temperature, componentProperties, a, b );
        return b;
      }, absTol, relTol );

      real64 const dT = 1.0e-6 * temperature;
      internal::testNumericalDerivative( temperature, dT, daCoefficient_dt,
                                         [&]( real64 t ) -> real64 {
        real64 a = 0.0, b = 0.0;
        SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::computePureCoefficients( ic, pressure, t, componentProperties, a, b );
        return a;
      }, absTol, relTol );
      internal::testNumericalDerivative( temperature, dT, dbCoefficient_dt,
                                         [&]( real64 t ) -> real64 {
        real64 a = 0.0, b = 0.0;
        SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::computePureCoefficients( ic, pressure, t, componentProperties, a, b );
        return b;
      }, absTol, relTol );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

template< int NC >
struct FluidData {};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > create()
  {
    return TestFluid< 4 >::create( {Fluid::N2, Fluid::C1, Fluid::CO2, Fluid::H2O} );
  }
};

template< integer NC, SoreideWhitsonPhaseType PHASE_TYPE, EquationOfStateType EOS_TYPE >
SoreideWhitsonEOSModelTestFixture< NC, PHASE_TYPE, EOS_TYPE >::SoreideWhitsonEOSModelTestFixture():
  m_fluid( FluidData< NC >::create() )
{}

using PengRobinsonAqueous4 = SoreideWhitsonEOSModelTestFixture< 4, SoreideWhitsonPhaseType::Aqueous, EquationOfStateType::PengRobinson >;
using PengRobinsonVapour4 = SoreideWhitsonEOSModelTestFixture< 4, SoreideWhitsonPhaseType::Vapour, EquationOfStateType::PengRobinson >;

TEST_P( PengRobinsonAqueous4, testPureCoefficients )
{
  testPureCoefficients( GetParam() );
}
TEST_P( PengRobinsonVapour4, testPureCoefficients )
{
  testPureCoefficients( GetParam() );
}

template< int NC >
struct TestFeed {};

template<>
struct TestFeed< 2 >
{
  static std::array< Feed< 2 >, 3 > constexpr feeds = {
    Feed< 2 >{0.100000, 0.900000},
    Feed< 2 >{0.500000, 0.500000},
    Feed< 2 >{0.001000, 0.999000}
  };
};

template<>
struct TestFeed< 4 >
{
  static std::array< Feed< 4 >, 3 > constexpr feeds = {
    Feed< 4 >{0.030933, 0.319683, 0.637861, 0.011523},
    Feed< 4 >{0.000000, 0.349686, 0.637891, 0.012423},
    Feed< 4 >{0.000000, 0.349686, 0.650314, 0.000000}
  };
};

template< int NC >
std::vector< TestData< NC > > generateTestData()
{
  std::array< real64 const, 2 > pressures( {1.83959e+06, 1.83959e+08} );
  std::array< real64 const, 2 > temperatures( {2.97150e+02, 3.63000e+02} );
  std::vector< TestData< NC > > testData;
  for( const auto & composition : TestFeed< NC >::feeds )
  {
    for( const real64 pressure : pressures )
    {
      for( const real64 temperature : temperatures )
      {
        testData.emplace_back( pressure, temperature, composition );
      }
    }
  }
  return testData;
}

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSModelTest, PengRobinsonAqueous4, ::testing::ValuesIn( generateTestData< 4 >()) );
INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSModelTest, PengRobinsonVapour4, ::testing::ValuesIn( generateTestData< 4 >()) );

} // namespace testing

} // namespace geos
