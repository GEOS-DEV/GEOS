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
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

static constexpr real64 relTol = 1.0e-5;
static constexpr real64 absTol = 1.0e-8;

TEST( CubicEOSTest, testCubicEOSTwoComponentsSRK )
{
  constexpr integer numComps = 2;

  auto fluid = TestFluid< numComps >::create( {"C1", "C5"} );

  auto criticalPressure = fluid->getCriticalPressure();
  auto criticalTemperature = fluid->getCriticalTemperature();
  auto omega = fluid->getAcentricFactor();
  real64 binaryInteractionCoefficients = 0.0; // not implemented yet

  real64 pressure = 0.0;
  real64 temperature = 0.0;
  array1d< real64 > composition( numComps );
  array1d< real64 > logFugacityCoefficients( numComps );
  array1d< real64 > expectedLogFugacityCoefficients( numComps );

  ///////////////////////////////////////////

  pressure = 1e6;
  temperature = 350;
  composition[0] = 0.1;
  composition[1] = 0.9;

  expectedLogFugacityCoefficients[0] = 0.0126163;
  expectedLogFugacityCoefficients[1] = -0.00820777;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

  ///////////////////////////////////////////

  pressure = 5e7;
  temperature = 350;
  composition[0] = 0.1;
  composition[1] = 0.9;

  expectedLogFugacityCoefficients[0] = 0.481514;
  expectedLogFugacityCoefficients[1] = -0.0701117;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

  ///////////////////////////////////////////

  pressure = 1e6;
  temperature = 350;
  composition[0] = 0.5;
  composition[1] = 0.5;

  expectedLogFugacityCoefficients[0] = 0.00721367;
  expectedLogFugacityCoefficients[1] = -0.00589892;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

  ///////////////////////////////////////////

  pressure = 5e7;
  temperature = 350;
  composition[0] = 0.5;
  composition[1] = 0.5;

  expectedLogFugacityCoefficients[0] = 0.334027;
  expectedLogFugacityCoefficients[1] = -0.00629384;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );

}

TEST( CubicEOSTest, testCubicEOSFourComponentsPR )
{
  constexpr integer numComps = 4;

  auto fluid = TestFluid< numComps >::create( {"N2", "C8", "C10", "H2O"} );

  auto criticalPressure = fluid->getCriticalPressure();
  auto criticalTemperature = fluid->getCriticalTemperature();
  auto omega = fluid->getAcentricFactor();
  real64 binaryInteractionCoefficients = 0.0; // not implemented yet

  real64 pressure = 0.0;
  real64 temperature = 0.0;
  array1d< real64 > composition( numComps );
  array1d< real64 > logFugacityCoefficients( numComps );
  array1d< real64 > expectedLogFugacityCoefficients( numComps );

  ///////////////////////////////////////////

  pressure = 1e7;
  temperature = 297.15;
  composition[0] = 0.0569514;
  composition[1] = 0.104818;
  composition[2] = 0.104822;
  composition[3] = 0.733409;

  expectedLogFugacityCoefficients[0] = 2.8298;
  expectedLogFugacityCoefficients[1] = -8.88628;
  expectedLogFugacityCoefficients[2] = -17.0201;
  expectedLogFugacityCoefficients[3] = -5.33003;

  CubicEOSPhaseModel< PengRobinsonEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 1e5;
  temperature = 297.15;
  composition[0] = 0.00185559;
  composition[1] = 0.332324;
  composition[2] = 0.664862;
  composition[3] = 0.000958244;

  expectedLogFugacityCoefficients[0] = 6.28652;
  expectedLogFugacityCoefficients[1] = -5.83771;
  expectedLogFugacityCoefficients[2] = -16.638;
  expectedLogFugacityCoefficients[3] = 0.361984;

  CubicEOSPhaseModel< PengRobinsonEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 4.78429e+06;
  temperature = 297.15;
  composition[0] = 0.0566196;
  composition[1] = 0.31411;
  composition[2] = 0.628223;
  composition[3] = 0.001047;

  expectedLogFugacityCoefficients[0] = 2.49484;
  expectedLogFugacityCoefficients[1] = -9.36508;
  expectedLogFugacityCoefficients[2] = -19.8123;
  expectedLogFugacityCoefficients[3] = -3.42481;

  CubicEOSPhaseModel< PengRobinsonEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

}

TEST( CubicEOSTest, testCubicEOSFourComponentsSRK )
{
  constexpr integer numComps = 4;

  auto fluid = TestFluid< numComps >::create( {"N2", "C8", "C10", "H2O"} );

  auto criticalPressure = fluid->getCriticalPressure();
  auto criticalTemperature = fluid->getCriticalTemperature();
  auto omega = fluid->getAcentricFactor();
  real64 binaryInteractionCoefficients = 0.0; // not implemented yet

  real64 pressure = 0.0;
  real64 temperature = 0.0;
  array1d< real64 > composition( 4 );
  array1d< real64 > logFugacityCoefficients( 4 );
  array1d< real64 > expectedLogFugacityCoefficients( 4 );

  ///////////////////////////////////////////

  pressure = 1e7;
  temperature = 297.15;
  composition[0] = 0.994214;
  composition[1] = 6.05198e-05;
  composition[2] = 5.98122e-08;
  composition[3] = 0.00572563;

  expectedLogFugacityCoefficients[0] = 0.00588361;
  expectedLogFugacityCoefficients[1] = -1.44445;
  expectedLogFugacityCoefficients[2] = -2.83086;
  expectedLogFugacityCoefficients[3] = -0.618972;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 1e5;
  temperature = 297.15;
  composition[0] = 0.997965;
  composition[1] = 0.000851981;
  composition[2] = 2.89283e-08;
  composition[3] = 0.00118249;

  expectedLogFugacityCoefficients[0] = -5.94544e-05;
  expectedLogFugacityCoefficients[1] = -0.0168209;
  expectedLogFugacityCoefficients[2] = -0.0334318;
  expectedLogFugacityCoefficients[3] = -0.00664411;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

  ///////////////////////////////////////////

  pressure = 1.83959e+06;
  temperature = 297.15;
  composition[0] = 0.0309329;
  composition[1] = 0.319683;
  composition[2] = 0.637861;
  composition[3] = 0.011523;

  expectedLogFugacityCoefficients[0] = 3.47428;
  expectedLogFugacityCoefficients[1] = -8.75355;
  expectedLogFugacityCoefficients[2] = -19.6075;
  expectedLogFugacityCoefficients[3] = -2.69792;

  CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
  compute( numComps,
           pressure, temperature, composition,
           criticalPressure, criticalTemperature, omega, binaryInteractionCoefficients,
           logFugacityCoefficients );

  checkRelativeError( logFugacityCoefficients[0], expectedLogFugacityCoefficients[0], relTol );
  checkRelativeError( logFugacityCoefficients[1], expectedLogFugacityCoefficients[1], relTol );
  checkRelativeError( logFugacityCoefficients[2], expectedLogFugacityCoefficients[2], relTol );
  checkRelativeError( logFugacityCoefficients[3], expectedLogFugacityCoefficients[3], relTol );

}

// -----------------------------------------------------------------
// Derivative tests
// -----------------------------------------------------------------

template< int NC >
using TestData = std::tuple< real64 const, real64 const, Feed< NC > const >;

template< typename EOS, int NC >
class DerivativeTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  using ParamType = std::tuple< real64 const, real64 const, Feed< NC > const >;
public:
  DerivativeTestFixture();
  ~DerivativeTestFixture() = default;

protected:
  void checkDerivative( real64 const a, real64 const b, string const & name ) const
  {
    checkRelativeError( a, b, relTol, absTol, name );
  }

  std::unique_ptr< TestFluid< NC > > m_fluid;
};

template<>
DerivativeTestFixture< PengRobinsonEOS, 2 >::DerivativeTestFixture()
  : m_fluid( TestFluid< 2 >::create( {"C1", "C5"} ))
{}
template<>
DerivativeTestFixture< SoaveRedlichKwongEOS, 2 >::DerivativeTestFixture()
  : m_fluid( TestFluid< 2 >::create( {"C1", "C5"} ))
{}

template<>
DerivativeTestFixture< PengRobinsonEOS, 4 >::DerivativeTestFixture()
  : m_fluid( TestFluid< 4 >::create( {"N2", "C8", "C10", "H2O"} ))
{}
template<>
DerivativeTestFixture< SoaveRedlichKwongEOS, 4 >::DerivativeTestFixture()
  : m_fluid( TestFluid< 4 >::create( {"N2", "C8", "C10", "H2O"} ))
{}

template< typename EOS, int NC >
class MixCoeffDerivativeTestFixture : public DerivativeTestFixture< EOS, NC >
{
public:
  using DerivativeTestFixture< EOS, NC >::numComps;
  using ParamType = typename DerivativeTestFixture< EOS, NC >::ParamType;
public:
  void testNumericalDerivatives( ParamType const & testData ) const
  {
    auto criticalPressure = this->m_fluid->getCriticalPressure();
    auto criticalTemperature = this->m_fluid->getCriticalTemperature();
    auto omega = this->m_fluid->getAcentricFactor();
    real64 binaryInteractionCoefficients = 0.0; // not implemented yet

    array1d< real64 > aPureCoefficient( numComps );
    array1d< real64 > bPureCoefficient( numComps );
    real64 aMixtureCoefficient = 0.0;
    real64 bMixtureCoefficient = 0.0;
    real64 currentAMixtureCoefficient = 0.0;
    real64 currentBMixtureCoefficient = 0.0;
    real64 fdDerivative = 0.0;

    real64 daMixtureCoefficient_dp = 0.0;
    real64 dbMixtureCoefficient_dp = 0.0;
    real64 daMixtureCoefficient_dt = 0.0;
    real64 dbMixtureCoefficient_dt = 0.0;
    array1d< real64 > daMixtureCoefficient_dz( numComps );
    array1d< real64 > dbMixtureCoefficient_dz( numComps );

    array1d< real64 > composition;
    real64 const pressure = std::get< 0 >( testData );
    real64 const temperature = std::get< 1 >( testData );
    TestFluid< NC >::createArray( composition, std::get< 2 >( testData ));

    auto computeCoefficients = [&]( real64 const p, real64 const t, auto const & z, real64 & a, real64 & b ){
      CubicEOSPhaseModel< EOS >::computeMixtureCoefficients(
        numComps,
        p, t, z,
        criticalPressure, criticalTemperature, omega,
        binaryInteractionCoefficients,
        aPureCoefficient,
        bPureCoefficient,
        a, b
        );
    };

    // Calculate values
    computeCoefficients( pressure, temperature, composition, aMixtureCoefficient, bMixtureCoefficient );
    // Calculate derivatives
    CubicEOSPhaseModel< EOS >::computeMixtureCoefficients(
      numComps,
      pressure,
      temperature,
      composition,
      criticalPressure,
      criticalTemperature,
      omega,
      binaryInteractionCoefficients,
      aPureCoefficient,
      bPureCoefficient,
      aMixtureCoefficient,
      bMixtureCoefficient,
      daMixtureCoefficient_dp,
      dbMixtureCoefficient_dp,
      daMixtureCoefficient_dt,
      dbMixtureCoefficient_dt,
      daMixtureCoefficient_dz,
      dbMixtureCoefficient_dz );
    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    computeCoefficients( pressure-dp, temperature, composition, currentAMixtureCoefficient, currentBMixtureCoefficient );
    fdDerivative = -(currentAMixtureCoefficient - aMixtureCoefficient) / dp;
    this->checkDerivative( daMixtureCoefficient_dp, fdDerivative, "Mixing Coeff A left pressure derivative" );
    fdDerivative = -(currentBMixtureCoefficient - bMixtureCoefficient) / dp;
    this->checkDerivative( dbMixtureCoefficient_dp, fdDerivative, "Mixing Coeff B left pressure derivative" );
    computeCoefficients( pressure+dp, temperature, composition, currentAMixtureCoefficient, currentBMixtureCoefficient );
    fdDerivative = (currentAMixtureCoefficient - aMixtureCoefficient) / dp;
    this->checkDerivative( daMixtureCoefficient_dp, fdDerivative, "Mixing Coeff A right pressure derivative" );
    fdDerivative = (currentBMixtureCoefficient - bMixtureCoefficient) / dp;
    this->checkDerivative( dbMixtureCoefficient_dp, fdDerivative, "Mixing Coeff B right pressure derivative" );
    // -- Temperature derivative
    real64 const dt = 1.0e-6 * temperature;
    computeCoefficients( pressure, temperature-dt, composition, currentAMixtureCoefficient, currentBMixtureCoefficient );
    fdDerivative = -(currentAMixtureCoefficient - aMixtureCoefficient) / dt;
    this->checkDerivative( daMixtureCoefficient_dt, fdDerivative, "Mixing Coeff A left temperature derivative" );
    fdDerivative = -(currentBMixtureCoefficient - bMixtureCoefficient) / dt;
    this->checkDerivative( dbMixtureCoefficient_dt, fdDerivative, "Mixing Coeff B left temperature derivative" );
    computeCoefficients( pressure, temperature+dt, composition, currentAMixtureCoefficient, currentBMixtureCoefficient );
    fdDerivative = (currentAMixtureCoefficient - aMixtureCoefficient) / dt;
    this->checkDerivative( daMixtureCoefficient_dt, fdDerivative, "Mixing Coeff A right temperature derivative" );
    fdDerivative = (currentBMixtureCoefficient - bMixtureCoefficient) / dt;
    this->checkDerivative( dbMixtureCoefficient_dt, fdDerivative, "Mixing Coeff B right temperature derivative" );
    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      composition[ic] -= dz;
      computeCoefficients( pressure, temperature, composition, currentAMixtureCoefficient, currentBMixtureCoefficient );
      fdDerivative = -(currentAMixtureCoefficient - aMixtureCoefficient) / dz;
      this->checkDerivative( daMixtureCoefficient_dz[ic], fdDerivative, "Mixing Coeff A left composition derivative" );
      fdDerivative = -(currentBMixtureCoefficient - bMixtureCoefficient) / dz;
      this->checkDerivative( dbMixtureCoefficient_dz[ic], fdDerivative, "Mixing Coeff B left composition derivative" );
      composition[ic] += 2.0*dz;
      computeCoefficients( pressure, temperature, composition, currentAMixtureCoefficient, currentBMixtureCoefficient );
      fdDerivative = (currentAMixtureCoefficient - aMixtureCoefficient) / dz;
      this->checkDerivative( daMixtureCoefficient_dz[ic], fdDerivative, "Mixing Coeff A right composition derivative" );
      fdDerivative = (currentBMixtureCoefficient - bMixtureCoefficient) / dz;
      this->checkDerivative( dbMixtureCoefficient_dz[ic], fdDerivative, "Mixing Coeff B right composition derivative" );
      composition[ic] -= dz;
    }
  }
};

using MixCoeffDerivativePR2TestFixture = MixCoeffDerivativeTestFixture< PengRobinsonEOS, 2 >;
using MixCoeffDerivativePR4TestFixture = MixCoeffDerivativeTestFixture< PengRobinsonEOS, 4 >;
using MixCoeffDerivativeSRK2TestFixture = MixCoeffDerivativeTestFixture< SoaveRedlichKwongEOS, 2 >;
using MixCoeffDerivativeSRK4TestFixture = MixCoeffDerivativeTestFixture< SoaveRedlichKwongEOS, 4 >;

TEST_P( MixCoeffDerivativePR2TestFixture, testNumericalDerivatives )
{
  testNumericalDerivatives( GetParam() );
}
TEST_P( MixCoeffDerivativePR4TestFixture, testNumericalDerivatives )
{
  testNumericalDerivatives( GetParam() );
}
TEST_P( MixCoeffDerivativeSRK2TestFixture, testNumericalDerivatives )
{
  testNumericalDerivatives( GetParam() );
}
TEST_P( MixCoeffDerivativeSRK4TestFixture, testNumericalDerivatives )
{
  testNumericalDerivatives( GetParam() );
}

// 2-component fluid test
INSTANTIATE_TEST_SUITE_P(
  CubicEOSTest,
  MixCoeffDerivativePR2TestFixture,
  ::testing::Values(
    TestData< 2 >{ 1.83959e+06, 2.97150e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+06, 2.97150e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+06, 2.97150e+02, {0.001000, 0.999000} },
    TestData< 2 >{ 1.83959e+08, 2.97150e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+08, 2.97150e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+08, 2.97150e+02, {0.001000, 0.999000} },
    TestData< 2 >{ 1.83959e+06, 3.63000e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+06, 3.63000e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+06, 3.63000e+02, {0.001000, 0.999000} },
    TestData< 2 >{ 1.83959e+08, 3.63000e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+08, 3.63000e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+08, 3.63000e+02, {0.001000, 0.999000} }
    )
  );
INSTANTIATE_TEST_SUITE_P(
  CubicEOSTest,
  MixCoeffDerivativeSRK2TestFixture,
  ::testing::Values(
    TestData< 2 >{ 1.83959e+06, 2.97150e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+06, 2.97150e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+06, 2.97150e+02, {0.001000, 0.999000} },
    TestData< 2 >{ 1.83959e+08, 2.97150e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+08, 2.97150e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+08, 2.97150e+02, {0.001000, 0.999000} },
    TestData< 2 >{ 1.83959e+06, 3.63000e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+06, 3.63000e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+06, 3.63000e+02, {0.001000, 0.999000} },
    TestData< 2 >{ 1.83959e+08, 3.63000e+02, {0.100000, 0.900000} },
    TestData< 2 >{ 1.83959e+08, 3.63000e+02, {0.500000, 0.500000} },
    TestData< 2 >{ 1.83959e+08, 3.63000e+02, {0.001000, 0.999000} }
    )
  );

// 4-component fluid test
INSTANTIATE_TEST_SUITE_P(
  CubicEOSTest,
  MixCoeffDerivativePR4TestFixture,
  ::testing::Values(
    TestData< 4 >{ 1.83959e+06, 2.97150e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+06, 2.97150e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+06, 2.97150e+02, {0.000000, 0.349686, 0.650314, 0.000000} },
    TestData< 4 >{ 1.83959e+08, 2.97150e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+08, 2.97150e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+08, 2.97150e+02, {0.000000, 0.349686, 0.650314, 0.000000} },
    TestData< 4 >{ 1.83959e+06, 3.63000e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+06, 3.63000e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+06, 3.63000e+02, {0.000000, 0.349686, 0.650314, 0.000000} },
    TestData< 4 >{ 1.83959e+08, 3.63000e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+08, 3.63000e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+08, 3.63000e+02, {0.000000, 0.349686, 0.650314, 0.000000} }
    )
  );

INSTANTIATE_TEST_SUITE_P(
  CubicEOSTest,
  MixCoeffDerivativeSRK4TestFixture,
  ::testing::Values(
    TestData< 4 >{ 1.83959e+06, 2.97150e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+06, 2.97150e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+06, 2.97150e+02, {0.000000, 0.349686, 0.650314, 0.000000} },
    TestData< 4 >{ 1.83959e+08, 2.97150e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+08, 2.97150e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+08, 2.97150e+02, {0.000000, 0.349686, 0.650314, 0.000000} },
    TestData< 4 >{ 1.83959e+06, 3.63000e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+06, 3.63000e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+06, 3.63000e+02, {0.000000, 0.349686, 0.650314, 0.000000} },
    TestData< 4 >{ 1.83959e+08, 3.63000e+02, {0.030933, 0.319683, 0.637861, 0.011523} },
    TestData< 4 >{ 1.83959e+08, 3.63000e+02, {0.000000, 0.349686, 0.637891, 0.012423} },
    TestData< 4 >{ 1.83959e+08, 3.63000e+02, {0.000000, 0.349686, 0.650314, 0.000000} }
    )
  );
