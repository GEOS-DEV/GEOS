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
#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
using TestData = std::tuple<
  real64 const,     // 0 - pressure
  real64 const,     // 1 - temperature
  Feed< NC > const, // 2 - composition
  real64 const,     // 3 - expected molar density
  real64 const      // 4 - expected mass density
  >;

template< integer NC >
std::unique_ptr< TestFluid< NC > > createFluid();

template<>
std::unique_ptr< TestFluid< 4 > > createFluid< 4 >()
{
  return TestFluid< 4 >::create( {Fluid::N2, Fluid::C8H18, Fluid::C10H22, Fluid::H2O} );
}

template< typename EOS_TYPE, integer NC >
class TestDataTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDof = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  TestDataTestFixture()
    : m_fluid( createFluid< NC >() )
  {}
  ~TestDataTestFixture() = default;

  // Compares the calculated molar density against the expected value from PVT package
  void testMolarDensity( TestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    array1d< real64 > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 const expectedMolarDensity = std::get< 3 >( data );

    real64 const molarDensity = computeMolarDensity( pressure,
                                                     temperature,
                                                     composition.toSliceConst() );
    checkRelativeError( molarDensity, 0*expectedMolarDensity+molarDensity, internal::relTol, internal::absTol );
  }

  // Compares the calculated molar density derivatives against numerical calculated
  // finite difference values
  void testMolarDensityDerivative( TestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    array1d< real64 > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 constexpr molarDensityScale = 1.0e-3;

    real64 molarDensity = 0.0;
    stackArray1d< real64, numDof > molarDensityDerivs( numDof );

    computeMolarDensity( pressure,
                         temperature,
                         composition.toSliceConst(),
                         molarDensity,
                         molarDensityDerivs.toSlice() );

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative( pressure, dp, molarDensityDerivs[Deriv::dP], [&]( real64 const p ) -> real64 {
      return computeMolarDensity( p, temperature, composition );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative( temperature, dT, molarDensityDerivs[Deriv::dT], [&]( real64 const t ) -> real64 {
      return computeMolarDensity( pressure, t, composition );
    } );

    // -- Composition derivatives derivative
    real64 constexpr dz = 1.0e-6;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative( 0.0, dz, molarDensityScale*molarDensityDerivs[Deriv::dC + ic], [&]( real64 const z ) -> real64 {
        real64 const z_orig = composition[ic];
        composition[ic] += z;
        real64 const density = computeMolarDensity( pressure, temperature, composition );
        composition[ic] = z_orig;
        return molarDensityScale*density;
      } );
    }
  }

  // Compares the calculated mass density against the expected value from PVT package
  void testMassDensity( TestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    array1d< real64 > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 const expectedMassDensity = std::get< 4 >( data );

    real64 const massDensity = computeMassDensity( pressure,
                                                   temperature,
                                                   composition );
    checkRelativeError( massDensity, expectedMassDensity, internal::relTol, internal::absTol );
  }

  // Compares the calculated mass density derivatives against numerical calculated
  // finite difference values
  void testMassDensityDerivative( TestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    array1d< real64 > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 constexpr massDensityScale = 1.0e-2;

    real64 massDensity = 0.0;
    array1d< real64 > massDensityDerivs( numDof );

    computeMassDensity( pressure,
                        temperature,
                        composition,
                        massDensity,
                        massDensityDerivs );

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative( pressure, dp, massDensityDerivs[Deriv::dP], [&]( real64 const p ) -> real64 {
      return computeMassDensity( p, temperature, composition );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative( temperature, dT, massDensityDerivs[Deriv::dT], [&]( real64 const t ) -> real64 {
      return computeMassDensity( pressure, t, composition );
    } );

    // -- Composition derivatives
    real64 constexpr dz = 1.0e-6;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative( 0.0, dz, massDensityScale*massDensityDerivs[Deriv::dC + ic], [&]( real64 const z ) -> real64 {
        real64 const z0 = composition[ic];
        composition[ic] += z;
        real64 const density = computeMassDensity( pressure, temperature, composition );
        composition[ic] = z0;
        return massDensityScale*density;
      } );
    }
  }

private:
  real64 computeMolarDensity( real64 const pressure, real64 const temperature,
                              arraySlice1d< real64 const > const & composition ) const
  {
    real64 molarDensity = 0.0;
    stackArray1d< real64, numDof > molarDensityDerivs( numDof );
    computeMolarDensity( pressure, temperature, composition, molarDensity, molarDensityDerivs.toSlice() );
    return molarDensity;
  }

  void computeMolarDensity( real64 const pressure, real64 const temperature,
                            arraySlice1d< real64 const > const & composition,
                            real64 & molarDensity,
                            arraySlice1d< real64 > const molarDensityDerivs ) const
  {
    auto const componentProperties = this->m_fluid->createKernelWrapper();
    auto const volumeShift = componentProperties.m_componentVolumeShift;

    real64 compressibilityFactor = 0.0;
    stackArray1d< real64, numDof > compressibilityFactorDerivs( numDof );

    CubicEOSPhaseModel< EOS_TYPE >::
    computeCompressibilityFactorAndDerivs( numComps,
                                           pressure,
                                           temperature,
                                           composition,
                                           componentProperties,
                                           compressibilityFactor,
                                           compressibilityFactorDerivs.toSlice() );

    CompositionalProperties::computeMolarDensity( numComps,
                                                  pressure,
                                                  temperature,
                                                  composition,
                                                  volumeShift,
                                                  compressibilityFactor,
                                                  compressibilityFactorDerivs.toSliceConst(),
                                                  molarDensity,
                                                  molarDensityDerivs );
  }

  real64 computeMassDensity( real64 const pressure, real64 const temperature,
                             arraySlice1d< real64 const > const & composition ) const
  {
    real64 massDensity = 0.0;
    stackArray1d< real64, numDof > massDensityDerivs( numDof );
    computeMassDensity( pressure, temperature, composition, massDensity, massDensityDerivs.toSlice() );
    return massDensity;
  }

  void computeMassDensity( real64 const pressure, real64 const temperature,
                           arraySlice1d< real64 const > const & composition,
                           real64 & massDensity,
                           arraySlice1d< real64 > const massDensityDerivs ) const
  {
    auto const componentProperties = this->m_fluid->createKernelWrapper();
    auto const molecularWeight = componentProperties.m_componentMolarWeight;

    real64 molarDensity = 0.0;
    stackArray1d< real64, numDof > molarDensityDerivs( numDof );

    computeMolarDensity( pressure, temperature,
                         composition,
                         molarDensity,
                         molarDensityDerivs.toSlice() );

    CompositionalProperties::computeMassDensity( numComps,
                                                 composition,
                                                 molecularWeight,
                                                 molarDensity,
                                                 molarDensityDerivs.toSliceConst(),
                                                 massDensity,
                                                 massDensityDerivs );
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid;
};

using PR4Comp = TestDataTestFixture< PengRobinsonEOS, 4 >;

TEST_P( PR4Comp, testMolarDensity )
{
  testMolarDensity( GetParam() );
}

TEST_P( PR4Comp, testMolarDensityDerivative )
{
  testMolarDensityDerivative( GetParam() );
}

TEST_P( PR4Comp, testMassDensity )
{
  testMassDensity( GetParam() );
}

TEST_P( PR4Comp, testMassDensityDerivative )
{
  testMassDensityDerivative( GetParam() );
}

// Test data generated from PVT package
// All compositions are single phase

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(CompositionalPropertiesTest, PR4Comp,
  ::testing::ValuesIn<TestData< 4 >>({
    {1.00000e+05, 277.15, {0.000000, 0.495099, 0.495118, 0.009783}, 5.362153e+03, 6.819441e+02},
    {1.00000e+05, 277.15, {0.000652, 0.128231, 0.128281, 0.742836}, 1.504981e+04, 6.968129e+02},
    {1.00000e+05, 277.15, {0.855328, 0.000205, 0.000000, 0.144467}, 4.348763e+01, 1.156192e+00},
    {1.00000e+05, 288.65, {0.000507, 0.112984, 0.113029, 0.773480}, 1.617810e+04, 6.946347e+02},
    {1.00000e+05, 288.65, {0.777870, 0.000520, 0.000000, 0.221610}, 4.177800e+01, 1.079650e+00},
    {1.00000e+05, 288.65, {0.985235, 0.000000, 0.000000, 0.014765}, 4.169505e+01, 1.161865e+00},
    {1.00000e+05, 298.15, {0.653033, 0.000901, 0.000000, 0.346066}, 4.049439e+01, 9.974206e-01},
    {1.00000e+05, 298.15, {0.000506, 0.143046, 0.143326, 0.713122}, 1.366205e+04, 6.775557e+02},
    {1.00000e+05, 298.15, {0.582848, 0.000748, 0.000000, 0.416404}, 4.053152e+01, 9.692966e-01},
    {1.00000e+05, 298.15, {0.972848, 0.000000, 0.000000, 0.027152}, 4.036551e+01, 1.119817e+00},
    {1.00000e+05, 333.15, {0.000000, 0.477146, 0.477164, 0.045690}, 5.283543e+03, 6.510324e+02},
    {1.00000e+05, 333.15, {0.210877, 0.008984, 0.000001, 0.780137}, 3.640648e+01, 7.641050e-01},
    {1.00000e+05, 333.15, {0.818043, 0.000000, 0.000000, 0.181957}, 3.614886e+01, 9.468895e-01},
    {1.00000e+05, 372.15, {0.000000, 0.104818, 0.104822, 0.790360}, 3.285426e+01, 1.351168e+00},
    {1.00000e+05, 372.15, {0.000117, 0.356347, 0.549688, 0.093848}, 5.169992e+03, 6.235520e+02},
    {1.01325e+05, 277.15, {0.000000, 0.495099, 0.495118, 0.009783}, 5.362161e+03, 6.819450e+02},
    {1.01325e+05, 288.65, {0.000516, 0.112950, 0.112994, 0.773540}, 1.618129e+04, 6.946496e+02},
    {1.01325e+05, 288.65, {0.780815, 0.000514, 0.000000, 0.218672}, 4.233161e+01, 1.095180e+00},
    {1.01325e+05, 288.65, {0.000599, 0.132633, 0.132752, 0.734016}, 1.453988e+04, 6.874317e+02},
    {1.01325e+05, 288.65, {0.743752, 0.000426, 0.000000, 0.255823}, 4.235003e+01, 1.079604e+00},
    {1.01325e+05, 298.15, {0.000448, 0.114592, 0.114680, 0.770280}, 1.586896e+04, 6.870624e+02},
    {1.01325e+05, 298.15, {0.657746, 0.000890, 0.000000, 0.341364}, 4.103058e+01, 1.012517e+00},
    {1.01325e+05, 298.15, {0.000516, 0.142463, 0.142736, 0.714285}, 1.370089e+04, 6.777102e+02},
    {1.01325e+05, 333.15, {0.000000, 0.477146, 0.477164, 0.045690}, 5.283557e+03, 6.510339e+02},
    {1.01325e+05, 333.15, {0.000000, 0.000000, 0.000000, 1.000000}, 4.595576e+04, 8.279069e+02},
    {1.01325e+05, 333.15, {0.000000, 0.000000, 0.000000, 1.000000}, 4.595576e+04, 8.279069e+02},
    {1.01325e+05, 333.15, {0.820407, 0.000000, 0.000000, 0.179593}, 3.662781e+01, 9.603010e-01},
    {1.01325e+05, 372.15, {0.080408, 0.000000, 0.000000, 0.919592}, 3.300001e+01, 6.210347e-01},
    {5.00000e+06, 277.15, {0.000000, 0.495216, 0.495235, 0.009549}, 5.386778e+03, 6.852148e+02},
    {5.00000e+06, 277.15, {0.000000, 0.000000, 0.000000, 1.000000}, 4.775353e+04, 8.602942e+02},
    {5.00000e+06, 277.15, {0.029212, 0.107914, 0.107918, 0.754956}, 1.676639e+04, 7.058780e+02},
    {5.00000e+06, 288.65, {0.000000, 0.493093, 0.493111, 0.013796}, 5.365411e+03, 6.799847e+02},
    {5.00000e+06, 288.65, {0.000000, 0.000000, 0.000000, 1.000000}, 4.742029e+04, 8.542907e+02},
    {5.00000e+06, 288.65, {0.999549, 0.000000, 0.000000, 0.000451}, 2.123111e+03, 5.946598e+01},
    {5.00000e+06, 298.15, {0.000000, 0.490855, 0.490873, 0.018272}, 5.350468e+03, 6.754508e+02},
    {5.00000e+06, 298.15, {0.000000, 0.000000, 0.000000, 1.000000}, 4.713460e+04, 8.491439e+02},
    {5.00000e+06, 298.15, {0.030096, 0.107834, 0.107839, 0.754231}, 1.641282e+04, 6.908503e+02},
    {5.00000e+06, 333.15, {0.031858, 0.107709, 0.107720, 0.752713}, 1.573686e+04, 6.622529e+02},
    {5.00000e+06, 333.15, {0.964491, 0.000264, 0.000000, 0.035245}, 1.825541e+03, 5.053788e+01},
    {5.00000e+06, 333.15, {0.035989, 0.120533, 0.120574, 0.722904}, 1.460777e+04, 6.566977e+02},
    {5.00000e+06, 333.15, {0.961874, 0.000247, 0.000000, 0.037879}, 1.826826e+03, 5.052265e+01},
    {5.00000e+06, 372.15, {0.037408, 0.121753, 0.121923, 0.718916}, 1.364945e+04, 6.177002e+02},
    {5.00000e+06, 372.15, {0.889694, 0.001014, 0.000003, 0.109288}, 1.644155e+03, 4.440610e+01},
    {1.00000e+07, 277.15, {0.000004, 0.000000, 0.000000, 0.999996}, 4.777672e+04, 8.607139e+02},
    {1.00000e+07, 277.15, {0.999835, 0.000000, 0.000000, 0.000165}, 4.471992e+03, 1.252683e+02},
    {1.00000e+07, 288.65, {0.000000, 0.493260, 0.493279, 0.013460}, 5.390317e+03, 6.833402e+02},
    {1.00000e+07, 288.65, {0.000000, 0.000000, 0.000000, 1.000000}, 4.744610e+04, 8.547556e+02},
    {1.00000e+07, 288.65, {0.055934, 0.104932, 0.104936, 0.734198}, 1.678938e+04, 7.002932e+02},
    {1.00000e+07, 298.15, {0.000000, 0.491084, 0.491102, 0.017814}, 5.377368e+03, 6.791183e+02},
    {1.00000e+07, 298.15, {0.000000, 0.000000, 0.000000, 1.000000}, 4.716269e+04, 8.496500e+02},
    {1.00000e+07, 298.15, {0.056951, 0.104818, 0.104822, 0.733409}, 1.662695e+04, 6.932690e+02},
    {1.00000e+07, 298.15, {0.065451, 0.116299, 0.116311, 0.701940}, 1.552571e+04, 6.879882e+02},
    {1.00000e+07, 298.15, {0.991760, 0.000074, 0.000000, 0.008166}, 4.105988e+03, 1.147137e+02},
    {1.00000e+07, 333.15, {0.000000, 0.478424, 0.478442, 0.043135}, 5.359052e+03, 6.618465e+02},
    {1.00000e+07, 333.15, {0.000000, 0.000000, 0.000000, 1.000000}, 4.603212e+04, 8.292824e+02},
    {1.00000e+07, 333.15, {0.996818, 0.000000, 0.000000, 0.003182}, 3.591935e+03, 1.005080e+02},
    {1.00000e+07, 372.15, {0.000000, 0.453197, 0.453215, 0.093588}, 5.412579e+03, 6.383520e+02},
    {1.00000e+08, 277.15, {0.000000, 0.496608, 0.496627, 0.006766}, 5.621343e+03, 7.167775e+02},
    {1.00000e+08, 288.65, {0.156951, 0.104818, 0.104822, 0.633409}, 1.719322e+04, 7.340702e+02},
    {1.00000e+08, 288.65, {0.000031, 0.000000, 0.000000, 0.999969}, 4.786341e+04, 8.622886e+02},
    {1.00000e+08, 288.65, {0.999447, 0.000000, 0.000000, 0.000553}, 2.330021e+04, 6.525893e+02},
    {1.00000e+08, 298.15, {0.000000, 0.493623, 0.493642, 0.012735}, 5.616788e+03, 7.125000e+02},
    {1.00000e+08, 333.15, {0.010000, 0.000000, 0.000000, 0.990000}, 4.641206e+04, 8.407676e+02},
    {1.00000e+08, 372.15, {0.156951, 0.104818, 0.104822, 0.633409}, 1.615558e+04, 6.897677e+02},
    {1.00000e+08, 372.15, {0.000671, 0.000000, 0.000000, 0.999329}, 4.540878e+04, 8.183575e+02},
    {1.00000e+08, 372.15, {0.991955, 0.000000, 0.000000, 0.008045}, 1.979968e+04, 5.530638e+02}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
