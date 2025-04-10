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
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"

#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"
#include "functions/FunctionManager.hpp"
#include "common/PhysicsConstants.hpp"

#include <conduit.hpp>

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using DensityData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // phase composition
  real64 const,       // expected molar density
  real64 const        // expected mass density
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    return TestFluid< 4 >::create( {Fluid::CO2, Fluid::H2O, Fluid::C1, Fluid::N2} );
  }
};

template< int NC, EquationOfStateType EOS_TYPE, int SALINITY = 0 >
class CompositionalPhillipsBrineDensityTestFixture :  public ::testing::TestWithParam< DensityData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  CompositionalPhillipsBrineDensityTestFixture()
    : m_parent( "parent", m_node ),
    m_fluid( FluidData< NC >::createFluid() )
  {
    m_functionManager = std::make_unique< FunctionManager >( "FunctionManager", &m_parent );
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();
    m_parameters = PhillipsBrineDensity::createParameters( std::make_unique< ModelParameters >() );

    auto equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EOS_TYPE );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto coordinates = const_cast< PressureTemperatureCoordinates * >(m_parameters->get< PressureTemperatureCoordinates >());
    createPressurePoints( coordinates->m_pressureCoordinates );
    createTemperaturePoints( coordinates->m_temperatureCoordinates );

    auto brineSalinity = const_cast< BrineSalinity * >(m_parameters->get< BrineSalinity >());
    real64 const massFraction = 1.0e-6*SALINITY;
    brineSalinity->m_salinity = massFraction / 58.44e-3;

    string const name = GEOS_FMT( "PhaseDensity{}{}", eosName, SALINITY );
    m_density = std::make_unique< PhillipsBrineDensity >( name, componentProperties, 0, *m_parameters );
  }

  ~CompositionalPhillipsBrineDensityTestFixture() = default;

  void testDensityValues( DensityData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));
    real64 const expectedMolarDensity = std::get< 3 >( data );
    real64 const expectedMassDensity = std::get< 4 >( data );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = m_density->createKernelWrapper();

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           phaseComposition.toSliceConst(),
                           molarDensity,
                           tempDerivs.toSlice(),
                           massDensity,
                           tempDerivs.toSlice(),
                           false );

    checkRelativeError( molarDensity, expectedMolarDensity, relTol, absTol );
    checkRelativeError( massDensity, expectedMassDensity, relTol, absTol );
  }

  void testDensityDerivatives( DensityData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = m_density->createKernelWrapper();

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    stackArray2d< real64, 3*numDofs > derivSpace( 3, numDofs );
    arraySlice1d< real64 > molarDensityDerivs = derivSpace[0];
    arraySlice1d< real64 > massDensityDerivs = derivSpace[1];
    arraySlice1d< real64 > tempDerivs = derivSpace[2];

    integer constexpr numValues = 2;
    stackArray1d< real64, numValues > derivatives( numValues );

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           phaseComposition.toSliceConst(),
                           molarDensity,
                           molarDensityDerivs,
                           massDensity,
                           massDensityDerivs,
                           false );

    auto const concactDerivatives = [&]( integer const idof ){
      derivatives[0] = molarDensityDerivs[idof];
      derivatives[1] = massDensityDerivs[idof];
    };

    // Compare against numerical derivatives
    // -- Pressure derivative
    concactDerivatives( Deriv::dP );
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative< numValues >(
      pressure, dp, derivatives.toSliceConst(),
      [&]( real64 const p, auto & values ){
      kernelWrapper.compute( componentProperties,
                             p,
                             temperature,
                             phaseComposition.toSliceConst(),
                             values[0],
                             tempDerivs,
                             values[1],
                             tempDerivs,
                             false );
    } );

    // -- Temperature derivative
    concactDerivatives( Deriv::dT );
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative< numValues >(
      temperature, dT, derivatives.toSliceConst(),
      [&]( real64 const t, auto & values ){
      kernelWrapper.compute( componentProperties,
                             pressure,
                             t,
                             phaseComposition.toSliceConst(),
                             values[0],
                             tempDerivs,
                             values[1],
                             tempDerivs,
                             false );
    } );

    // -- Composition derivatives
    real64 constexpr dz = 1.0e-6;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      concactDerivatives( Deriv::dC + ic );
      internal::testNumericalDerivative< numValues >(
        0.0, dz, derivatives.toSliceConst(),
        [&]( real64 const z, auto & values ){
        real64 const z_old = phaseComposition[ic];
        phaseComposition[ic] += z;
        kernelWrapper.compute( componentProperties,
                               pressure,
                               temperature,
                               phaseComposition.toSliceConst(),
                               values[0],
                               tempDerivs,
                               values[1],
                               tempDerivs,
                               false );
        phaseComposition[ic] = z_old;
      } );
    }
  }

private:
  static void createPressurePoints( array1d< real64 > & pressureCoordinates )
  {
    int constexpr n = 20;
    real64 constexpr minPressure = 0.995e5;
    real64 constexpr maxPressure = 1000.005e5;
    real64 const r = pow( maxPressure/minPressure, 1.0 / n );
    pressureCoordinates.resize( n+1 );
    pressureCoordinates[0] = minPressure;
    for( integer i = 1; i <= n; i++ )
    {
      pressureCoordinates[i] = pressureCoordinates[i-1]*r;
    }
  }

  static void createTemperaturePoints( array1d< real64 > & temperatureCoordinates )
  {
    int constexpr n = 20;
    real64 constexpr minTemperature = 12.0 + constants::zeroDegreesCelsiusInKelvin;
    real64 constexpr maxTemperature = 120.0 + constants::zeroDegreesCelsiusInKelvin;
    real64 constexpr dT = (maxTemperature - minTemperature) / n;
    temperatureCoordinates.resize( n+1 );
    for( integer i = 0; i <= n; i++ )
    {
      temperatureCoordinates[i] = minTemperature + i*dT + 0.01;
    }
  }

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  std::unique_ptr< FunctionManager > m_functionManager{};
  std::unique_ptr< PhillipsBrineDensity > m_density{};
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using PengRobinson = CompositionalPhillipsBrineDensityTestFixture< 4, EquationOfStateType::PengRobinson >;
using PengRobinson100K = CompositionalPhillipsBrineDensityTestFixture< 4, EquationOfStateType::PengRobinson, 100000 >;

TEST_P( PengRobinson, testDensityValues )
{
  testDensityValues( GetParam() );
}

TEST_P( PengRobinson100K, testDensityValues )
{
  testDensityValues( GetParam() );
}

TEST_P( PengRobinson, testDensityDerivatives )
{
  testDensityDerivatives( GetParam() );
}

TEST_P( PengRobinson100K, testDensityDerivatives )
{
  testDensityDerivatives( GetParam() );
}

/* UNCRUSTIFY-OFF */

// Test data

INSTANTIATE_TEST_SUITE_P(
  CompositionalPhillipsBrineDensityTest, PengRobinson,
  ::testing::ValuesIn( {
    // Pure water: Mass density compared against CO2-Brine implementation
    DensityData<4>{1.0e+05, 288.15, {0.000, 1.000, 0.000, 0.000}, 5.54544e+04, 9.99012e+02},
    DensityData<4>{1.0e+06, 288.15, {0.000, 1.000, 0.000, 0.000}, 5.54769e+04, 9.99416e+02},
    DensityData<4>{5.0e+07, 288.15, {0.000, 1.000, 0.000, 0.000}, 5.67137e+04, 1.02170e+03},
    DensityData<4>{1.0e+05, 293.15, {0.000, 1.000, 0.000, 0.000}, 5.54069e+04, 9.98155e+02},
    DensityData<4>{1.0e+06, 293.15, {0.000, 1.000, 0.000, 0.000}, 5.54293e+04, 9.98559e+02},
    DensityData<4>{5.0e+07, 293.15, {0.000, 1.000, 0.000, 0.000}, 5.66650e+04, 1.02082e+03},
    DensityData<4>{1.0e+05, 373.15, {0.000, 1.000, 0.000, 0.000}, 5.27943e+04, 9.51089e+02}, // Undersaturated value
    DensityData<4>{1.0e+06, 373.15, {0.000, 1.000, 0.000, 0.000}, 5.32169e+04, 9.58703e+02},
    DensityData<4>{5.0e+07, 373.15, {0.000, 1.000, 0.000, 0.000}, 5.44032e+04, 9.80073e+02},
    // Mixtures
    DensityData<4>{1.0e+05, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.49382e+04, 9.99502e+02},
    DensityData<4>{1.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.49604e+04, 9.99906e+02},
    DensityData<4>{5.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.50593e+04, 1.00170e+03},
    DensityData<4>{5.0e+07, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.61822e+04, 1.02213e+03},
    DensityData<4>{1.0e+08, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.74545e+04, 1.04528e+03},
    DensityData<4>{1.0e+05, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.48866e+04, 9.98562e+02},
    DensityData<4>{1.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.49088e+04, 9.98966e+02},
    DensityData<4>{5.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.50076e+04, 1.00076e+03},
    DensityData<4>{5.0e+07, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.61300e+04, 1.02118e+03},
    DensityData<4>{1.0e+08, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.74018e+04, 1.04432e+03},
    DensityData<4>{1.0e+05, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.33815e+04, 9.71180e+02},
    DensityData<4>{1.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.34034e+04, 9.71578e+02},
    DensityData<4>{5.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.35006e+04, 9.73346e+02},
    DensityData<4>{5.0e+07, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.46032e+04, 9.93407e+02},
    DensityData<4>{1.0e+08, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.58496e+04, 1.01608e+03},
    DensityData<4>{1.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.26439e+04, 9.57762e+02},
    DensityData<4>{5.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.27404e+04, 9.59516e+02},
    DensityData<4>{5.0e+07, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.38331e+04, 9.79396e+02},
    DensityData<4>{1.0e+08, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.50664e+04, 1.00183e+03},
  } )
);

INSTANTIATE_TEST_SUITE_P(
  CompositionalPhillipsBrineDensityTest, PengRobinson100K,
  ::testing::ValuesIn( {
    // Pure water: Mass density compared against CO2-Brine implementation
    DensityData<4>{1.0e+05, 288.15, {0.000, 1.000, 0.000, 0.000}, 5.48094e+04, 1.06077e+03},
    DensityData<4>{1.0e+06, 288.15, {0.000, 1.000, 0.000, 0.000}, 5.48347e+04, 1.06126e+03},
    DensityData<4>{5.0e+07, 288.15, {0.000, 1.000, 0.000, 0.000}, 5.62329e+04, 1.08832e+03},
    DensityData<4>{1.0e+05, 293.15, {0.000, 1.000, 0.000, 0.000}, 5.46468e+04, 1.05762e+03},
    DensityData<4>{1.0e+06, 293.15, {0.000, 1.000, 0.000, 0.000}, 5.46721e+04, 1.05811e+03},
    DensityData<4>{5.0e+07, 293.15, {0.000, 1.000, 0.000, 0.000}, 5.60653e+04, 1.08507e+03},
    DensityData<4>{1.0e+05, 373.15, {0.000, 1.000, 0.000, 0.000}, 5.15482e+04, 9.97652e+02}, // Undersaturated value
    DensityData<4>{1.0e+06, 373.15, {0.000, 1.000, 0.000, 0.000}, 5.19594e+04, 1.00561e+03},
    DensityData<4>{5.0e+07, 373.15, {0.000, 1.000, 0.000, 0.000}, 5.34550e+04, 1.03456e+03},
    // Mixtures
    DensityData<4>{1.0e+05, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.43082e+04, 1.06038e+03},
    DensityData<4>{1.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.43332e+04, 1.06087e+03},
    DensityData<4>{5.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.44439e+04, 1.06303e+03},
    DensityData<4>{5.0e+07, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.57127e+04, 1.08780e+03},
    DensityData<4>{1.0e+08, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.72372e+04, 1.11757e+03},
    DensityData<4>{1.0e+05, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.41443e+04, 1.05718e+03},
    DensityData<4>{1.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.41692e+04, 1.05767e+03},
    DensityData<4>{5.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.42801e+04, 1.05983e+03},
    DensityData<4>{5.0e+07, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.55444e+04, 1.08452e+03},
    DensityData<4>{1.0e+08, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.70525e+04, 1.11396e+03},
    DensityData<4>{1.0e+05, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.21182e+04, 1.01762e+03},
    DensityData<4>{1.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.21454e+04, 1.01815e+03},
    DensityData<4>{5.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.22658e+04, 1.02050e+03},
    DensityData<4>{5.0e+07, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.35742e+04, 1.04605e+03},
    DensityData<4>{1.0e+08, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.49938e+04, 1.07377e+03},
    DensityData<4>{1.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.14190e+04, 1.00397e+03},
    DensityData<4>{5.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.15461e+04, 1.00645e+03},
    DensityData<4>{5.0e+07, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.29090e+04, 1.03306e+03},
    DensityData<4>{1.0e+08, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.43435e+04, 1.06107e+03},
  } )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
