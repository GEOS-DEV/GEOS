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
#include "constitutive/fluid/multifluid/compositional/models/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PressureTemperatureCoordinates.hpp"

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
    auto kernelWrapper = this->m_density->createKernelWrapper();

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
    stackArray1d< real64, numDofs > molarDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           phaseComposition.toSliceConst(),
                           molarDensity,
                           molarDensityDerivs.toSlice(),
                           massDensity,
                           massDensityDerivs.toSlice(),
                           false );

    real64 constexpr scale = 1.0;

    molarDensity *= scale;
    massDensity *= scale;
    LvArray::forValuesInSlice( molarDensityDerivs.toSlice(), []( real64 & a ){ a *= scale; } );
    LvArray::forValuesInSlice( massDensityDerivs.toSlice(), []( real64 & a ){ a *= scale; } );

    auto calculateDensity = [&]( real64 const p, real64 const t, auto const & zmf ) -> std::pair< real64, real64 > {
      real64 densityMolar = 0.0;
      real64 densityMass = 0.0;
      stackArray1d< real64, numDofs > tempDerivs( numDofs );
      kernelWrapper.compute( componentProperties, p, t,
                             zmf.toSliceConst(),
                             densityMolar,
                             tempDerivs.toSlice(),
                             densityMass,
                             tempDerivs.toSlice(),
                             false );
      return {scale * densityMolar, scale * densityMass};
    };

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, molarDensityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return calculateDensity( p, temperature, phaseComposition ).first;
    } );
    internal::testNumericalDerivative(
      pressure, dp, massDensityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return calculateDensity( p, temperature, phaseComposition ).second;
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, molarDensityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return calculateDensity( pressure, t, phaseComposition ).first;
    } );
    internal::testNumericalDerivative(
      temperature, dT, massDensityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return calculateDensity( pressure, t, phaseComposition ).second;
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-6;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, molarDensityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = phaseComposition[jc];
        }
        zmf[ic] += z;
        return calculateDensity( pressure, temperature, zmf ).first;
      } );
      internal::testNumericalDerivative(
        0.0, dz, massDensityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = phaseComposition[jc];
        }
        zmf[ic] += z;
        return calculateDensity( pressure, temperature, zmf ).second;
      } );
    }
  }

private:
  static void createPressurePoints( array1d< real64 > & pressureCoordinates )
  {
    int constexpr n = 20;
    real64 constexpr minPressure = 0.995e5;
    real64 constexpr maxPressure = 1000.005e5;
    real64 constexpr r = pow( maxPressure/minPressure, 1.0 / n );
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

INSTANTIATE_TEST_SUITE_P(
  CompositionalPhillipsBrineDensityTest, PengRobinson,
  ::testing::ValuesIn( {
    DensityData<4>{1.0e+05, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.5944615e+04, 1.0178113e+03},
    DensityData<4>{1.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.5973705e+04, 1.0183406e+03},
    DensityData<4>{5.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.6102402e+04, 1.0206820e+03},
    DensityData<4>{5.0e+07, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.7501010e+04, 1.0461271e+03},
    DensityData<4>{1.0e+08, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.9016111e+04, 1.0736917e+03},
    DensityData<4>{1.0e+05, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.5752244e+04, 1.0143115e+03},
    DensityData<4>{1.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.5781710e+04, 1.0148476e+03},
    DensityData<4>{5.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.5912023e+04, 1.0172184e+03},
    DensityData<4>{5.0e+07, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.7323318e+04, 1.0428944e+03},
    DensityData<4>{1.0e+08, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.8840359e+04, 1.0704942e+03},
    DensityData<4>{1.0e+05, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.3210272e+04, 9.6806491e+02},
    DensityData<4>{1.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.3246283e+04, 9.6872005e+02},
    DensityData<4>{5.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.3405091e+04, 9.7160928e+02},
    DensityData<4>{5.0e+07, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.5071382e+04, 1.0019244e+03},
    DensityData<4>{1.0e+08, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.6727179e+04, 1.0320487e+03},
    DensityData<4>{1.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.2276940e+04, 9.5108462e+02},
    DensityData<4>{5.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.2448976e+04, 9.5421451e+02},
    DensityData<4>{5.0e+07, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.4242231e+04, 9.8683954e+02},
    DensityData<4>{1.0e+08, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.5991229e+04, 1.0186594e+03}
  } )
);

INSTANTIATE_TEST_SUITE_P(
  CompositionalPhillipsBrineDensityTest, PengRobinson100K,
  ::testing::ValuesIn( {
    DensityData<4>{1.0e+05, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.4308247e+04, 1.0603832e+03},
    DensityData<4>{1.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.4333160e+04, 1.0608696e+03},
    DensityData<4>{5.0e+06, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.4443949e+04, 1.0630328e+03},
    DensityData<4>{5.0e+07, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.5712671e+04, 1.0878050e+03},
    DensityData<4>{1.0e+08, 288.15, {0.003, 0.995, 0.005, 0.002}, 5.7237189e+04, 1.1175716e+03},
    DensityData<4>{1.0e+05, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.4144306e+04, 1.0571822e+03},
    DensityData<4>{1.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.4169246e+04, 1.0576692e+03},
    DensityData<4>{5.0e+06, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.4280108e+04, 1.0598338e+03},
    DensityData<4>{5.0e+07, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.5544412e+04, 1.0845197e+03},
    DensityData<4>{1.0e+08, 293.15, {0.003, 0.995, 0.005, 0.002}, 5.7052472e+04, 1.1139649e+03},
    DensityData<4>{1.0e+05, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.2118214e+04, 1.0176222e+03},
    DensityData<4>{1.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.2145414e+04, 1.0181533e+03},
    DensityData<4>{5.0e+06, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.2265760e+04, 1.0205031e+03},
    DensityData<4>{5.0e+07, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.3574220e+04, 1.0460511e+03},
    DensityData<4>{1.0e+08, 353.15, {0.003, 0.995, 0.005, 0.002}, 5.4993761e+04, 1.0737680e+03},
    DensityData<4>{1.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.1419016e+04, 1.0039702e+03},
    DensityData<4>{5.0e+06, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.1546068e+04, 1.0064509e+03},
    DensityData<4>{5.0e+07, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.2909038e+04, 1.0330633e+03},
    DensityData<4>{1.0e+08, 373.15, {0.003, 0.995, 0.005, 0.002}, 5.4343460e+04, 1.0610708e+03}
  } )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
