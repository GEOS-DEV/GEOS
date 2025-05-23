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
#include "constitutive/fluid/multifluid/compositional/models/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"
#include "functions/FunctionManager.hpp"

#include <conduit.hpp>

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using ViscosityData = std::tuple<
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< NC > const,         // phase composition
  real64 const              // expected viscosity
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 2 >
{
  static std::unique_ptr< TestFluid< 2 > > createFluid()
  {
    auto fluid = TestFluid< 2 >::create( {Fluid::H2O, Fluid::CO2} );
    return fluid;
  }
};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    auto fluid = TestFluid< 4 >::create( {Fluid::H2O, Fluid::CO2, Fluid::C1, Fluid::C5} );
    return fluid;
  }
};

template< int NC, int SALINITY = 0 >
class PhillipsBrineViscosityTestFixture :  public ::testing::TestWithParam< ViscosityData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  PhillipsBrineViscosityTestFixture()
    : m_parent( "parent", m_node ),
    m_fluid( FluidData< NC >::createFluid() )
  {
    m_functionManager = std::make_unique< FunctionManager >( "FunctionManager", &m_parent );

    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = PhillipsBrineViscosity::createParameters( std::move( m_parameters ) );

    auto coordinates = const_cast< PressureTemperatureCoordinates * >(m_parameters->get< PressureTemperatureCoordinates >());
    createTemperaturePoints( coordinates->m_temperatureCoordinates );

    auto brineSalinity = const_cast< BrineSalinity * >(m_parameters->get< BrineSalinity >());
    real64 const massFraction = 1.0e-6*SALINITY;
    brineSalinity->m_salinity = massFraction / 58.44e-3;

    m_viscosity = std::make_unique< PhillipsBrineViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~PhillipsBrineViscosityTestFixture() = default;

  void testViscosity( ViscosityData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));
    real64 const expectedViscosity = std::get< 3 >( data );

    real64 const massDensity = 1000.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    viscosityKernelWrapper.compute( componentProperties,
                                    pressure,
                                    temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    tempDerivs.toSliceConst(),
                                    viscosity,
                                    tempDerivs.toSlice(),
                                    false );

    checkRelativeError( viscosity, expectedViscosity, relTol, absTol );
  }

  void testViscosityDerivatives( ViscosityData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));

    auto componentProperties = m_fluid->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    real64 const massDensity = 1000.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > viscosityDerivs( numDofs );
    LvArray::forValuesInSlice( massDensityDerivs.toSlice(), []( real64 & a ){ a = 0.0; } );

    viscosityKernelWrapper.compute( componentProperties,
                                    pressure,
                                    temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    massDensityDerivs.toSliceConst(),
                                    viscosity,
                                    viscosityDerivs.toSlice(),
                                    false );

    auto calculateViscosity = [&]( real64 const p, real64 const t, auto const & zmf ) -> real64 {
      real64 phaseViscosity = 0.0;
      stackArray1d< real64, numDofs > tempDerivs( numDofs );
      viscosityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                      massDensity, tempDerivs.toSliceConst(), phaseViscosity, tempDerivs.toSlice(), false );
      return phaseViscosity;
    };

    // Viscosity values are very small so we will inflate the values to avoid false positives due
    // to the absolute value check
    real64 constexpr scale = 1.0e6;

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, scale*viscosityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return scale*calculateViscosity( p, temperature, phaseComposition );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, scale*viscosityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return scale*calculateViscosity( pressure, t, phaseComposition );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < 1; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, scale*viscosityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = phaseComposition[jc];
        }
        zmf[ic] += z;
        return scale*calculateViscosity( pressure, temperature, zmf );
      } );
    }
  }

private:
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
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< PhillipsBrineViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using PhillipsBrineViscosity4 = PhillipsBrineViscosityTestFixture< 4 >;
using PhillipsBrineViscosity2_100K = PhillipsBrineViscosityTestFixture< 2, 100000 >;

TEST_P( PhillipsBrineViscosity4, testViscosity )
{
  testViscosity( GetParam() );
}

TEST_P( PhillipsBrineViscosity2_100K, testViscosity )
{
  testViscosity( GetParam() );
}

TEST_P( PhillipsBrineViscosity4, testViscosityDerivatives )
{
  testViscosityDerivatives( GetParam() );
}

TEST_P( PhillipsBrineViscosity2_100K, testViscosityDerivatives )
{
  testViscosityDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  PhillipsBrineViscosity,
  PhillipsBrineViscosity4,
  ::testing::ValuesIn({
    ViscosityData<4>{ 1.0e+07, 288.15, {1.0000, 0.0000, 0.0000, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {1.0000, 0.0000, 0.0000, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {1.0000, 0.0000, 0.0000, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {1.0000, 0.0000, 0.0000, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9990, 0.0000, 0.0010, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9990, 0.0000, 0.0010, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9990, 0.0000, 0.0010, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9990, 0.0000, 0.0010, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9900, 0.0000, 0.0100, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9900, 0.0000, 0.0100, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9900, 0.0000, 0.0100, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9900, 0.0000, 0.0100, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9990, 0.0010, 0.0000, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9990, 0.0010, 0.0000, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9990, 0.0010, 0.0000, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9990, 0.0010, 0.0000, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9980, 0.0010, 0.0010, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9980, 0.0010, 0.0010, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9980, 0.0010, 0.0010, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9980, 0.0010, 0.0010, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9890, 0.0010, 0.0100, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9890, 0.0010, 0.0100, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9890, 0.0010, 0.0100, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9890, 0.0010, 0.0100, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9900, 0.0100, 0.0000, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9900, 0.0100, 0.0000, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9900, 0.0100, 0.0000, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9900, 0.0100, 0.0000, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9890, 0.0100, 0.0010, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9890, 0.0100, 0.0010, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9890, 0.0100, 0.0010, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9890, 0.0100, 0.0010, 0.0000}, 2.822282e-04 },
    ViscosityData<4>{ 1.0e+07, 288.15, {0.9800, 0.0100, 0.0100, 0.0000}, 1.153800e-03 },
    ViscosityData<4>{ 1.0e+07, 293.15, {0.9800, 0.0100, 0.0100, 0.0000}, 1.012544e-03 },
    ViscosityData<4>{ 1.0e+07, 353.15, {0.9800, 0.0100, 0.0100, 0.0000}, 3.552664e-04 },
    ViscosityData<4>{ 1.0e+07, 373.15, {0.9800, 0.0100, 0.0100, 0.0000}, 2.822282e-04 }
  })
);

INSTANTIATE_TEST_SUITE_P(
  PhillipsBrineViscosity,
  PhillipsBrineViscosity2_100K,
  ::testing::ValuesIn({
    ViscosityData<2>{ 1.0e+07, 288.15, {1.0000, 0.0000}, 1.364366e-03 },
    ViscosityData<2>{ 1.0e+07, 293.15, {1.0000, 0.0000}, 1.199555e-03 },
    ViscosityData<2>{ 1.0e+07, 353.15, {1.0000, 0.0000}, 4.302578e-04 },
    ViscosityData<2>{ 1.0e+07, 373.15, {1.0000, 0.0000}, 3.442847e-04 },
    ViscosityData<2>{ 1.0e+07, 288.15, {0.9990, 0.0010}, 1.364366e-03 },
    ViscosityData<2>{ 1.0e+07, 293.15, {0.9990, 0.0010}, 1.199555e-03 },
    ViscosityData<2>{ 1.0e+07, 353.15, {0.9990, 0.0010}, 4.302578e-04 },
    ViscosityData<2>{ 1.0e+07, 373.15, {0.9990, 0.0010}, 3.442847e-04 },
    ViscosityData<2>{ 1.0e+07, 288.15, {0.9900, 0.0100}, 1.364366e-03 },
    ViscosityData<2>{ 1.0e+07, 293.15, {0.9900, 0.0100}, 1.199555e-03 },
    ViscosityData<2>{ 1.0e+07, 353.15, {0.9900, 0.0100}, 4.302578e-04 },
    ViscosityData<2>{ 1.0e+07, 373.15, {0.9900, 0.0100}, 3.442847e-04 }
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
