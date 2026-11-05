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
#include "constitutive/fluid/multifluid/compositional/functions/StabilityTest.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr std::array< real64, 4 > testPressures{1.0, 50.0, 100.0, 500.0}; // Pressure in bars
static constexpr std::array< real64, 4 > testTemperature{10.0, 25.0, 80.0, 150.0}; // Temperature in C

template< integer NC >
struct FluidData
{};

template<>
struct FluidData< 2 >
{
  static constexpr real64 e = 1.0e-6;

  static constexpr std::array< real64, 8 > testCompositions = {
    0.5, 0.5,
    0.0, 1.0,
    e, 1.0-e,
    1.0-e, e
  };

  static std::unique_ptr< TestFluid< 2 > > createFluid()
  {
    auto fluid = TestFluid< 2 >::create( {Fluid::CO2, Fluid::H2O} );
    fluid->setBinaryCoefficients( Feed< 1 >{0.1896} );
    return fluid;
  }
};

template<>
struct FluidData< 4 >
{
  static constexpr real64 e = 1.0e-6;

  static constexpr std::array< real64, 16 > testCompositions = {
    0.3, 0.2, 0.2, 0.3,
    0.0, 0.0, 0.0, 1.0,
    0.0, 0.5, 0.0, 0.5,
    e, e, e, 1.0-3*e
  };

  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    auto fluid = TestFluid< 4 >::create( {Fluid::CH4, Fluid::CO2, Fluid::H2S, Fluid::H2O} );
    const std::array< real64 const, 6 > bics = { 0.0, 0.0, 0.0, 0.4850, 0.1896, 0.1353 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template<>
struct FluidData< 6 >
{
  static constexpr real64 e = 1.0e-7;
  static constexpr real64 f = 0.16666;

  static constexpr std::array< real64, 30 > testCompositions = {
    f, f, f, f, f, 1.0-5*f,
    e, e, e, e, e, 1.0-5*e,
    e, 1.0-5*e, e, e, e, e,
    0.260871, 0.0663935, 0.0906286, 0.0524321, 0.52755293, 0.00212187,
    6.233e-08, 2.99791e-09, 8.00657e-06, 8.187700789e-05, 5.10942e-08, 0.99991
  };

  static std::unique_ptr< TestFluid< 6 > > createFluid()
  {
    auto fluid = TestFluid< 6 >::create( {Fluid::CH4, Fluid::C2H6, Fluid::CO2, Fluid::H2S, Fluid::N2, Fluid::H2O} );
    const std::array< real64 const, 15 > bics = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4850, 0.4920, 0.1896, 0.1353, 0.4778 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< integer NC >
using FlashData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const    // total composition
  >;

template< integer NC, EquationOfStateType EOS_TYPE >
class StabilityTestDerivativesTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr integer numComps = NC;
  static constexpr integer numDofs = NC + 2;

  using Deriv = geos::constitutive::multifluid::DerivativeOffset;

public:
  StabilityTestDerivativesTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {}

  ~StabilityTestDerivativesTestFixture() = default;

  void testStabilityTestDerivatives( FlashData< NC > const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    // Extract test point
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    StackArray< real64, 2, 4*numComps > workSpace( 4, numComps );
    StackArray< real64, 3, 2*numComps *numDofs > derivativeWorkSpace( 2, numComps, numDofs );
    arraySlice1d< real64 > incipientComposition = workSpace[0];
    arraySlice1d< real64 > kValues = workSpace[1];
    arraySlice1d< real64 > derivatives = workSpace[2];
    arraySlice1d< real64 > displacedIncipientComposition = workSpace[3];
    arraySlice2d< real64 > sampleDerivatives = derivativeWorkSpace[0];
    arraySlice2d< real64 > incipientDerivatives = derivativeWorkSpace[1];

    auto const setZero = []( real64 & value ){ value = 0.0; };

    LvArray::forValuesInSlice( workSpace.toSlice(), setZero );

    EquationOfStateType constexpr equationOfState = EOS_TYPE;

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = equationOfState;
    if constexpr (equationOfState == EquationOfStateType::SoreideWhitson)
    {
      flashData.vapourEos = EquationOfStateType::PengRobinson;
    }
    else
    {
      flashData.vapourEos = equationOfState;
    }

    // Extract column
    auto const extractColumn = [&]( integer const idof, real64 const scale = 1.0 )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        derivatives[ic] = scale * incipientDerivatives( ic, idof );
      }
    };

    if constexpr (equationOfState == EquationOfStateType::SoreideWhitson)
    {
      KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                          pressure,
                                                          temperature,
                                                          componentProperties,
                                                          kValues );
    }
    else
    {
      KValueInitialization::computeSoreideWhitsonKvalue( numComps,
                                                         pressure,
                                                         temperature,
                                                         componentProperties,
                                                         composition.toSliceConst(),
                                                         kValues );
    }

    auto const parameters = FlashParameters::create( std::make_unique< ModelParameters >() );
    auto const * flashParameters = parameters->get< FlashParameters >();
    flashParameters->m_discreteParameters[FlashParameters::STABILITY_MAX_ITERATIONS] = 20;

    bool unstableMixture = false;
    EquationOfStateType outputEos = EOS_TYPE;
    bool const stabilityStatus = StabilityTest::compute( numComps,
                                                         pressure,
                                                         temperature,
                                                         composition.toSliceConst(),
                                                         componentProperties,
                                                         flashData,
                                                         kValues.toSliceConst(),
                                                         flashParameters->m_continuousParameters,
                                                         flashParameters->m_discreteParameters,
                                                         unstableMixture,
                                                         outputEos,
                                                         incipientComposition );
    EXPECT_EQ( stabilityStatus, true );

    EquationOfStateType const sampleEos = (outputEos == flashData.vapourEos) ? flashData.liquidEos : flashData.vapourEos;
    StabilityTest::computeDerivatives( numComps,
                                       pressure,
                                       temperature,
                                       composition.toSliceConst(),
                                       componentProperties,
                                       sampleEos,
                                       outputEos,
                                       flashData,
                                       incipientComposition.toSliceConst(),
                                       incipientDerivatives,
                                       sampleDerivatives );

    // -- Pressure derivative
    extractColumn( Deriv::dP );
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative< numComps >(
      pressure, dp, derivatives,
      [&]( real64 const p, auto & values ){
      bool displacedUnstableMixture = false;
      EquationOfStateType displacedOutputEos = EOS_TYPE;
      StabilityTest::compute( numComps,
                              p,
                              temperature,
                              composition.toSliceConst(),
                              componentProperties,
                              flashData,
                              kValues.toSliceConst(),
                              flashParameters->m_continuousParameters,
                              flashParameters->m_discreteParameters,
                              displacedUnstableMixture,
                              displacedOutputEos,
                              displacedIncipientComposition );
      for( integer ic = 0; ic < numComps; ++ic )
      {
        values[ic] = displacedIncipientComposition[ic];
      }
    } );

    // -- Temperature derivative
    extractColumn( Deriv::dT );
    real64 const dT = 1.0e-5 * temperature;
    geos::testing::internal::testNumericalDerivative< numComps >(
      temperature, dT, derivatives,
      [&]( real64 const t, auto & values ){
      bool displacedUnstableMixture = false;
      EquationOfStateType displacedOutputEos = EOS_TYPE;
      StabilityTest::compute( numComps,
                              pressure,
                              t,
                              composition.toSliceConst(),
                              componentProperties,
                              flashData,
                              kValues.toSliceConst(),
                              flashParameters->m_continuousParameters,
                              flashParameters->m_discreteParameters,
                              displacedUnstableMixture,
                              displacedOutputEos,
                              displacedIncipientComposition );
      for( integer ic = 0; ic < numComps; ++ic )
      {
        values[ic] = displacedIncipientComposition[ic];
      }
    } );

    // -- Composition derivatives
    real64 constexpr dz = 1.0e-9;
    for( integer kc = 0; kc < numComps; ++kc )
    {
      // If the mixture is stable and the component doesn't exist then checking the derivative here will
      // not work because the incipient phase composition is not continuous in the component dimension.
      if( composition[kc] < MultiFluidConstants::minForSpeciesPresence && !unstableMixture )
      {
        continue;
      }
      extractColumn( Deriv::dC + kc );
      geos::testing::internal::testNumericalDerivative< numComps >(
        0.0, dz, derivatives,
        [&]( real64 const z, auto & values ){
        real64 const zOld = composition[kc];
        composition[kc] += z;
        bool displacedUnstableMixture = false;
        EquationOfStateType displacedOutputEos = EOS_TYPE;
        StabilityTest::compute( numComps,
                                pressure,
                                temperature,
                                composition.toSliceConst(),
                                componentProperties,
                                flashData,
                                kValues.toSliceConst(),
                                flashParameters->m_continuousParameters,
                                flashParameters->m_discreteParameters,
                                displacedUnstableMixture,
                                displacedOutputEos,
                                displacedIncipientComposition );

        for( integer ic = 0; ic < numComps; ++ic )
        {
          values[ic] = displacedIncipientComposition[ic];
        }
        composition[kc] = zOld;
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

using SoaveRedlichKwong2 = StabilityTestDerivativesTestFixture< 2, EquationOfStateType::SoaveRedlichKwong >;
using PengRobinson4 = StabilityTestDerivativesTestFixture< 4, EquationOfStateType::PengRobinson >;
using SoreideWhitson6 = StabilityTestDerivativesTestFixture< 6, EquationOfStateType::SoreideWhitson >;

TEST_P( SoaveRedlichKwong2, testDerivatives )
{
  testStabilityTestDerivatives( GetParam() );
}

TEST_P( PengRobinson4, testDerivatives )
{
  testStabilityTestDerivatives( GetParam() );
}

TEST_P( SoreideWhitson6, testDerivatives )
{
  testStabilityTestDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

template< integer NC >
stdVector< FlashData< NC > > generateData()
{
  constexpr integer dataSize = static_cast< integer >(FluidData< NC >::testCompositions.size());
  integer constexpr numSamples =  dataSize / NC;
  static_assert( numSamples * NC == dataSize );

  stdVector< FlashData< NC > > data;
  for( integer sampleIndex = 0; sampleIndex < numSamples; sampleIndex++ )
  {
    Feed< NC > sample;
    for( integer ic = 0; ic < NC; ic++ )
    {
      sample[ic] = FluidData< NC >::testCompositions[sampleIndex*NC + ic];
    }
    for( real64 const p : testPressures ) // Pressure in bars
    {
      for( real64 const t : testTemperature ) // Temperature in C
      {
        data.emplace_back( 1.0e5*p, t + 273.15, sample );
      }
    }
  }
  return data;
}

INSTANTIATE_TEST_SUITE_P( StabilityTestDerivatives, SoaveRedlichKwong2, ::testing::ValuesIn( generateData< 2 >()) );
INSTANTIATE_TEST_SUITE_P( StabilityTestDerivatives, PengRobinson4, ::testing::ValuesIn( generateData< 4 >()) );
INSTANTIATE_TEST_SUITE_P( StabilityTestDerivatives, SoreideWhitson6, ::testing::ValuesIn( generateData< 6 >()) );

} // testing
} // geos
