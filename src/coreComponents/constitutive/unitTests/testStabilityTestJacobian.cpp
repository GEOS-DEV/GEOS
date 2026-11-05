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
#include "constitutive/fluid/multifluid/compositional/functions/KValueInitialization.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/Utilities.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
struct FluidData
{};

template<>
struct FluidData< 2 >
{
  static constexpr real64 e = 1.0e-8;

  static constexpr std::array< real64, 8 > testCompositions = {
    0.5, 0.5,
    0.0, 1.0,
    e, 1.0-e,
    1.0-e, e
  };

  static std::unique_ptr< TestFluid< 2 > > createFluid()
  {
    auto fluid = TestFluid< 2 >::create( {Fluid::CO2, Fluid::C5H12} );
    fluid->setBinaryCoefficients( Feed< 1 >{0.025} );
    return fluid;
  }
};

template<>
struct FluidData< 4 >
{
  static constexpr real64 e = 1.0e-7;

  static constexpr std::array< real64, 16 > testCompositions = {
    0.3, 0.2, 0.2, 0.3,
    0.0, 0.0, 0.0, 1.0,
    0.0, 0.5, 0.0, 0.5,
    e, e, 0.0, 1.0-2*e
  };

  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    auto fluid = TestFluid< 4 >::create( {Fluid::CH4, Fluid::CO2, Fluid::H2S, Fluid::H2O} );
    const std::array< real64 const, 6 > bics = { 0.0, 0.0, 0.0, 0.4850, 0.1896, 0.1353 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< integer NC >
using FlashData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // total composition
  real64 const        // power to raise k-values
  >;

template< integer NC, EquationOfStateType EOS_TYPE >
class StabilityTestJacobianTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr integer numComps = NC;
  static constexpr integer numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  StabilityTestJacobianTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {}

  ~StabilityTestJacobianTestFixture() = default;

  void testJacobian( FlashData< NC > const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

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

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 const alpha = std::get< 3 >( data );

    StackArray< real64, 2, 7*numComps > workSpace( 7, numComps );
    StackArray< real64, 3, numComps *numDofs > derivativeWorkSpace( 1, numComps, numDofs );
    StackArray< real64, 3, numComps *numComps > jacobianSpace( 1, numComps, numComps );
    StackArray< integer, 1, numComps > nonZeroComponents( numComps );
    arraySlice1d< real64 > testComposition = workSpace[0];
    arraySlice1d< real64 > logTestComposition = workSpace[1];
    arraySlice1d< real64 > testFugacity = workSpace[2];
    arraySlice1d< real64 > kValues = workSpace[3];
    arraySlice2d< real64 > logTestFugacityDerivs = derivativeWorkSpace[0];
    arraySlice2d< real64 > jacobian = jacobianSpace[0];
    arraySlice1d< real64 > residual = workSpace[4];
    arraySlice1d< real64 > displacedResidual = workSpace[5];
    arraySlice1d< real64 > derivatives = workSpace[6];

    auto const setZero = []( real64 & value ){ value = 0.0; };

    LvArray::forValuesInSlice( workSpace.toSlice(), setZero );

    // Extract column
    auto const extractColumn = [&]( integer const kc, real64 const scale = 1.0 )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        derivatives[ic] = scale * jacobian( ic, kc );
      }
    };

    calculatePresentComponents( numComps, composition.toSliceConst(), nonZeroComponents );
    auto presentComponents = nonZeroComponents.toSliceConst();

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

    for( integer const ic: presentComponents )
    {
      real64 const logK = LvArray::math::log( kValues[ic] );
      logTestComposition[ic] = LvArray::math::log( composition[ic] ) + alpha * logK;
    }

    real64 tpd = 0.0;

    StabilityTest::calculateResidualAndJacobian( numComps,
                                                 presentComponents,
                                                 pressure,
                                                 temperature,
                                                 componentProperties,
                                                 equationOfState,
                                                 flashData,
                                                 logTestComposition.toSliceConst(),
                                                 tpd,
                                                 testComposition,
                                                 testFugacity,
                                                 logTestFugacityDerivs,
                                                 residual,
                                                 jacobian );

    // composition derivatives
    // Derivatives are wrt to log(Y)
    real64 constexpr deltaLogY = 1.0e-3;
    for( integer kc = 0; kc < numComps; ++kc )
    {
      extractColumn( kc );
      geos::testing::internal::testNumericalDerivative< numComps >(
        0.0, deltaLogY, derivatives,
        [&]( real64 const dY, auto & values )
      {
        LvArray::forValuesInSlice( displacedResidual, setZero );
        real64 const logYOld = logTestComposition[kc];
        logTestComposition[kc] += dY;
        StabilityTest::calculateResidual( numComps,
                                          presentComponents,
                                          pressure,
                                          temperature,
                                          componentProperties,
                                          equationOfState,
                                          flashData,
                                          logTestComposition.toSliceConst(),
                                          tpd,
                                          testComposition,
                                          testFugacity,
                                          displacedResidual );
        for( integer ic = 0; ic < numComps; ++ic )
        {
          values[ic] = displacedResidual[ic];
        }
        logTestComposition[kc] = logYOld;
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

using PengRobinson4 = StabilityTestJacobianTestFixture< 4, EquationOfStateType::PengRobinson >;
using SoreideWhitson4 = StabilityTestJacobianTestFixture< 4, EquationOfStateType::SoreideWhitson >;
using SoaveRedlichKwong2 = StabilityTestJacobianTestFixture< 2, EquationOfStateType::SoaveRedlichKwong >;

TEST_P( PengRobinson4, testJacobian )
{
  testJacobian( GetParam() );
}

TEST_P( SoreideWhitson4, testJacobian )
{
  testJacobian( GetParam() );
}

TEST_P( SoaveRedlichKwong2, testJacobian )
{
  testJacobian( GetParam() );
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
    for( real64 const p : {1.0, 50.0, 100.0, 500.0} ) // Pressure in bars
    {
      for( real64 const t : {10.0, 25.0, 80.0, 150.0} ) // Temperature in C
      {
        for( real64 const a : {-1.0, 0.0, 1.0} )
        {
          data.emplace_back( 1.0e5*p, t + 273.15, sample, a );
        }
      }
    }
  }
  return data;
}

INSTANTIATE_TEST_SUITE_P( StabilityTestJacobian, PengRobinson4, ::testing::ValuesIn( generateData< 4 >()) );
INSTANTIATE_TEST_SUITE_P( StabilityTestJacobian, SoreideWhitson4, ::testing::ValuesIn( generateData< 4 >()) );
INSTANTIATE_TEST_SUITE_P( StabilityTestJacobian, SoaveRedlichKwong2, ::testing::ValuesIn( generateData< 2 >()) );

} // testing
} // geos
