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
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
struct FluidData
{};

template<>
struct FluidData< 2 >
{
  static std::unique_ptr< TestFluid< 2 > > createFluid()
  {
    auto fluid = TestFluid< 2 >::create( {Fluid::CO2, Fluid::C5H12} );
    return fluid;
  }
};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    return TestFluid< 4 >::create( {Fluid::N2, Fluid::CH4, Fluid::CO2, Fluid::H2O} );
  }
};

template< int NC >
using FlashData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // total composition
  real64 const        // power to raise k-values
  >;

template< int NC, EquationOfStateType EOS_TYPE >
class StabilityTestJacobianTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  StabilityTestJacobianTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {}

  ~StabilityTestJacobianTestFixture() = default;

  void testJacobian( FlashData< NC > const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    EquationOfStateType const equationOfState = EOS_TYPE;

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = equationOfState;
    flashData.vapourEos = equationOfState;

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

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kValues );
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

    std::cout << "RESIDUAL " << std::scientific << std::setprecision( 5 ) << residual << "\n";
    std::cout << "RESIDUAL " << std::scientific << std::setprecision( 5 ) << displacedResidual << "\n";
    std::cout << "JACOBIAN\n";
    for( integer kc = 0; kc < numComps; ++kc )
    {
      std::cout << "  " << jacobian[kc] << "\n";
    }
    std::cout << "=============================================\n";

    // k-value derivatives
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
        real64 const kvOld = logTestComposition[kc];
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
        logTestComposition[kc] = kvOld;
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

using PengRobinson4 = StabilityTestJacobianTestFixture< 4, EquationOfStateType::PengRobinson >;
//using SoreideWhitson4 = StabilityTestJacobianTestFixture< 4, EquationOfStateType::SoreideWhitson >;

TEST_P( PengRobinson4, testJacobian )
{
  testJacobian( GetParam() );
}

//TEST_P( SoreideWhitson4, testJacobian )
//{
//  testJacobian( GetParam() );
//}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(StabilityTestJacobian, PengRobinson4,
  ::testing::ValuesIn<FlashData< 4 >>({
    {1.0e+07, 193.15, {0.3, 0.2, 0.2, 0.3},  0.0}
//    {1.0e+07, 193.15, {0.3, 0.2, 0.2, 0.3},  1.0},
//    {1.0e+07, 193.15, {0.3, 0.2, 0.2, 0.3}, -1.0},
//    {1.0e+07, 193.15, {0.0, 0.2, 0.0, 0.8}, -1.0}
  })
);

/*
INSTANTIATE_TEST_SUITE_P(StabilityTestJacobian, SoreideWhitson4,
  ::testing::ValuesIn<FlashData< 4 >>({
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.0},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.5},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 1.0},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.998},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 1.002}
  })
);
*/

/* UNCRUSTIFY-ON */

} // testing
} // geos
