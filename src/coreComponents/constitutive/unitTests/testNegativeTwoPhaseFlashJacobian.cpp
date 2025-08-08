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
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/KValueInitialization.hpp"
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

template<>
struct FluidData< 9 >
{
static std::unique_ptr< TestFluid< 9 > > createFluid()
{
  std::unique_ptr< TestFluid< 9 > > fluid = TestFluid< 9 >::create( {0, 0, 0, 0, 0, 0, 0, 0, 0} );
  // Manually populate
  TestFluid< 9 >::populateArray( fluid->criticalPressure, Feed< 9 >{73.8659e5, 33.9439e5, 46.0421e5, 48.8387e5, 42.4552e5, 37.47e5, 33.5892e5, 30.1037e5, 20.549e5} );
  TestFluid< 9 >::populateArray( fluid->criticalTemperature, Feed< 9 >{304.7, 126.2, 190.6, 305.43, 369.8, 419.5, 465.9, 507.5, 678.8} );
  TestFluid< 9 >::populateArray( fluid->criticalVolume, Feed< 9 >{9.3999e-05, 9.0001e-05, 9.7999e-05, 1.4800e-04, 2.0000e-04, 2.5800e-04, 3.1000e-04, 3.5100e-04, 6.8243e-04} );
  TestFluid< 9 >::populateArray( fluid->acentricFactor, Feed< 9 >{0.225, 0.04, 0.013, 0.0986, 0.1524, 0.1956, 0.2413, 0.299, 0.5618} );
  TestFluid< 9 >::populateArray( fluid->molecularWeight, Feed< 9 >{44.01e-3, 28.01e-3, 16.04e-3, 30.07e-3, 44.1e-3, 58.12e-3, 72.15e-3, 84e-3, 173e-3} );
  TestFluid< 9 >::populateArray( fluid->volumeShift, Feed< 9 >{ -0.04958, -0.136012, -0.1486264, -0.10863408, -0.08349872, -0.06331568, -0.04196464, -0.0150072, 0.0000 } );
  fluid->setBinaryCoefficients( Feed< 36 >{
        1.0000e-02,
        0.0000e+00, 3.7320e-03,
        0.0000e+00, 1.0000e-02, 0.0000e+00,
        0.0000e+00, 1.0000e-02, 0.0000e+00, 0.0000e+00,
        0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
        0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
        1.0000e-02, 0.0000e+00, 2.8000e-02, 1.0000e-02, 1.0000e-02, 0.0000e+00, 0.0000e+00,
        1.0000e-02, 0.0000e+00, 4.5320e-02, 1.0000e-02, 1.0000e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00
      } );
  return fluid;
}
};

template< int NC >
using FlashData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // total composition
  real64 const        // vapour fraction
  >;

template< int NC, EquationOfStateType LIQUID_EOS_TYPE, EquationOfStateType VAPOUR_EOS_TYPE >
class NegativeTwoPhaseFlashJacobianTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  NegativeTwoPhaseFlashJacobianTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {}

  ~NegativeTwoPhaseFlashJacobianTestFixture() = default;

  void testJacobian( FlashData< NC > const & data )
  {
    // Number of values in each column
    constexpr integer numValues = 1 + numComps;

    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = LIQUID_EOS_TYPE;
    flashData.vapourEos = VAPOUR_EOS_TYPE;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 const vapourFraction = std::get< 3 >( data );

    StackArray< real64, 2, 5*numComps > workSpace( 5, numComps );
    StackArray< real64, 3, 2*numComps *numDofs > derivativeWorkSpace( 2, numComps, numDofs );
    StackArray< real64, 3, numValues *numValues > jacobianSpace( 1, numValues, numValues );
    StackArray< real64, 2, 3*numValues > residualSpace( 3, numValues );
    arraySlice1d< real64 > liquidComposition = workSpace[0];
    arraySlice1d< real64 > vapourComposition = workSpace[1];
    arraySlice1d< real64 > logLiquidFugacity = workSpace[2];
    arraySlice1d< real64 > logVapourFugacity = workSpace[3];
    arraySlice1d< real64 > kValues = workSpace[4];
    arraySlice2d< real64 > logLiquidFugacityDerivs = derivativeWorkSpace[0];
    arraySlice2d< real64 > logVapourFugacityDerivs = derivativeWorkSpace[1];
    arraySlice2d< real64 > jacobian = jacobianSpace[0];
    arraySlice1d< real64 > residual = residualSpace[0];
    arraySlice1d< real64 > displacedResidual = residualSpace[1];
    arraySlice1d< real64 > derivatives = residualSpace[2];

    // Extract column
    auto const extractColumn = [&]( integer const kc, real64 const scale = 1.0 )
    {
      for( integer ic = 0; ic < numValues; ++ic )
      {
        derivatives[ic] = scale * jacobian( ic, kc );
      }
    };

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kValues );

    std::array<real64,9> kv{ 5.50659394e-01, 6.30017745e+00, 2.74421602e+00, 4.80354013e-01, 1.56281500e-01, 5.21905616e-02, 1.78835233e-02, 5.75547635e-03, 9.74533928e-06 };
    for( integer ic = 0; ic < numComps; ++ic )
    {
      kValues[ic] = kv[ic];
    }


    NegativeTwoPhaseFlash::calculateResidualAndJacobian( numComps,
                                                         pressure,
                                                         temperature,
                                                         composition.toSliceConst(),
                                                         componentProperties,
                                                         flashData,
                                                         kValues.toSliceConst(),
                                                         vapourFraction,
                                                         liquidComposition,
                                                         vapourComposition,
                                                         logLiquidFugacity,
                                                         logVapourFugacity,
                                                         logLiquidFugacityDerivs,
                                                         logVapourFugacityDerivs,
                                                         residual,
                                                         jacobian );

    std::cout << "RESIDUAL " << std::scientific << std::setprecision( 5 ) << residual << "\n";
    std::cout << "JACOBIAN\n";
    for( integer kc = 0; kc <= numComps; ++kc )
    {
      std::cout << "  " << jacobian[kc] << "\n";
    }
    std::cout << "=============================================\n";

    // k-value derivatives
    // Derivatives are wrt to log(K)
    real64 constexpr deltaLogK = 1.0e-3;
    for( integer kc = 0; kc < numComps; ++kc )
    {
      extractColumn( kc );
      geos::testing::internal::testNumericalDerivative< numValues >(
        0.0, deltaLogK, derivatives,
        [&]( real64 const dK, auto & values )
      {
        real64 const kvOld = kValues[kc];
        kValues[kc] *= LvArray::math::exp( dK );
        NegativeTwoPhaseFlash::calculateResidual( numComps,
                                                  pressure,
                                                  temperature,
                                                  composition.toSliceConst(),
                                                  componentProperties,
                                                  flashData,
                                                  kValues.toSliceConst(),
                                                  vapourFraction,
                                                  liquidComposition,
                                                  vapourComposition,
                                                  logLiquidFugacity,
                                                  logVapourFugacity,
                                                  displacedResidual );
        for( integer ic = 0; ic < numValues; ++ic )
        {
          values[ic] = displacedResidual[ic];
        }
        kValues[kc] = kvOld;
      } );
    }

    // Vapour fraction derivative
    real64 constexpr deltaV = 1.0e-5;
    extractColumn( numComps );
    geos::testing::internal::testNumericalDerivative< numValues >(
      0.0, deltaV, derivatives,
      [&]( real64 const dV, auto & values )
    {
      NegativeTwoPhaseFlash::calculateResidual( numComps,
                                                pressure,
                                                temperature,
                                                composition.toSliceConst(),
                                                componentProperties,
                                                flashData,
                                                kValues.toSliceConst(),
                                                vapourFraction + dV,
                                                liquidComposition,
                                                vapourComposition,
                                                logLiquidFugacity,
                                                logVapourFugacity,
                                                displacedResidual );
      for( integer ic = 0; ic < numValues; ++ic )
      {
        values[ic] = displacedResidual[ic];
      }
    } );
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

//using PengRobinson4 = NegativeTwoPhaseFlashJacobianTestFixture< 4, EquationOfStateType::PengRobinson, EquationOfStateType::PengRobinson >;
//using SoreideWhitson4 = NegativeTwoPhaseFlashJacobianTestFixture< 4, EquationOfStateType::SoreideWhitson, EquationOfStateType::PengRobinson >;
using PengRobinson9 = NegativeTwoPhaseFlashJacobianTestFixture< 9, EquationOfStateType::PengRobinson, EquationOfStateType::PengRobinson >;

/*
TEST_P( PengRobinson4, testJacobian )
{
  testJacobian( GetParam() );
}

TEST_P( SoreideWhitson4, testJacobian )
{
  testJacobian( GetParam() );
}
*/
TEST_P( PengRobinson9, testJacobian )
{
  testJacobian( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

/*
INSTANTIATE_TEST_SUITE_P(NegativeTwoPhaseFlashJacobian, PengRobinson4,
  ::testing::ValuesIn<FlashData< 4 >>({
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.5}
  })
);

INSTANTIATE_TEST_SUITE_P(NegativeTwoPhaseFlashJacobian, SoreideWhitson4,
  ::testing::ValuesIn<FlashData< 4 >>({
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.0},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.5},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 1.0},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 0.998},
    {1.00000e+06, 193.15, {0.0, 0.0, 0.0, 1.0}, 1.002}
  })
);
*/

INSTANTIATE_TEST_SUITE_P(NegativeTwoPhaseFlashJacobian, PengRobinson9,
  ::testing::ValuesIn<FlashData< 9 >>({
    {1.00000e+07, 278.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, -0.18866442}
  })
);

/* UNCRUSTIFY-ON */

} // testing
} // geos
