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
    {
      std::vector< real64 > x{0.03512106, 0.10723110, 0.10724847, 0.75039936};
      std::vector< real64 > y{0.99661534, 0.00094661, 0.00037503, 0.00206303};
      for( integer kc = 0; kc < numComps; ++kc )
      {
        kValues[kc] = y[kc]/x[kc];
      }
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

using PengRobinson4 = NegativeTwoPhaseFlashJacobianTestFixture< 4, EquationOfStateType::PengRobinson, EquationOfStateType::PengRobinson >;
using SoreideWhitson4 = NegativeTwoPhaseFlashJacobianTestFixture< 4, EquationOfStateType::SoreideWhitson, EquationOfStateType::SoreideWhitson >;

TEST_P( PengRobinson4, testJacobian )
{
  testJacobian( GetParam() );
}

TEST_P( SoreideWhitson4, testJacobian )
{
  testJacobian( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

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

/* UNCRUSTIFY-ON */

} // testing
} // geos
