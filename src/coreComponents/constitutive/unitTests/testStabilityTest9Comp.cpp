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
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/StabilityTest.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr integer numComps = 9;
static constexpr integer targetComponent = 8;

using StabilityData = std::tuple<
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< numComps > const,   // total composition
  bool const,               // expected instability status of mixture
  real64 const              // expected mole fraction of targetComponent in incipient phase
  >;

template< EquationOfStateType EOS_TYPE >
class StabilityTestTest9CompFixture :  public ::testing::TestWithParam< StabilityData >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
public:
  StabilityTestTest9CompFixture()
    : m_fluid( createFluid() )
  {}

  ~StabilityTestTest9CompFixture() = default;

  void testStability( StabilityData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EOS_TYPE;
    flashData.vapourEos = EOS_TYPE;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

    bool const expectedUnstableMixture = std::get< 3 >( data );
    real64 const expectedZTarget = std::get< 4 >( data );

    bool unstableMixture = false;
    EquationOfStateType incipientEos = EOS_TYPE;

    stackArray1d< real64, numComps > kValues( numComps );
    stackArray1d< real64, numComps > incipientComposition( numComps );

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kValues.toSlice() );

    auto const parameters = FlashParameters::create( std::make_unique< ModelParameters >() );
    auto const * flashParameters = parameters->get< FlashParameters >();
    setMaximumNumberOfIterations( flashParameters, 20 );

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
                                                         incipientEos,
                                                         incipientComposition.toSlice() );

    // Expect this to succeed
    ASSERT_EQ( stabilityStatus, true );

    // Check the stability label
    ASSERT_EQ( unstableMixture, expectedUnstableMixture );

    // Check the incipient methane fraction
    real64 const incipientZTarget = incipientComposition[targetComponent];
    checkRelativeError( expectedZTarget, incipientZTarget, relTol, absTol );
  }

  void testFlash( StabilityData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EOS_TYPE;
    flashData.vapourEos = EOS_TYPE;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

    bool unstableMixture = false;
    EquationOfStateType incipientEos = EOS_TYPE;

    StackArray< real64, 2, 3*numComps > phaseCompositions( 3, numComps );
    real64 vapourFraction = -1.0;
    stackArray2d< real64, numComps > kValues( 1, numComps );
    arraySlice1d< real64 > liquidComposition = phaseCompositions[0];
    arraySlice1d< real64 > vapourComposition = phaseCompositions[1];
    arraySlice1d< real64 > incipientComposition = phaseCompositions[2];

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kValues[0] );

    auto const parameters = FlashParameters::create( std::make_unique< ModelParameters >() );
    auto const * flashParameters = parameters->get< FlashParameters >();
    setMaximumNumberOfIterations( flashParameters, 20 );

    StabilityTest::compute( numComps,
                            pressure,
                            temperature,
                            composition.toSliceConst(),
                            componentProperties,
                            flashData,
                            kValues[0].toSliceConst(),
                            flashParameters->m_continuousParameters,
                            flashParameters->m_discreteParameters,
                            unstableMixture,
                            incipientEos,
                            incipientComposition );

    // Now perform the nagative flash
    bool const status = NegativeTwoPhaseFlash::compute(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      flashData,
      flashParameters->m_continuousParameters,
      flashParameters->m_discreteParameters,
      kValues.toSlice(),
      vapourFraction,
      liquidComposition,
      vapourComposition );

    // Expect the flash to succeed for an unstable mixture
    ASSERT_EQ( status || !unstableMixture, true );

    if( status )
    {
      constexpr real64 epsilon = MultiFluidConstants::epsilon;
      real64 const V = vapourFraction;
      real64 const L = 1.0 - V;

      // Check stability vs flash
      if( unstableMixture )
      {
        // Unstable mixture: both L and V should be positive
        ASSERT_GT( V, epsilon );
        ASSERT_GT( L, epsilon );
      }
      else
      {
        // Stable mixture: one of L or V should be zero
        ASSERT_TRUE( V < epsilon || L < epsilon );
      }
    }
  }

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
private:
  static void setMaximumNumberOfIterations( FlashParameters const * parameters, integer const maxIterations );
  static std::unique_ptr< TestFluid< numComps > > createFluid();
};

template< EquationOfStateType EOS_TYPE >
void StabilityTestTest9CompFixture< EOS_TYPE >::setMaximumNumberOfIterations( FlashParameters const * parameters,
                                                                              integer const maxIterations )
{
  auto * flashParameters = const_cast< FlashParameters * >(parameters);
  flashParameters->m_discreteParameters[FlashParameters::STABILITY_MAX_ITERATIONS] = maxIterations;
}

template< EquationOfStateType EOS_TYPE >
std::unique_ptr< TestFluid< numComps > > StabilityTestTest9CompFixture< EOS_TYPE >::createFluid()
{
  std::unique_ptr< TestFluid< numComps > > fluid = TestFluid< numComps >::create( {0, 0, 0, 0, 0, 0, 0, 0, 0} );
  // Manually populate
  TestFluid< numComps >::populateArray( fluid->criticalPressure, Feed< 9 >{73.8659e5, 33.9439e5, 46.0421e5, 48.8387e5, 42.4552e5, 37.47e5, 33.5892e5, 30.1037e5, 20.549e5} );
  TestFluid< numComps >::populateArray( fluid->criticalTemperature, Feed< 9 >{304.7, 126.2, 190.6, 305.43, 369.8, 419.5, 465.9, 507.5, 678.8} );
  TestFluid< numComps >::populateArray( fluid->criticalVolume, Feed< 9 >{9.3999e-05, 9.0001e-05, 9.7999e-05, 1.4800e-04, 2.0000e-04, 2.5800e-04, 3.1000e-04, 3.5100e-04, 6.8243e-04} );
  TestFluid< numComps >::populateArray( fluid->acentricFactor, Feed< 9 >{0.225, 0.04, 0.013, 0.0986, 0.1524, 0.1956, 0.2413, 0.299, 0.5618} );
  TestFluid< numComps >::populateArray( fluid->molecularWeight, Feed< 9 >{44.01e-3, 28.01e-3, 16.04e-3, 30.07e-3, 44.1e-3, 58.12e-3, 72.15e-3, 84e-3, 173e-3} );
  TestFluid< numComps >::populateArray( fluid->volumeShift, Feed< 9 >{ -0.04958, -0.136012, -0.1486264, -0.10863408, -0.08349872, -0.06331568, -0.04196464, -0.0150072, 0.0000 } );
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

using PengRobinson = StabilityTestTest9CompFixture< EquationOfStateType::PengRobinson >;
using SoaveRedlichKwong = StabilityTestTest9CompFixture< EquationOfStateType::SoaveRedlichKwong >;

TEST_P( PengRobinson, testStabilityTest )
{
  testStability( GetParam() );
}
TEST_P( SoaveRedlichKwong, testStabilityTest )
{
  testStability( GetParam() );
}
TEST_P( PengRobinson, testStabilityWithFlash )
{
  testFlash( GetParam() );
}
TEST_P( SoaveRedlichKwong, testStabilityWithFlash )
{
  testFlash( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, PengRobinson,
  ::testing::ValuesIn<StabilityData>({
    {1.0e+05, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  9.99930e-01},
    {1.0e+05, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  9.99889e-01},
    {1.0e+05, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  9.99068e-01},
    {1.0e+05, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {5.0e+06, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  3.31496e-06},
    {5.0e+06, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  6.19514e-06},
    {5.0e+06, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  1.68840e-04},
    {5.0e+06, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  7.60017e-01},
    {1.0e+07, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  4.69115e-05},
    {1.0e+07, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  6.37278e-05},
    {1.0e+07, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  4.72565e-04},
    {1.0e+07, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  5.69337e-01},
    {1.0e+08, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+08, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+08, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+08, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+05, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  9.99931e-01},
    {1.0e+05, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  9.99891e-01},
    {1.0e+05, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  9.99113e-01},
    {1.0e+05, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {5.0e+06, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  3.10425e-06},
    {5.0e+06, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  5.93265e-06},
    {5.0e+06, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  2.11049e-04},
    {5.0e+06, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  7.78486e-01},
    {1.0e+07, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  3.30848e-05},
    {1.0e+07, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  4.73576e-05},
    {1.0e+07, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  4.31414e-04},
    {1.0e+07, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  6.09508e-01},
    {1.0e+08, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+08, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+08, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+08, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+05, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  9.99380e-01},
    {1.0e+05, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  9.99048e-01},
    {1.0e+05, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  9.93558e-01},
    {1.0e+05, 573.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {5.0e+06, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  1.20037e-06},
    {5.0e+06, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  2.41217e-06},
    {5.0e+06, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  6.49045e-05},
    {5.0e+06, 573.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+07, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.43584e-05},
    {1.0e+07, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  2.41543e-05},
    {1.0e+07, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  3.04002e-04},
    {1.0e+08, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+08, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+08, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+08, 573.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+05, 573.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.0000000e+00}, false, 1.00000e+00},
    {1.0e+08, 288.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.0000000e+00}, false, 1.00000e+00},
    {1.0e+05, 573.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 1.3300000e-05, 9.9997340e-01}, false, 9.99973e-01},
    {1.0e+08, 288.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 1.3300000e-05, 9.9997340e-01}, false, 9.99973e-01},
    {1.0e+05, 573.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 9.9997340e-01, 1.3300000e-05}, false, 1.33000e-05},
    {1.0e+08, 288.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 9.9997340e-01, 1.3300000e-05}, false, 1.33000e-05}
  })
);

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, SoaveRedlichKwong,
  ::testing::ValuesIn<StabilityData>({
    {1.0e+05, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  9.99943e-01},
    {1.0e+05, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  9.99907e-01},
    {1.0e+05, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  9.99162e-01},
    {1.0e+05, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {5.0e+06, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  2.35604e-06},
    {5.0e+06, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  4.52623e-06},
    {5.0e+06, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  1.34208e-04},
    {5.0e+06, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  7.61380e-01},
    {1.0e+07, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  3.24422e-05},
    {1.0e+07, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  4.49843e-05},
    {1.0e+07, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  3.62686e-04},
    {1.0e+07, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, true,  5.80899e-01},
    {1.0e+08, 288.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+08, 298.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+08, 353.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+08, 573.15, {9.0000000e-03, 3.0000000e-03, 5.3470000e-01, 1.1460000e-01, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6920000e-01}, false, 1.69200e-01},
    {1.0e+05, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  9.99944e-01},
    {1.0e+05, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  9.99909e-01},
    {1.0e+05, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  9.99202e-01},
    {1.0e+05, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {5.0e+06, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  2.21223e-06},
    {5.0e+06, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  4.33855e-06},
    {5.0e+06, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  1.63405e-04},
    {5.0e+06, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  7.79335e-01},
    {1.0e+07, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  2.28707e-05},
    {1.0e+07, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  3.33621e-05},
    {1.0e+07, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  3.27432e-04},
    {1.0e+07, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, true,  6.18577e-01},
    {1.0e+08, 288.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+08, 298.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+08, 353.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+08, 573.15, {9.0000000e-03, 0.0000000e+00, 6.5230000e-01, 1.3300000e-05, 8.7900000e-02, 4.5600000e-02, 2.0900000e-02, 1.5100000e-02, 1.6918670e-01}, false, 1.69187e-01},
    {1.0e+05, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  9.99489e-01},
    {1.0e+05, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  9.99201e-01},
    {1.0e+05, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  9.94186e-01},
    {1.0e+05, 573.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {5.0e+06, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  8.06893e-07},
    {5.0e+06, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  1.68319e-06},
    {5.0e+06, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  5.24670e-05},
    {5.0e+06, 573.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+07, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  8.46960e-06},
    {1.0e+07, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  1.49116e-05},
    {1.0e+07, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, true,  2.25256e-04},
    {1.0e+08, 288.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+08, 298.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+08, 353.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+08, 573.15, {1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111100e-01, 1.1111200e-01}, false, 1.11112e-01},
    {1.0e+05, 573.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.0000000e+00}, false, 1.00000e+00},
    {1.0e+08, 288.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.0000000e+00}, false, 1.00000e+00},
    {1.0e+05, 573.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 1.3300000e-05, 9.9997340e-01}, false, 9.99973e-01},
    {1.0e+08, 288.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 1.3300000e-05, 9.9997340e-01}, false, 9.99973e-01},
    {1.0e+05, 573.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 9.9997340e-01, 1.3300000e-05}, false, 1.33000e-05},
    {1.0e+08, 288.15, {0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.3300000e-05, 9.9997340e-01, 1.3300000e-05}, false, 1.33000e-05}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
