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

using FlashData = std::tuple<
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< numComps > const,   // total composition
  real64 const              // expected tangent plane distance
  >;

template< EquationOfStateType EOS_TYPE >
class StabilityTestTest9CompFixture :  public ::testing::TestWithParam< FlashData >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr real64 tpdTol = MultiFluidConstants::fugacityTolerance;
public:
  StabilityTestTest9CompFixture()
    : m_fluid( createFluid() )
  {}

  ~StabilityTestTest9CompFixture() = default;

  void testStability( FlashData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EOS_TYPE;
    flashData.vapourEos = EOS_TYPE;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

    real64 const expectedTangentPlaneDistance = std::get< 3 >( data );

    real64 tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
    stackArray1d< real64, numComps > kValues( numComps );

    bool const stabilityStatus = StabilityTest::compute( numComps,
                                                         pressure,
                                                         temperature,
                                                         composition.toSliceConst(),
                                                         componentProperties,
                                                         EOS_TYPE,
                                                         flashData,
                                                         tangentPlaneDistance,
                                                         kValues.toSlice() );

    // Expect this to succeed
    ASSERT_EQ( stabilityStatus, true );

    // Check the tanget plane distance
    checkRelativeError( expectedTangentPlaneDistance, tangentPlaneDistance, relTol, absTol );
  }

  void testFlash( FlashData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EOS_TYPE;
    flashData.vapourEos = EOS_TYPE;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

    real64 tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
    stackArray2d< real64, numComps > kValues( 1, numComps );

    StabilityTest::compute( numComps,
                            pressure,
                            temperature,
                            composition.toSliceConst(),
                            componentProperties,
                            EOS_TYPE,
                            flashData,
                            tangentPlaneDistance,
                            kValues[0] );

    // Now perform the nagative flash
    real64 vapourFraction = -1.0;
    stackArray1d< real64, numComps > liquidComposition( numComps );
    stackArray1d< real64, numComps > vapourComposition( numComps );

    bool const status = NegativeTwoPhaseFlash::compute(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      flashData,
      kValues.toSlice(),
      vapourFraction,
      liquidComposition.toSlice(),
      vapourComposition.toSlice() );

    // Expect this to succeed
    ASSERT_EQ( status, true );

    constexpr real64 epsilon = MultiFluidConstants::epsilon;
    real64 const V = vapourFraction;
    real64 const L = 1.0 - V;

    // Check stability vs flash
    if( tangentPlaneDistance < -tpdTol )
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

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
private:
  static std::unique_ptr< TestFluid< numComps > > createFluid();
};

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
  ::testing::Values(
    FlashData(1.0e+05, 2.8815e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -1.654557e+03),
    FlashData(1.0e+05, 2.9715e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -8.451746e+02),
    FlashData(1.0e+05, 3.5315e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -2.648191e+01),
    FlashData(1.0e+05, 3.9315e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -3.373702e+00),
    FlashData(1.0e+05, 5.7315e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -2.862294e-16)
  )
);

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, SoaveRedlichKwong,
  ::testing::Values(
    FlashData(1.0e+05, 2.8815e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -1.814994e+03),
    FlashData(1.0e+05, 2.9715e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -9.241583e+02),
    FlashData(1.0e+05, 3.5315e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -2.828694e+01),
    FlashData(1.0e+05, 3.9315e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -3.557644e+00),
    FlashData(1.0e+05, 5.7315e+02, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, -2.736526e-16)
  )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
