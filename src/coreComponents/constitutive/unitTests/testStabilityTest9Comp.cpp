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
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr integer numComps = 9;
static constexpr integer targetComponent = 8;

template< integer NC >
struct StabilityTestData : public TestDataBase< NC >
{
  bool expectedUnstable;         ///< Expected instability status of the mixture
  real64 expectedTargetFraction; ///< Expected mole fraction of targetComponent in the incipient phase
};

using StabilityData = StabilityTestData< numComps >;

template< EquationOfStateType EOS_TYPE >
class StabilityTestTest9CompFixture : public TestFluidFixture< TestFluidName::VolatileOil, StabilityTestData >
{
public:
  ~StabilityTestTest9CompFixture() = default;

  void testStability( StabilityData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EOS_TYPE;
    flashData.vapourEos = EOS_TYPE;

    real64 const pressure = data.pressure;
    real64 const temperature = data.temperature;
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, data.composition );

    bool const expectedUnstableMixture = data.expectedUnstable;
    real64 const expectedZTarget = data.expectedTargetFraction;

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

    // Check the stability labe;
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

    real64 const pressure = data.pressure;
    real64 const temperature = data.temperature;
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, data.composition );

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

    // Expect this to succeed
    ASSERT_EQ( status, true );

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

protected:
private:
};

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
    {{1.0e+05, 288.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99992981},
    {{1.0e+05, 297.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99989343},
    {{1.0e+05, 353.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99906841},
    {{1.0e+05, 393.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99685491},
    {{1.0e+05, 573.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, false, 0.16920000}
  })
);

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, SoaveRedlichKwong,
  ::testing::ValuesIn<StabilityData>({
    {{1.0e+05, 288.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99994255},
    {{1.0e+05, 297.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99991141},
    {{1.0e+05, 353.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99916193},
    {{1.0e+05, 393.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, true,  0.99704625},
    {{1.0e+05, 573.15, {0.0090, 0.0030, 0.5347, 0.1146, 0.0879, 0.0456, 0.0209, 0.0151, 0.1692}}, false, 0.16920000}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
