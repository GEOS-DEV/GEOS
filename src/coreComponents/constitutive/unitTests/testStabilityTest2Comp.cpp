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
#include "TestFluid.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr integer numComps = 2;

using StabilityData = std::tuple<
  real64 const,   // pressure
  real64 const,   // temperature
  real64 const,   // CH4 mole fraction
  real64 const    // expected tangent plane distance
  >;

template< EquationOfStateType EOS_TYPE >
class StabilityTestTest2CompFixture :  public ::testing::TestWithParam< StabilityData >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;

public:
  StabilityTestTest2CompFixture()
    : m_fluid( createFluid() )
  {}

  ~StabilityTestTest2CompFixture() = default;

  void testStability( StabilityData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    real64 const zCH4 = std::get< 2 >( data );
    real64 const expectedTangentPlaneDistance = std::get< 3 >( data );

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EOS_TYPE;
    flashData.vapourEos = EOS_TYPE;

    stackArray1d< real64, numComps > composition( numComps );
    composition[0] = zCH4;
    composition[1] = 1.0 - zCH4;

    real64 tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
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
                                                         EOS_TYPE,
                                                         flashData,
                                                         kValues.toSliceConst(),
                                                         flashParameters->m_continuousParameters,
                                                         flashParameters->m_discreteParameters,
                                                         tangentPlaneDistance,
                                                         incipientComposition.toSlice() );

    // Expect this to succeed
    ASSERT_EQ( stabilityStatus, true );

    // Check the tanget plane distance
    checkRelativeError( expectedTangentPlaneDistance, tangentPlaneDistance, relTol, absTol );
  }

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};

private:
  static std::unique_ptr< TestFluid< numComps > > createFluid()
  {
    std::unique_ptr< TestFluid< numComps > > fluid = TestFluid< numComps >::create( {Fluid::CH4, Fluid::C3H8} );
    fluid->setBinaryCoefficients( Feed< 1 >{ 0.1 } );
    return fluid;
  }
};

using PengRobinson = StabilityTestTest2CompFixture< EquationOfStateType::PengRobinson >;
using SoaveRedlichKwong = StabilityTestTest2CompFixture< EquationOfStateType::SoaveRedlichKwong >;
TEST_P( PengRobinson, testStabilityTest )
{
  testStability( GetParam() );
}
TEST_P( SoaveRedlichKwong, testStabilityTest )
{
  testStability( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(StabilityTest, PengRobinson,
  ::testing::ValuesIn<StabilityData>({
    {1.0e+06, 297.15, 0.400, -2.7755576e-16},
    {1.0e+06, 353.15, 0.400, -1.1102230e-16},
    {5.0e+06, 297.15, 0.400, -1.8699196e-01},
    {5.0e+06, 353.15, 0.400,  1.6653345e-16},
    {2.0e+07, 297.15, 0.400, -1.2767565e-15},
    {2.0e+07, 353.15, 0.400, -1.6653345e-16}
  })
);

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, SoaveRedlichKwong,
  ::testing::ValuesIn<StabilityData>({
    {1.0e+06, 297.15, 0.350, -2.7755576e-16},
    {1.0e+06, 353.15, 0.350, -3.8857806e-16},
    {5.0e+06, 297.15, 0.350, -2.0115518e-01},
    {5.0e+06, 353.15, 0.350, -2.7755576e-16},
    {2.0e+07, 297.15, 0.350, -6.6613381e-16},
    {2.0e+07, 353.15, 0.350, -7.7715612e-16}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
