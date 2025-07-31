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
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/InverseCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"

#include "codingUtilities/UnitTestUtilities.hpp"

#include "dataRepository/xmlWrapper.hpp"
#include "common/initializeEnvironment.hpp"
#include "functions/FunctionManager.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace geos::constitutive;
using namespace geos::dataRepository;

namespace geos
{
namespace testing
{

template< integer NP >
using TestData = std::tuple<
  std::array< real64 const, NP > const,    // capillary pressures
  std::array< real64 const, NP > const,    // phase saturations
  real64 const                             // JFunction multiplier (if required)
  >;

template< typename CAP_PRESSURE, integer NUM_PHASE=2 >
class InverseCapillaryPressureTestFixture : public ::testing::TestWithParam< TestData< NUM_PHASE > >
{
  using CapPressureModel = CAP_PRESSURE;
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;

public:
  InverseCapillaryPressureTestFixture();

  void testInversion( TestData< NUM_PHASE > const & testData );

private:
  conduit::Node m_node;
  std::unique_ptr< Group > m_parent{};
  std::unique_ptr< FunctionManager > m_functionManager{};
  std::unique_ptr< ConstitutiveManager > m_constitutiveManager{};
};

template< typename CAP_PRESSURE, integer NUM_PHASE >
InverseCapillaryPressureTestFixture< CAP_PRESSURE, NUM_PHASE >::InverseCapillaryPressureTestFixture()
  : m_parent( std::make_unique< Group >( "parent", m_node )),
  m_functionManager( std::make_unique< FunctionManager >( FunctionManager::catalogName(), m_parent.get() )),
  m_constitutiveManager( std::make_unique< ConstitutiveManager >( ConstitutiveManager::groupKeyStruct::constitutiveModelsString(), m_parent.get() ))
{
  /* UNCRUSTIFY-OFF */
  string const inputStream =
    "<Constitutive>"
    "<BrooksCoreyCapillaryPressure"
    "    name=\"BrooksCoreyCapillaryPressure3\""
    "    phaseNames=\"{ gas, oil, water }\""
    "    phaseMinVolumeFraction=\"{ 0.05, 0.1, 0.2 }\""
    "    phaseCapPressureExponentInv=\"{ 4, 0, 2 }\""
    "    phaseEntryPressure=\"{ 0.75e6, 2.2e5, 1.2e5 }\""
    "    capPressureEpsilon=\"1e-8\" />"
    "<TableCapillaryPressure"
    "    name=\"TableCapillaryPressure2\""
    "    phaseNames=\"{ water, gas }\""
    "    wettingNonWettingCapPressureTableName=\"PCOW\" />"
    "<TableCapillaryPressure"
    "    name=\"TableCapillaryPressure3\""
    "    phaseNames=\"{ oil, water, gas }\""
    "    wettingIntermediateCapPressureTableName=\"PCOW\""
    "    nonWettingIntermediateCapPressureTableName=\"PCGO\" />"
    "</Constitutive>"
    "<Functions>"
    "<TableFunction"
    "  name=\"PCOW\""
    "  coordinates=\"{ 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90 }\""
    "  values=\"{ 2.08e+05, 1.82e+05, 1.65e+05, 1.53e+05, 1.44e+05, 1.37e+05, 1.31e+05, 1.26e+05, 1.22e+05, 1.18e+05, 1.14e+05, 1.11e+05, 1.09e+05, 1.06e+05, 1.04e+05 }\" />"
    "<TableFunction"
    "  name=\"PCGO\""
    "  coordinates=\"{ 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80 }\""
    "  values=\"{ 0.00e+00, 3.30e+03, 7.19e+03, 1.18e+04, 1.72e+04, 2.35e+04, 3.10e+04, 3.99e+04, 5.03e+04, 6.25e+04, 7.70e+04, 9.40e+04, 1.14e+05, 1.38e+05, 1.65e+05, 1.98e+05, 2.37e+05 }\" />"
    "</Functions>";
  /* UNCRUSTIFY-ON */

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( inputStream );
  if( !xmlResult )
  {
    GEOS_LOG_RANK_0( "XML parsed with errors!" );
    GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.getChild( "Constitutive" );
  m_constitutiveManager->processInputFileRecursive( xmlDocument, xmlConstitutiveNode );
  m_constitutiveManager->postInputInitializationRecursive();

  xmlWrapper::xmlNode xmlFunctionNode = xmlDocument.getChild( "Functions" );
  m_functionManager->processInputFileRecursive( xmlDocument, xmlFunctionNode );
  m_functionManager->postInputInitializationRecursive();
}

template< typename CAP_PRESSURE, integer NUM_PHASE >
void
InverseCapillaryPressureTestFixture< CAP_PRESSURE, NUM_PHASE >::testInversion( TestData< NUM_PHASE > const & testData )
{
  StackArray< real64, 2, NUM_PHASE, compflow::LAYOUT_PHASE > phaseVolumeFraction( 1, NUM_PHASE );
  StackArray< real64, 3, NUM_PHASE, constitutive::cappres::LAYOUT_CAPPRES > capillaryPressure( 1, 1, NUM_PHASE );
  StackArray< real64, 1, NUM_PHASE > jFunctionMultiplier( NUM_PHASE );
  StackArray< real64, 2, NUM_PHASE *NUM_PHASE > dPhaseCapPres_dSaturation( NUM_PHASE, NUM_PHASE );

  auto const & capPres = std::get< 0 >( testData );
  auto const & expectedSaturation = std::get< 1 >( testData );
  auto const & jFunctionScale = std::get< 2 >( testData );
  for( integer ip = 0; ip < NUM_PHASE; ++ip )
  {
    jFunctionMultiplier[ip] = jFunctionScale;
    capillaryPressure[0][0][ip] = capPres[ip];
    phaseVolumeFraction[0][ip] = 1.0/NUM_PHASE;
  }

  string const modelName = GEOS_FMT( "{}{}", CapPressureModel::catalogName(), NUM_PHASE );
  CapPressureModel & model = m_constitutiveManager->getConstitutiveRelation< CapPressureModel >( modelName );

  InverseCapillaryPressure< CapPressureModel > inverse( model );
  auto kernelWrapper = inverse.createKernelWrapper();
  bool status = kernelWrapper.compute( capillaryPressure[0][0].toSliceConst(),
                                       jFunctionMultiplier.toSliceConst(),
                                       phaseVolumeFraction[0] );

  ASSERT_EQ( status, true );

  for( integer ip = 0; ip < NUM_PHASE; ++ip )
  {
    checkRelativeError( expectedSaturation[ip], phaseVolumeFraction[0][ip], relTol, absTol );
  }
}

using TableCapillaryPressure2Inversion = InverseCapillaryPressureTestFixture< TableCapillaryPressure, 2 >;
using TableCapillaryPressure3Inversion = InverseCapillaryPressureTestFixture< TableCapillaryPressure, 3 >;
using BrooksCoreyCapillary3PressureInversion = InverseCapillaryPressureTestFixture< BrooksCoreyCapillaryPressure, 3 >;

TEST_P( TableCapillaryPressure2Inversion, testInversion )
{
  testInversion( GetParam() );
}

TEST_P( TableCapillaryPressure3Inversion, testInversion )
{
  testInversion( GetParam() );
}

TEST_P( BrooksCoreyCapillary3PressureInversion, testInversion )
{
  testInversion( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  CapillaryPressureInversion, TableCapillaryPressure3Inversion,
  ::testing::ValuesIn<TestData< 3 >>({
    {{ 0.00e+00,  2.50e+05,  5.00e+04}, { 0.80000000,  0.20000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  2.50e+05, -0.00e+00}, { 0.80000000,  0.20000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  2.50e+05, -1.00e+03}, { 0.78484848,  0.20000000,  0.01515152}, 1.00e+00},
    {{ 0.00e+00,  2.50e+05, -1.00e+05}, { 0.23500000,  0.20000000,  0.56500000}, 1.00e+00},
    {{ 0.00e+00,  2.50e+05, -2.37e+05}, { 0.00000000,  0.20000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  2.50e+05, -3.00e+05}, { 0.00000000,  0.20000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  2.08e+05,  5.00e+04}, { 0.80000000,  0.20000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  2.08e+05, -0.00e+00}, { 0.80000000,  0.20000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  2.08e+05, -1.00e+03}, { 0.78484848,  0.20000000,  0.01515152}, 1.00e+00},
    {{ 0.00e+00,  2.08e+05, -1.00e+05}, { 0.23500000,  0.20000000,  0.56500000}, 1.00e+00},
    {{ 0.00e+00,  2.08e+05, -2.37e+05}, { 0.00000000,  0.20000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  2.08e+05, -3.00e+05}, { 0.00000000,  0.20000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  2.07e+05,  5.00e+04}, { 0.79807692,  0.20192308,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  2.07e+05, -0.00e+00}, { 0.79807692,  0.20192308,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  2.07e+05, -1.00e+03}, { 0.78292541,  0.20192308,  0.01515152}, 1.00e+00},
    {{ 0.00e+00,  2.07e+05, -1.00e+05}, { 0.23307692,  0.20192308,  0.56500000}, 1.00e+00},
    {{ 0.00e+00,  2.07e+05, -2.37e+05}, {-0.00192308,  0.20192308,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  2.07e+05, -3.00e+05}, {-0.00192308,  0.20192308,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  1.72e+05,  5.00e+04}, { 0.72058824,  0.27941176,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  1.72e+05, -0.00e+00}, { 0.72058824,  0.27941176,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  1.72e+05, -1.00e+03}, { 0.70543672,  0.27941176,  0.01515152}, 1.00e+00},
    {{ 0.00e+00,  1.72e+05, -1.00e+05}, { 0.15558824,  0.27941176,  0.56500000}, 1.00e+00},
    {{ 0.00e+00,  1.72e+05, -2.37e+05}, {-0.07941176,  0.27941176,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  1.72e+05, -3.00e+05}, {-0.07941176,  0.27941176,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  1.04e+05,  5.00e+04}, { 0.00000000,  1.00000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  1.04e+05, -0.00e+00}, { 0.00000000,  1.00000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  1.04e+05, -1.00e+03}, { 0.08484848,  0.90000000,  0.01515152}, 1.00e+00},
    {{ 0.00e+00,  1.04e+05, -1.00e+05}, {-0.46500000,  0.90000000,  0.56500000}, 1.00e+00},
    {{ 0.00e+00,  1.04e+05, -2.37e+05}, {-0.70000000,  0.90000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  1.04e+05, -3.00e+05}, {-0.70000000,  0.90000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  8.60e+04,  5.00e+04}, { 0.00000000,  1.00000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  8.60e+04, -0.00e+00}, { 0.00000000,  1.00000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00,  8.60e+04, -1.00e+03}, { 0.08484848,  0.90000000,  0.01515152}, 1.00e+00},
    {{ 0.00e+00,  8.60e+04, -1.00e+05}, {-0.46500000,  0.90000000,  0.56500000}, 1.00e+00},
    {{ 0.00e+00,  8.60e+04, -2.37e+05}, {-0.70000000,  0.90000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00,  8.60e+04, -3.00e+05}, {-0.70000000,  0.90000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00, -1.20e+05,  5.00e+04}, { 0.00000000,  1.00000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00, -1.20e+05, -0.00e+00}, { 0.00000000,  1.00000000,  0.00000000}, 1.00e+00},
    {{ 0.00e+00, -1.20e+05, -1.00e+03}, { 0.08484848,  0.90000000,  0.01515152}, 1.00e+00},
    {{ 0.00e+00, -1.20e+05, -1.00e+05}, {-0.46500000,  0.90000000,  0.56500000}, 1.00e+00},
    {{ 0.00e+00, -1.20e+05, -2.37e+05}, {-0.70000000,  0.90000000,  0.80000000}, 1.00e+00},
    {{ 0.00e+00, -1.20e+05, -3.00e+05}, {-0.70000000,  0.90000000,  0.80000000}, 1.00e+00}
  })
);

INSTANTIATE_TEST_SUITE_P(
  CapillaryPressureInversion, TableCapillaryPressure2Inversion,
  ::testing::ValuesIn<TestData<2>>({
    {{ 2.25e+05,  0.00e+00}, { 0.20000000,  0.80000000}, 1.00e+00},
    {{ 2.08e+05,  0.00e+00}, { 0.20000000,  0.80000000}, 1.00e+00},
    {{ 2.05e+05,  0.00e+00}, { 0.20576923,  0.79423077}, 1.00e+00},
    {{ 1.86e+05,  0.00e+00}, { 0.24230769,  0.75769231}, 1.00e+00},
    {{ 1.32e+05,  0.00e+00}, { 0.49166667,  0.50833333}, 1.00e+00},
    {{ 1.05e+05,  0.00e+00}, { 0.87500000,  0.12500000}, 1.00e+00},
    {{ 1.04e+05,  0.00e+00}, { 0.90000000,  0.10000000}, 1.00e+00},
    {{ 6.80e+04,  0.00e+00}, { 0.90000000,  0.10000000}, 1.00e+00}
  })
);

INSTANTIATE_TEST_SUITE_P(
  CapillaryPressureInversion, BrooksCoreyCapillary3PressureInversion,
  ::testing::ValuesIn<TestData<3>>({
    {{-7.16e+05,  0.00e+00,  1.09e+05}, { 0.05000000,  0.10000000,  0.85000000}, 1.00e+00},
    {{-7.50e+05,  0.00e+00,  1.09e+05}, { 0.05000000,  0.10000000,  0.85000000}, 1.00e+00},
    {{-7.85e+05,  0.00e+00,  1.09e+05}, { 0.15839859, -0.00839859,  0.85000000}, 1.00e+00},
    {{-7.50e+06,  0.00e+00,  1.09e+05}, { 0.69993500, -0.54993500,  0.85000000}, 1.00e+00},
    {{-7.50e+07,  0.00e+00,  1.09e+05}, { 0.70000000, -0.55000000,  0.85000000}, 1.00e+00},
    {{-7.85e+07,  0.00e+00,  1.09e+05}, { 0.70000000, -0.55000000,  0.85000000}, 1.00e+00},
    {{-7.16e+05,  0.00e+00,  1.20e+05}, { 0.05000000,  0.10000000,  0.85000000}, 1.00e+00},
    {{-7.50e+05,  0.00e+00,  1.20e+05}, { 0.05000000,  0.10000000,  0.85000000}, 1.00e+00},
    {{-7.85e+05,  0.00e+00,  1.20e+05}, { 0.15839859, -0.00839859,  0.85000000}, 1.00e+00},
    {{-7.50e+06,  0.00e+00,  1.20e+05}, { 0.69993500, -0.54993500,  0.85000000}, 1.00e+00},
    {{-7.50e+07,  0.00e+00,  1.20e+05}, { 0.70000000, -0.55000000,  0.85000000}, 1.00e+00},
    {{-7.85e+07,  0.00e+00,  1.20e+05}, { 0.70000000, -0.55000000,  0.85000000}, 1.00e+00},
    {{-7.16e+05,  0.00e+00,  1.32e+05}, { 0.05000000,  0.21280992,  0.73719008}, 1.00e+00},
    {{-7.50e+05,  0.00e+00,  1.32e+05}, { 0.05000000,  0.21280992,  0.73719008}, 1.00e+00},
    {{-7.85e+05,  0.00e+00,  1.32e+05}, { 0.15839859,  0.10441132,  0.73719008}, 1.00e+00},
    {{-7.50e+06,  0.00e+00,  1.32e+05}, { 0.69993500, -0.43712508,  0.73719008}, 1.00e+00},
    {{-7.50e+07,  0.00e+00,  1.32e+05}, { 0.70000000, -0.43719008,  0.73719008}, 1.00e+00},
    {{-7.85e+07,  0.00e+00,  1.32e+05}, { 0.70000000, -0.43719008,  0.73719008}, 1.00e+00},
    {{-7.16e+05,  0.00e+00,  1.20e+07}, { 0.05000000,  0.74999961,  0.20000039}, 1.00e+00},
    {{-7.50e+05,  0.00e+00,  1.20e+07}, { 0.05000000,  0.74999961,  0.20000039}, 1.00e+00},
    {{-7.85e+05,  0.00e+00,  1.20e+07}, { 0.15839859,  0.64153641,  0.20006500}, 1.00e+00},
    {{-7.50e+06,  0.00e+00,  1.20e+07}, { 0.69999904,  0.10000055,  0.20000041}, 1.00e+00},
    {{-7.50e+07,  0.00e+00,  1.20e+07}, { 0.70000000,  0.09999986,  0.20000014}, 1.00e+00},
    {{-7.85e+07,  0.00e+00,  1.20e+07}, { 0.70000000,  0.09999986,  0.20000014}, 1.00e+00},
    {{-7.16e+05,  0.00e+00,  1.20e+09}, { 0.05000000,  0.75000000,  0.20000000}, 1.00e+00},
    {{-7.50e+05,  0.00e+00,  1.20e+09}, { 0.05000000,  0.75000000,  0.20000000}, 1.00e+00},
    {{-7.85e+05,  0.00e+00,  1.20e+09}, { 0.15839859,  0.64160141,  0.20000000}, 1.00e+00},
    {{-7.50e+06,  0.00e+00,  1.20e+09}, { 0.69993500,  0.10006500,  0.20000000}, 1.00e+00},
    {{-7.50e+07,  0.00e+00,  1.20e+09}, { 0.70000000,  0.10000000,  0.20000000}, 1.00e+00},
    {{-7.85e+07,  0.00e+00,  1.20e+09}, { 0.70000000,  0.10000000,  0.20000000}, 1.00e+00},
    {{-7.16e+05,  0.00e+00,  1.32e+09}, { 0.05000000,  0.75000000,  0.20000000}, 1.00e+00},
    {{-7.50e+05,  0.00e+00,  1.32e+09}, { 0.05000000,  0.75000000,  0.20000000}, 1.00e+00},
    {{-7.85e+05,  0.00e+00,  1.32e+09}, { 0.15839859,  0.64160141,  0.20000000}, 1.00e+00},
    {{-7.50e+06,  0.00e+00,  1.32e+09}, { 0.69993500,  0.10006500,  0.20000000}, 1.00e+00},
    {{-7.50e+07,  0.00e+00,  1.32e+09}, { 0.70000000,  0.10000000,  0.20000000}, 1.00e+00},
    {{-7.85e+07,  0.00e+00,  1.32e+09}, { 0.70000000,  0.10000000,  0.20000000}, 1.00e+00}
  })
);

/* UNCRUSTIFY-ON */

} // testing
} // geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::setupEnvironment( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::cleanupEnvironment();

  return result;
}
