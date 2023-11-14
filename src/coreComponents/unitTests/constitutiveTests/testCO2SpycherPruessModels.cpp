/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
//#include "common/DataTypes.hpp"
//#include "common/TimingMacros.hpp"
//#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineViscosity.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/EzrokhiBrineViscosity.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/FenghourCO2Viscosity.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineDensity.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/BrineEnthalpy.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Enthalpy.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "functions/FunctionManager.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
//using namespace geos::constitutive;
//using namespace geos::dataRepository;
using namespace geos::stringutilities;
using namespace geos::constitutive::PVTProps;
//using namespace geos::constitutive::multifluid;

class CO2SolubilitySpycherPruessTest : public ::testing::Test
{
public:
  CO2SolubilitySpycherPruessTest() = default;
  ~CO2SolubilitySpycherPruessTest() override = default;

protected:
  static std::unique_ptr< CO2Solubility > makeFlashModel( string const & fileContent );

  std::unique_ptr< CO2Solubility > flashModel{nullptr};
};

std::unique_ptr< CO2Solubility >
CO2SolubilitySpycherPruessTest::makeFlashModel( string const & fileContent )
{
  // Define phase names
  string_array phaseNames;
  phaseNames.resize( 2 );
  phaseNames[0] = "gas";
  phaseNames[1] = "liquid";

  // Define component names and molar weight
  string_array componentNames;
  componentNames.resize( 2 );
  componentNames[0] = "co2";
  componentNames[1] = "water";

  array1d< real64 > componentMolarWeight;
  componentMolarWeight.resize( 2 );
  componentMolarWeight[0] = 44.0e-3;
  componentMolarWeight[1] = 18.0e-3;

  // Read file parameters
  array1d< string > const strs = stringutilities::tokenizeBySpaces< array1d >( fileContent );

  return std::make_unique< CO2Solubility >( strs[1],
                                            strs,
                                            phaseNames,
                                            componentNames,
                                            componentMolarWeight );
}

TEST_F( CO2SolubilitySpycherPruessTest, flashCO2SolubilitySpycherPruessTest )
{
  std::cout << std::scientific << std::setprecision( 10 );
  //constexpr char const * fileContent = "FlashModel CO2Solubility 1e5 7.5e7 5e4 285.15 369.15 4.0 0.15";
  constexpr char const * fileContent = "FlashModel CO2Solubility 1e5 2e5 1e5 293.15 303.15 10.0 0.15 1.0e-8 SpycherPruess";
  flashModel = makeFlashModel( fileContent );
  auto const * co2Function = FunctionManager::getInstance().getGroupPointer< TableFunction >( "CO2Solubility_co2Solubility_table" );
  auto const * h2oFunction = FunctionManager::getInstance().getGroupPointer< TableFunction >( "CO2Solubility_waterVaporization_table" );
  const auto & coords = co2Function->getCoordinates();
  const auto & co2Values = co2Function->getValues();
  const auto & h2oValues = h2oFunction->getValues();
  integer const np = coords[0].size();
  integer const nt = coords[1].size();
  //integer const nv = values.size();
  for( integer it=0; it<nt; it++ )
  {
    for( integer ip=0; ip<np; ip++ )
    {
      std::cout << std::setw( 18 ) << coords[1][it];
      std::cout << std::setw( 18 ) << coords[0][ip];
      std::cout << std::setw( 18 ) << h2oValues[it*np+ip];
      std::cout << std::setw( 18 ) << co2Values[it*np+ip];
      std::cout << "\n";
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
