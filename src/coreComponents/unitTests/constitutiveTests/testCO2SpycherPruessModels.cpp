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
#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"
//#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp"
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

template< typename MODEL >
std::unique_ptr< MODEL > makePVTFunction( string const & fileContent )
{
  // define component names and molar weight
  string_array componentNames;
  componentNames.resize( 2 );
  componentNames[0] = "co2";
  componentNames[1] = "water";

  array1d< real64 > componentMolarWeight;
  componentMolarWeight.resize( 2 );
  componentMolarWeight[0] = 44e-3;
  componentMolarWeight[1] = 18e-3;

  // read parameters from file
  std::unique_ptr< MODEL > pvtFunction = nullptr;
  array1d< string > const strs = stringutilities::tokenizeBySpaces< array1d >( fileContent );

  pvtFunction = std::make_unique< MODEL >( strs[1],
                                           strs,
                                           componentNames,
                                           componentMolarWeight );

  return pvtFunction;
}

class SpanWagnerCO2DensityTest : public ::testing::Test
{
public:
  SpanWagnerCO2DensityTest()
    : pvtFunction ( makePVTFunction< SpanWagnerCO2Density >( fileContent ))
  {}

  ~SpanWagnerCO2DensityTest() override = default;

protected:
  string const fileContent = "DensityFun SpanWagnerCO2Density 1e6 5e7 1e6 340 360 10";
  std::unique_ptr< SpanWagnerCO2Density > pvtFunction;
};

TEST_F( SpanWagnerCO2DensityTest, spanWagnerCO2DensityMassValuesAndDeriv )
{
  std::cout << pvtFunction.get() << std::endl;
  auto const * tableFunction = FunctionManager::getInstance().getGroupPointer< TableFunction >( "SpanWagnerCO2Density_table" );
  std::cout << tableFunction << std::endl;
  std::cout << tableFunction->getName() << std::endl;
  const auto & coords = tableFunction->getCoordinates();
  const auto & values = tableFunction->getValues();
  integer const np = coords[0].size();
  integer const nt = coords[1].size();
  //integer const nv = values.size();
  std::cout << std::scientific << std::setprecision( 8 );
  for( integer ip=0; ip<np; ip++ )
  {
    std::cout << std::setw( 18 ) << coords[0][ip];
    for( integer it=0; it<nt; it++ )
    {
      std::cout << std::setw( 18 ) << values[it*np+ip];
    }
    std::cout << "\n";
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
