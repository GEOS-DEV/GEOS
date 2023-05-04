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
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

// Tests the Group::getGroup() and getPath() methods
TEST( testGroupPath, testGlobalPaths )
{
  using namespace geos;
  using namespace dataRepository;

  char const *  xmlInput =
    "<?xml version=\"1.0\" ?>\n"
    "<Problem>\n"
    "  <Solvers>\n"
    "    <SolidMechanics_LagrangianFEM\n"
    "      name=\"lagsolve\"\n"
    "      cflFactor=\"0.25\"\n"
    "      discretization=\"FE1\"\n"
    "      targetRegions=\"{ Region2 }\"/>\n"
    "  </Solvers>\n"
    "  <Mesh>\n"
    "    <InternalMesh\n"
    "      name=\"mesh1\"\n"
    "      elementTypes=\"{ C3D8 }\"\n"
    "      xCoords=\"{ 0, 3 }\"\n"
    "      yCoords=\"{ 0, 1 }\"\n"
    "      zCoords=\"{ 0, 1 }\"\n"
    "      nx=\"{ 4 }\"\n"
    "      ny=\"{ 1 }\"\n"
    "      nz=\"{ 1 }\"\n"
    "      cellBlockNames=\"{ cb1 }\"/>\n"
    "  </Mesh>\n"
    "  <Events\n"
    "    maxTime=\"1.0e-3\">\n"
    "    <PeriodicEvent\n"
    "      name=\"solverApplications\"\n"
    "      forceDt=\"1.0e-3\"\n"
    "      target=\"/Solvers/lagsolve\"/>\n"
    "  </Events>\n"
    "  <NumericalMethods>\n"
    "    <FiniteElements>\n"
    "      <FiniteElementSpace\n"
    "        name=\"FE1\"\n"
    "        order=\"1\"/>\n"
    "    </FiniteElements>\n"
    "  </NumericalMethods>\n"
    "  <ElementRegions>\n"
    "    <CellElementRegion\n"
    "      name=\"Region2\"\n"
    "      cellBlocks=\"{ cb1 }\"\n"
    "      materialList=\"{ shale }\"/>\n"
    "  </ElementRegions>\n"
    "  <Constitutive>\n"
    "    <ElasticIsotropic\n"
    "      name=\"shale\"\n"
    "      defaultDensity=\"2700\"\n"
    "      defaultBulkModulus=\"5.5556e9\"\n"
    "      defaultShearModulus=\"4.16667e9\"/>\n"
    "  </Constitutive>\n"
    "</Problem>";

  std::vector< string > const groupPaths{
    "/Mesh/mesh1",
    "/domain/MeshBodies/mesh1/meshLevels/Level0/ElementRegions/elementRegionsGroup/Region2",
    "/domain/Constitutive/shale",
    "/Events/solverApplications",
    "/NumericalMethods/FiniteElements/FE1",
    "/Solvers/lagsolve",
  };

  ProblemManager & problem = getGlobalState().getProblemManager();
  problem.parseInputString( xmlInput );

  for( string const & path : groupPaths )
  {
    Group const & group = problem.getGroupByPath( path );
    ASSERT_STREQ( path.c_str(), group.getPath().c_str() );
  }

  // test for a wrong path given to getGroupByPath()
  bool trowHappened = false;
  try
  {
    problem.getGroupByPath( "/Mesh/mesh2" );
  }
  catch( const std::domain_error & e )
  {
    static constexpr auto expectedMsg = "***** Controlling expression (should be false): child == nullptr\n"
                                        "***** Rank 0: Group /Mesh has no child named mesh2\n"
                                        "The children of Mesh are: { mesh1 }";
    // checks if the exception contains the expected message
    ASSERT_TRUE( string( e.what() ).find( expectedMsg ) != string::npos );
    trowHappened = true;
  }
  // checks if the exception has been thrown as expected
  ASSERT_TRUE( trowHappened );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv, false ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
