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

  char const * xmlInput =
    R"xml(
    <Problem>
      <Solvers>
        <SolidMechanics_LagrangianFEM
          name="lagsolve"
          cflFactor="0.25"
          discretization="FE1"
          targetRegions="{ Region2 }"/>
      </Solvers>
      <Mesh>
        <InternalMesh
          name="mesh1"
          elementTypes="{ C3D8 }"
          xCoords="{ 0, 3 }"
          yCoords="{ 0, 1 }"
          zCoords="{ 0, 1 }"
          nx="{ 4 }"
          ny="{ 1 }"
          nz="{ 1 }"
          cellBlockNames="{ cb1 }"/>
      </Mesh>
      <Events
        maxTime="1.0e-3">
        <PeriodicEvent
          name="solverApplications"
          forceDt="1.0e-3"
          target="/Solvers/lagsolve"/>
      </Events>
      <NumericalMethods>
        <FiniteElements>
          <FiniteElementSpace
            name="FE1"
            order="1"/>
        </FiniteElements>
      </NumericalMethods>
      <ElementRegions>
        <CellElementRegion
          name="Region2"
          cellBlocks="{ cb1 }"
          materialList="{ shale }"/>
      </ElementRegions>
      <Constitutive>
        <ElasticIsotropic
          name="shale"
          defaultDensity="2700"
          defaultBulkModulus="5.5556e9"
          defaultShearModulus="4.16667e9"/>
      </Constitutive>
    </Problem>
    )xml";

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
                                        "***** Rank 0: Group Mesh (CodeIncludedXML0, l.10) has no child named mesh2\n"
                                        "The children of Mesh are: { mesh1 }";
    // checks if the exception contains the expected message
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" << e.what() <<
                          "  expected sequence :\n" << expectedMsg );
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
