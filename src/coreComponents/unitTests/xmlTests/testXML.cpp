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

TEST( testXML, testXMLString )
{
  char const * xmlInput =
    R"xml(
    <?xml version="1.0" ?>
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

  geosx::getGlobalState().getProblemManager().parseInputString( xmlInput );
}

TEST( testXML, testXMLStringExpectedFail )
{
  char const * xmlInput =
    R"xml(
    <?xml version="1.0" ?>
    <Problem>
      <Solvers>
    </Problem>
    )xml";

  EXPECT_THROW( geosx::getGlobalState().getProblemManager().parseInputString( xmlInput ), geosx::InputError );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv, false ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
