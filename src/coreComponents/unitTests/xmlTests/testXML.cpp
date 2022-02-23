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

  geosx::getGlobalState().getProblemManager().parseInputString( xmlInput );
}

TEST( testXML, testXMLStringExpectedFail )
{
  char const *  xmlInput =
    "<?xml version=\"1.0\" ?>\n"
    "<Problem>\n"
    "  <Solvers>\n"
    "</Problem>";

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
