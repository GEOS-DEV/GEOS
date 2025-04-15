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

#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "fileIO/Outputs/MemoryStatsOutput.hpp"
#include "codingUtilities/Parsing.hpp"
#include "common/format/table/TableFormatter.hpp"

#include <gtest/gtest.h>

using namespace geos;

static const string basicXml =
  R"xml(
<Problem>
  <Mesh>
    <InternalMesh name="mesh1"
                  elementTypes="{ C3D8 }"
                  xCoords="{ 0, 3 }" yCoords="{ 0, 1 }" zCoords="{ 0, 1 }"
                  nx="{ 4 }" ny="{ 1 }" nz="{ 1 }"
                  cellBlockNames="{ cb1 }"/>
  </Mesh>

  <ElementRegions>
    <CellElementRegion name="Region2"
                       cellBlocks="{ * }"
                       materialList="{ shale }"/>
  </ElementRegions>

  <Outputs>
    <MemoryStats name="memoryOutput"
                 writeCSV="1"/>
  </Outputs>

</Problem>
  )xml";
static const string basicSimXml =
  R"xml(
<Problem>
  <!-- <Solvers>
    <SolidMechanicsLagrangianFEM name="lagsolve" 
                                 cflFactor="0.25" 
                                 discretization="FE1" 
                                 targetRegions="{ Region2 }"/>
  </Solvers>-->

  <Mesh>
    <InternalMesh name="mesh1"
                  elementTypes="{ C3D8 }"
                  xCoords="{ 0, 3 }" yCoords="{ 0, 1 }" zCoords="{ 0, 1 }"
                  nx="{ 4 }" ny="{ 1 }" nz="{ 1 }"
                  cellBlockNames="{ cb1 }"/>
  </Mesh>

  <!-- <Events
    maxTime="1.0e-3">
    <PeriodicEvent name="solverApplications"
                   forceDt="1.0e-3"
                   target="/Solvers/lagsolve"/>

    <PeriodicEvent name="outputs"
                   target="/Outputs/siloOutput"/>
  </Events>

  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace name="FE1"
                          order="1"/>
    </FiniteElements>
  </NumericalMethods>-->

  <ElementRegions>
    <CellElementRegion name="Region2"
                       cellBlocks="{ * }"
                       materialList="{ shale }"/>
  </ElementRegions>

  <!-- <Constitutive>
    <ElasticIsotropic name="shale"
                      defaultDensity="2700"
                      defaultBulkModulus="5.5556e9"
                      defaultShearModulus="4.16667e9"/>
  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification name="v0"
                        component="0"
                        fieldName="velocity"
                        functionName="timeFunction"
                        objectPath="nodeManager"
                        scale="1.0"
                        setNames="{ xneg }"/>

    <FieldSpecification name="xconstraint"
                        objectPath="nodeManager"
                        fieldName="velocity"
                        component="0"
                        scale="0.0"
                        setNames="{ xpos }"/>

    <FieldSpecification name="yconstraint"
                        objectPath="nodeManager"
                        fieldName="velocity"
                        component="1"
                        scale="0.0"
                        setNames="{ yneg, ypos }"/>

    <FieldSpecification name="zconstraint"
                        objectPath="nodeManager"
                        fieldName="velocity"
                        component="2"
                        scale="0.0"
                        setNames="{ zneg, zpos }"/>
  </FieldSpecifications>

  <Functions>
    <TableFunction name="timeFunction"
                   inputVarNames="{ time }"
                   coordinates="{ 0.0, 1.0e-6, 2.0e-6, 1.0e9 }"
                   values="{ 0.0, 1.0, 1.0, 1.0 }"/>
  </Functions>-->

  <Outputs>
    <MemoryStats name="memoryOutput"
                 writeCSV="1"/>
  </Outputs>

  <!-- <Geometry>
    <Box name="perf"
         xMin="{ 0, 0, 0 }"
         xMax="{ 1, 1, 1 }"/>
  </Geometry> -->
</Problem>
  )xml";

static const string memOutputPath = "/Outputs/memoryOutput";
static const string memOutputFileName = "MemoryStats_umpireStats.csv";

CommandLineOptions g_commandLineOptions;

// Test if a CSV has been effectively output because of MemoryStatsOutput task execution.
// (This test could be extended to check more memory pools depending on the platform, but may need
// to run a simulation)
TEST( testXML, testMemoryCSVOutput )
{
  // setup
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  problem.parseInputString( basicXml );
  problem.problemSetup();

  // do a MemoryStats test output with a dummy entry
  // EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";
  integer const dummyCycle = 123456;
  string const dummyCycleStr = std::to_string( dummyCycle );
  MemoryStatsOutput & memOutput = problem.getGroupByPath< MemoryStatsOutput >( memOutputPath );
  memOutput.execute( 0.0, 0.0, dummyCycle, 0, 0.0, problem.getDomainPartition() );

  // read the CSV output (parseFile() will throw if no CSV is generated)
  std::vector< string > csvLines;
  {
    std::ifstream is( memOutputFileName );
    EXPECT_TRUE( is.good() );
    string line;
    while( std::getline( is, line ) )
    {
      csvLines.emplace_back( line );
    }
  }

  auto const findCSVMemoryEntry = [=]( string_view cycleStr, string_view memSpaceName ) -> bool {
    for( string const & csvLine : csvLines )
    {
      std::vector< string > const lineEntries = stringutilities::tokenize( csvLine,
                                                                           TableCSVFormatter::m_separator,
                                                                           false );
      EXPECT_GT( lineEntries.size(), 3 );
      if( lineEntries[0] == cycleStr && lineEntries[2] == memSpaceName )
      { // we found the line coresponding to the given cycle & memory-space
        return true;
      }
    }
    return false;
  };
  EXPECT_TRUE( findCSVMemoryEntry( dummyCycleStr, "HOST" ) );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  basicCleanup();
  return result;
}
