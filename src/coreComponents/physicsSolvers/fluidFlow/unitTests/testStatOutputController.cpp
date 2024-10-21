/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "unitTests/testingUtilities/TestingTasks.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatistics.hpp"
#include "dataRepository/Group.hpp"
// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

/**
 * @brief this struct is used to provide the input data of each flow tests
 */
struct TestInputs
{
  string xmlInput;
};

/**
 * @brief this struct computes from the test inputs the values to expect from the simulation.
 */
struct TestSet
{
  TestInputs const inputs;
};

namespace ComponentsGeneration
{


TEST( testStatOutputController, checkSinglePhaseFluxStatistics )
{
  string const testSet =
    R"xml(
<Problem>

  <Solvers>
    <SinglePhaseFVM name="testSolver"
                    discretization="singlePhaseTPFA"
                    targetRegions="{ reservoir }" >

      <NonlinearSolverParameters newtonMaxIter="40"
                                 allowNonConverged="1" />
      <LinearSolverParameters solverType="gmres"
                              preconditionerType="iluk"
                              krylovTol="1.0e-6" />

    </SinglePhaseFVM>
  </Solvers>

   <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="singlePhaseTPFA" />
    </FiniteVolume>
  </NumericalMethods>

  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{ C3D8 }"
                  xCoords="{   0, 10 }"
                  yCoords="{   0, 10 }"
                  zCoords="{ -10,  0 }"
                  nx="{ 10 }"
                  ny="{ 10 }"
                  nz="{ 10 }"
                  cellBlockNames="{ cellBlock }" />
  </Mesh>

  <ElementRegions>
    <CellElementRegion name="reservoir"
                       cellBlocks="{ cellBlock }"
                       materialList="{ water, rock }" />
  </ElementRegions>

   <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  defaultDensity="1000"
                                  defaultViscosity="0.001"
                                  referencePressure="0.0"
                                  compressibility="5e-10"
                                  viscosibility="0.0" />

    <CompressibleSolidConstantPermeability name="rock"
                                           solidModelName="nullSolid"
                                           porosityModelName="rockPorosity"
                                           permeabilityModelName="rockPerm" />
    <NullModel name="nullSolid" />
    <PressurePorosity name="rockPorosity"
                      defaultReferencePorosity="0.05"
                      referencePressure="0.0"
                      compressibility="1.0e-9" />
    <ConstantPermeability name="rockPerm"
                          permeabilityComponents="{ 1.0e-12, 1.0e-12, 1.0e-15 }" />
  </Constitutive>

  <FieldSpecifications>
    <!-- The rates are positive, but we need negative values for injection (scale = -1) -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="-1"
                functionName="FluxRate"
                setNames="{ sourceBox }" />
    <!-- sink is producing 3x source rate -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="3"
                functionName="FluxRate"
                setNames="{ sinkBox }"/>

    <HydrostaticEquilibrium name="equil"
                            objectPath="ElementRegions"
                            maxNumberOfEquilibrationIterations="100"
                            datumElevation="-5"
                            datumPressure="1.895e7" />
  </FieldSpecifications>

  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -10.01 }"
         xMax="{  2.01,  1.01,  -8.99 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{  4.99, 8.99, -1.01 }"
         xMax="{ 10.01, 10.01, 0.01 }" />
  </Geometry>

  <Events
    maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />
    <PeriodicEvent
      name="statistics"
      timeFrequency="500"
      target="/Tasks/statController"/>
  </Events>

  <Functions>
    <!-- Unscaled injection / production rate in mol/s -->
    <TableFunction
      name="FluxRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.767,  0.894,  0.561,  0.234,  0.194,  0.178,  0.162,  0.059,  0.000,  0.000 }"
    />
  </Functions>

  <Outputs></Outputs>

  <Tasks>

    <StatOutputController name="statController" >
      <SinglePhaseStatistics
        name="singFlowStatistics"
        flowSolverName="testSolver"
        logLevel="1"
      />
    </StatOutputController>

  </Tasks>
  
  <!-- SPHINX_MESH -->  
  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{ C3D8 }"
                  xCoords="{   0, 10 }"
                  yCoords="{   0, 10 }"
                  zCoords="{ -10,  0 }"
                  nx="{ 10 }"
                  ny="{ 10 }"
                  nz="{ 10 }"
                  cellBlockNames="{ cellBlock }" />
  </Mesh>
  <!-- SPHINX_MESH_END -->  
  
  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -10.01 }"
         xMax="{  2.01,  1.01,  -8.99 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{  4.99, 8.99, -1.01 }"
         xMax="{ 10.01, 10.01, 0.01 }" />
  </Geometry>
  
</Problem>

)xml";


  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  OutputBase::setOutputDirectory( "." );
  setupProblemFromXML( problem, testSet.data() );

  std::vector< string > const packCollectionPaths{
    "/Tasks/packCollectionreservoiraveragePressure",
    "/Tasks/packCollectionreservoirminPressure",
    "/Tasks/packCollectionreservoirmaxPressure",
    "/Tasks/packCollectionreservoirminDeltaPressure",
    "/Tasks/packCollectionreservoirmaxDeltaPressure",
    "/Tasks/packCollectionreservoirtotalMass",
    "/Tasks/packCollectionreservoiraverageTemperature",
    "/Tasks/packCollectionreservoirminTemperature",
    "/Tasks/packCollectionreservoirmaxTemperature",
    "/Tasks/packCollectionreservoirtotalPoreVolume",
    "/Tasks/packCollectionreservoirtotalUncompactedPoreVolume",
    "/Outputs/compFlowHistoryreservoir"
  };

  string const outputPath = " /Outputs/compFlowHistoryreservoir";

  for( string const & path : packCollectionPaths )
  {
    Group const & group = problem.getGroupByPath( path );
    ASSERT_STREQ( path.c_str(), group.getPath().c_str() );
  }


  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

}
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
