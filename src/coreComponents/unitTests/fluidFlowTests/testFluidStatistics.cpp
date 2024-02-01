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

#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "fieldSpecification/SourceFluxStatistics.hpp"

#include <gtest/gtest.h>


using namespace geos;
using namespace geos::dataRepository;
// using namespace geos::constitutive;
// using namespace geos::constitutive::multifluid;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

struct TestParams
{
  string xmlInput;

  string sourceFluxName;
  string sinkFluxName;

  // input values to test out
  integer sourceElementsCount = 0;
  integer sinkElementsCount = 0;
  integer totalElementsCount = 0;
  std::vector< real64 > timestepSourceMassProd;
  std::vector< real64 > timestepSinkMassProd;
  real64 totalSourceMassProd = 0.0;
  real64 sourceMeanRate = 0.0;
  real64 totalSinkMassProd = 0.0;
  real64 sinkMeanRate = 0.0;
  real64 totalMassProd = 0.0;
  real64 totalMeanRate = 0.0;
};


//////////////////////////////// SinglePhase Flux Statistics Test ////////////////////////////////

TestParams getSinglephaseXmlInput()
{
  TestParams testParams;

  testParams.xmlInput =
    R"xml(
<Problem>

  <Solvers>
    <SinglePhaseFVM name="testSolver"
                    discretization="singlePhaseTPFA"
                    targetRegions="{ reservoir }" >
      <NonlinearSolverParameters newtonTol="1.0e-6"
                                 newtonMaxIter="8" />
      <LinearSolverParameters solverType="gmres"
                              preconditionerType="amg"
                              krylovTol="1.0e-10" />
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
                  xCoords="{ 0, 10 }"
                  yCoords="{ 0, 10 }"
                  zCoords="{ 0, 10 }"
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
    <FieldSpecification name="initialPressure"
                        initialCondition="1"
                        setNames="{ all }"
                        objectPath="ElementRegions/reservoir/cellBlock"
                        fieldName="pressure"
                        scale="5e6" />
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="-1"
                functionName="FluxRate"
                setNames="{ sourceBox }"/>
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="4"
                functionName="FluxRate"
                setNames="{ sinkBox }"/>

  </FieldSpecifications>

  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -0.01 }"
         xMax="{ 2.01, 1.01, 1.01 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{ 5.99, 8.99, -0.01 }"
         xMax="{ 10.01, 10.01, 1.01 }" />
  </Geometry>

  <!-- We are adding 500s to the whole sim time to force the wholeSimStatsEvent to be executed -->
  <Events maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />
    <PeriodicEvent name="timestepsStatsEvent"
                   timeFrequency="500.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/timestepsStats" />
    <PeriodicEvent name="wholeSimStatsEvent"
                   timeFrequency="5000.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/wholeSimStats" />
  </Events>

  <Tasks>
    <SourceFluxStatistics name="timestepsStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver"
                          logLevel="4" />
    <SourceFluxStatistics name="wholeSimStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver"
                          logLevel="4" />
  </Tasks>

  <Functions>
    <TableFunction
      name="FluxRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.767,  0.894,  0.561,  0.234,  0.194,  0.178,  0.162,  0.059,  0.000,  0.000 }"
    />
  </Functions>

</Problem>
)xml";


  testParams.sourceFluxName = "sourceFlux";
  testParams.sinkFluxName = "sinkFlux";

  // compute the expected statistics
  {
    static real64 const dt = 500.0;

    //elements count
    testParams.sourceElementsCount = 2;
    testParams.sinkElementsCount = 4;
    testParams.totalElementsCount = testParams.sourceElementsCount + testParams.sinkElementsCount;

    // FluxRate table from 0.0s to 5000.0s
    std::vector< real64 > const rates = { 0.000, 0.000, 0.767, 0.894, 0.561, 0.234, 0.194, 0.178, 0.162, 0.059, 0.000 };
    testParams.timestepSourceMassProd.reserve( rates.size() );
    testParams.timestepSinkMassProd.reserve( rates.size() );
    for( size_t i = 0; i < rates.size(); ++i )
    {
      // mass injection / production calculation (sink is 4x source production)
      testParams.timestepSourceMassProd.push_back( rates[i] * dt * -1.0 );
      testParams.totalSourceMassProd += testParams.timestepSourceMassProd.back();
      testParams.timestepSinkMassProd.push_back( rates[i] * dt * 4.0 );
      testParams.totalSinkMassProd += testParams.timestepSinkMassProd.back();

      // rates accumulation
      testParams.sourceMeanRate += rates[i] * -1.0;
      testParams.sinkMeanRate += rates[i] * 4.0;
    }
    // mean rates calculation
    real64 const ratesMeanDivisor = 1.0 / double( rates.size() - 1 );
    testParams.sourceMeanRate *= ratesMeanDivisor;
    testParams.sinkMeanRate *= ratesMeanDivisor;
    // totals
    testParams.totalMassProd = testParams.totalSinkMassProd + testParams.totalSourceMassProd;
    testParams.totalMeanRate = testParams.sinkMeanRate + testParams.sourceMeanRate;
  }

  return testParams;
}

TEST( FluidStatisticsTest, checkSinglePhaseFluxStatistics )
{
  TestParams const testParams = getSinglephaseXmlInput();
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  setupProblemFromXML( problem, testParams.xmlInput.data() );

  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  {
    SourceFluxStatsAggregator & wholeSimStats = problem.getGroupByPath< SourceFluxStatsAggregator >( "/Tasks/wholeSimStats" );
    DomainPartition & domain = problem.getDomainPartition();

    // verification that the source flux statistics are correct over the whole simulation
    wholeSimStats.forMeshLevelStatsWrapper( domain,
                                            [&] ( MeshLevel & meshLevel,
                                                  SourceFluxStatsAggregator::WrappedStats & meshLevelStats )
    {
      wholeSimStats.forAllFluxStatsWrappers( meshLevel,
                                             [&] ( MeshLevel &,
                                                   SourceFluxStatsAggregator::WrappedStats & fluxStats )
      {
        if( fluxStats.getFluxName() == testParams.sourceFluxName )
        {
          EXPECT_DOUBLE_EQ( fluxStats.stats().m_producedMass, testParams.totalSourceMassProd ) << "The source flux did not inject the expected total mass.";
          EXPECT_DOUBLE_EQ( fluxStats.stats().m_productionRate, testParams.sourceMeanRate ) << "The source flux did not inject at the expected rate.";
          EXPECT_DOUBLE_EQ( fluxStats.stats().m_elementCount, testParams.sourceElementsCount ) << "The source flux did not target the expected elements.";
        }
        else if( fluxStats.getFluxName() == testParams.sinkFluxName )
        {
          EXPECT_DOUBLE_EQ( fluxStats.stats().m_producedMass, testParams.totalSinkMassProd ) << "The sink flux did not produce the expected total mass.";
          EXPECT_DOUBLE_EQ( fluxStats.stats().m_productionRate, testParams.sinkMeanRate ) << "The sink flux did not produce at the expected rate.";
          EXPECT_DOUBLE_EQ( fluxStats.stats().m_elementCount, testParams.sinkElementsCount ) << "The sink flux did not target the expected elements.";
        }
        else
        {
          FAIL() << "Unexpected SourceFlux found!";
        }
      } );

      EXPECT_DOUBLE_EQ( meshLevelStats.stats().m_producedMass, testParams.totalMassProd ) << "The fluxes did not produce the expected total mass.";
      EXPECT_DOUBLE_EQ( meshLevelStats.stats().m_productionRate, testParams.totalMeanRate ) << "The fluxes did not produce at the expected rate.";
    } );
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
