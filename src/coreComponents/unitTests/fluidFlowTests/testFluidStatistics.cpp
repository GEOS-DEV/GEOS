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
  std::vector< real64 > sourceRates;
  std::vector< real64 > sinkRates;
  std::vector< real64 > sourceMassProd;
  std::vector< real64 > sinkMassProd;
  real64 sourceMeanRate = 0.0;
  real64 sinkMeanRate = 0.0;
  real64 totalSourceMassProd = 0.0;
  real64 totalSinkMassProd = 0.0;
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
    <!-- The rates are positive, but we need negative values for injection (scale = -1) -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="-1"
                functionName="FluxRate"
                setNames="{ sourceBox }" />
    <!-- sink is producing 4x source rate -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="4"
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
                          flowSolverName="testSolver" />
    <SourceFluxStatistics name="wholeSimStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver" />
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
    testParams.sourceMassProd.reserve( rates.size() );
    testParams.sinkMassProd.reserve( rates.size() );
    for( size_t i = 0; i < rates.size(); ++i )
    {
      // mass injection / production calculation (sink is 4x source production)
      testParams.sourceRates.push_back( rates[i] * -1.0 );
      testParams.sourceMassProd.push_back( rates[i] * dt * -1.0 );
      testParams.totalSourceMassProd += testParams.sourceMassProd.back();
      testParams.sinkRates.push_back( rates[i] * 4.0 );
      testParams.sinkMassProd.push_back( rates[i] * dt * 4.0 );
      testParams.totalSinkMassProd += testParams.sinkMassProd.back();

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

  // run simulation
  setupProblemFromXML( problem, testParams.xmlInput.data() );
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";


  DomainPartition & domain = problem.getDomainPartition();
  // SourceFluxStatsAggregator & timestepsStats = problem.getGroupByPath< SourceFluxStatsAggregator >( "/Tasks/timestepsStats" );
  SourceFluxStatsAggregator & wholeSimStats = problem.getGroupByPath< SourceFluxStatsAggregator >( "/Tasks/wholeSimStats" );

  // check timestep statistics
  {
    //!\\ TODO
  }

  // check whole simulation statistics
  {

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


//////////////////////////////// Multiphase Flux Statistics Test ////////////////////////////////


TestParams getMultiphaseXmlInput()
{
  TestParams testParams;

  testParams.xmlInput =
    R"xml(
<Problem>

  <Solvers>
    <CompositionalMultiphaseFVM name="testSolver"
                                discretization="fluidTPFA"
                                targetRegions="{ reservoir }">
      <NonlinearSolverParameters newtonTol="1.0e-6"
                                 newtonMaxIter="8" />
      <LinearSolverParameters solverType="gmres"
                                 preconditionerType="amg"
                                 krylovTol="1.0e-10" />
    </CompositionalMultiphaseFVM>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="fluidTPFA" />
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
                       materialList="{ fluid, rock, relperm }" />
  </ElementRegions>

  <Constitutive>
    <CO2BrineEzrokhiFluid name="fluid"
                          phaseNames="{ gas, water }"
                          componentNames="{ co2, water }"
                          componentMolarWeight="{ 44e-3, 18e-3 }"
                          phasePVTParaFiles="{ pvtgas.txt, pvtliquidezrokhi.txt }"
                          flashModelParaFile="co2flash.txt" />

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
    <TableRelativePermeability name="relperm"
                               phaseNames="{ gas, water }"
                               wettingNonWettingRelPermTableNames="{ waterRelativePermeabilityTable, gasRelativePermeabilityTable }" />
  </Constitutive>

  <FieldSpecifications>
    <!-- The rates are positive, but we need negative values for injection (scale=-1) -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="-1"
                functionName="FluxRate"
                setNames="{ sourceBox }" />
    <!-- sink is producing 4x source rate -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="4"
                functionName="FluxRate"
                setNames="{ sinkBox }" />
      
    <HydrostaticEquilibrium name="equil"
                            objectPath="ElementRegions"
                            maxNumberOfEquilibrationIterations="100"
                            datumElevation="-5"
                            datumPressure="1.895e7"
                            initialPhaseName="water"
                            componentNames="{ co2, water }"
                            componentFractionVsElevationTableNames="{ initGasCompFracTable, initWaterCompFracTable }"
                            temperatureVsElevationTableName="initTempTable" />
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
                          flowSolverName="testSolver" />
    <SourceFluxStatistics name="wholeSimStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver" />
  </Tasks>

  <Functions>

    <TableFunction
      name="FluxRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.767,  0.894,  0.561,  0.234,  0.194,  0.178,  0.162,  0.059,  0.000,  0.000 }"
    />

    <TableFunction name="initGasCompFracTable"
                   coordinates="{ -10.0, -7.0, -3.0, -1.0 }"
                   values="{        0.0,  0.0,  0.0,  0.0 }" />
    <TableFunction name="initWaterCompFracTable"
                   coordinates="{ -10.0, -7.0, -3.0, -1.0 }"
                   values="{        1.0,  1.0,  1.0,  1.0 }" />
    <TableFunction name="initTempTable"
                   coordinates="{ -10.0,   -7.0,   -3.0,   -1.0 }"
                   values="{  395.15, 389.15, 382.15, 378.15 }" />
    <TableFunction name="waterRelativePermeabilityTable"
                   coordinates="{ 0.3000,          0.3175,         0.3350,         0.3525,        0.3700,        0.3875,       0.4050,       0.4225,       0.4400,       0.4575,       0.4750,      0.4925,      0.5100,      0.5275,      0.5450,      0.5625,      0.5800,      0.5975,     0.6150,     0.6325,     0.6500,     0.6675,     0.6850,     0.7025,     0.7200,     0.7375,     0.7550,     0.7725,    0.7900,    0.8054,    0.8209,    0.8404,    0.8600,    0.8775,    0.8950,    0.9125,    0.9300,    0.9475,    0.9650,    0.9825, 1.0000   }"
                   values="{         0.0, 0.0000001069690, 0.000001523818, 0.000007304599, 0.00002242961, 0.00005398050, 0.0001113999, 0.0002068239, 0.0003554932, 0.0005762517, 0.0008921512, 0.001331180, 0.001927144, 0.002720726, 0.003760776, 0.005105868, 0.006826186, 0.009005830, 0.01174561, 0.01516648, 0.01941368, 0.02466185, 0.03112128, 0.03904542, 0.04874017, 0.06057494, 0.07499593, 0.09254174, 0.1138611, 0.1364565, 0.1632363, 0.2042135, 0.2547712, 0.3097943, 0.3755964, 0.4536528, 0.5451093, 0.6502388, 0.7674166, 0.8909226, 1.000000 }" />
    <TableFunction name="gasRelativePermeabilityTable"
                   coordinateFiles="{ 0.0000,      0.0175,       0.0350,      0.0525,      0.0700,     0.0875,     0.1050,     0.1225,     0.1400,     0.1595,     0.1790,     0.1945,     0.2100,     0.2275,     0.2450,     0.2625,     0.2800,     0.2975,     0.3150,    0.3325,    0.3500,    0.3675,    0.3850,    0.4025,    0.4200,    0.4375,    0.4550,    0.4725,    0.4900,    0.5075,    0.5250,    0.5425,    0.5600,    0.5775,    0.5950,    0.6125,    0.6300,    0.6475,    0.6650,    0.6825,    0.7000 }"
                   values="{     0.000000000, 0.0008885248, 0.002483741, 0.004583224, 0.007135315, 0.01012132, 0.01353719, 0.01738728, 0.02168159, 0.02701850, 0.03295183, 0.03808925, 0.04363513, 0.05042783, 0.05779578, 0.06577020, 0.07438478, 0.08367565, 0.09368138, 0.1044429, 0.1160032, 0.1284076, 0.1417029, 0.1559376, 0.1711607, 0.1874214, 0.2047679, 0.2232459, 0.2428968, 0.2637550, 0.2858446, 0.3091747, 0.3337331, 0.3594782, 0.3863263, 0.4141347, 0.4426735, 0.4715782, 0.5002513, 0.5275887, 0.5500000 }" />
  </Functions>

</Problem>
)xml";


  testParams.sourceFluxName = "sourceFlux";
  testParams.sinkFluxName = "sinkFlux";

  // compute the expected statistics
  {
    //!\\ TODO
  }

  return testParams;
}



//////////////////////////////// Main ////////////////////////////////


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
