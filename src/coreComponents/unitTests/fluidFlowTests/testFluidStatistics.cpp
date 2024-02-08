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
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;


//////////////////////////////// Test base utilities ////////////////////////////////


/**
 * @brief this struct is used to provide the input data of each fluid tests
 */
struct TestInputs
{
  string xmlInput;
  std::map< string, string > tableFiles;

  string sourceFluxName;
  string sinkFluxName;
  string timeStepCheckerPath;
  string timeStepFluxStatsPath;
  string wholeSimFluxStatsPath;

  // rates for each timesteps, for each phases
  array2d< real64 > fluxRates;

  real64 dt;
  real64 sourceRateFactor;
  real64 sinkRateFactor;
  integer sourceElementsCount;
  integer sinkElementsCount;

  void setFluxRates( std::initializer_list< std::initializer_list< real64 > > timestepPhaseValues )
  {
    fluxRates.resize( timestepPhaseValues.size(), timestepPhaseValues.begin()->size() );
    integer timestepId = 0;
    for( auto const & phaseValues : timestepPhaseValues )
    {
      integer ip = 0;
      for( auto const & phaseValue : phaseValues )
      {
        fluxRates[timestepId][ip++] = phaseValue;
      }
      ++timestepId;
    }
  }
};

/**
 * @brief this struct computes from the test inputs the values to expect from the simulation.
 */
struct TestSet
{
  TestInputs const inputs;

  integer timestepCount;
  integer totalElementsCount;
  integer phaseCount;

  // stats for each timesteps, for each phases
  array2d< real64 > sourceRates;
  array2d< real64 > sinkRates;
  array2d< real64 > sourceMassProd;
  array2d< real64 > sinkMassProd;
  array2d< real64 > massDeltas;

  // overall simulation stats for each phases
  array1d< real64 > sourceMeanRate;
  array1d< real64 > sinkMeanRate;
  array1d< real64 > totalSourceMassProd;
  array1d< real64 > totalSinkMassProd;
  array1d< real64 > totalMassProd;
  array1d< real64 > totalMeanRate;

  /**
   * @brief Compute the expected statistics set for the tested simulation.
   * @param inputParams the test simulation input parameters
   */
  TestSet( TestInputs const & inputParams ):
    inputs( inputParams )
  {
    timestepCount = inputs.fluxRates.size( 0 );
    phaseCount = inputs.fluxRates.size( 1 );
    totalElementsCount = inputs.sourceElementsCount + inputs.sinkElementsCount;

    sourceRates.resize( timestepCount, phaseCount );
    sinkRates.resize( timestepCount, phaseCount );
    sourceMassProd.resize( timestepCount, phaseCount );
    sinkMassProd.resize( timestepCount, phaseCount );
    massDeltas.resize( timestepCount, phaseCount );
    sourceMeanRate.resize( phaseCount );
    sinkMeanRate.resize( phaseCount );
    totalSourceMassProd.resize( phaseCount );
    totalSinkMassProd.resize( phaseCount );
    totalMassProd.resize( phaseCount );
    totalMeanRate.resize( phaseCount );
    for( integer ip = 0; ip < phaseCount; ++ip )
    {
      for( integer timestepId = 0; timestepId < timestepCount; ++timestepId )
      {
        // mass production / injection calculation
        sourceRates[timestepId][ip] = inputs.fluxRates[timestepId][ip] * inputs.sourceRateFactor;
        sourceMassProd[timestepId][ip] = inputs.fluxRates[timestepId][ip] * inputs.dt * inputs.sourceRateFactor;
        totalSourceMassProd[ip] += sourceMassProd[timestepId][ip];
        sinkRates[timestepId][ip] = inputs.fluxRates[timestepId][ip] * inputs.sinkRateFactor;
        sinkMassProd[timestepId][ip] = inputs.fluxRates[timestepId][ip] * inputs.dt * inputs.sinkRateFactor;
        massDeltas[timestepId][ip] = -( sourceMassProd[timestepId][ip] + sinkMassProd[timestepId][ip] );
        totalSinkMassProd[ip] += sinkMassProd[timestepId][ip];
        // rates accumulations
        sourceMeanRate[ip] += inputs.fluxRates[timestepId][ip] * inputs.sourceRateFactor;
        sinkMeanRate[ip] += inputs.fluxRates[timestepId][ip] * inputs.sinkRateFactor;
      }
      // mean rates calculation
      real64 const ratesMeanDivisor = 1.0 / double( timestepCount - 1 );
      sourceMeanRate[ip] *= ratesMeanDivisor;
      sinkMeanRate[ip] *= ratesMeanDivisor;
      // totals
      totalMassProd[ip] = totalSinkMassProd[ip] + totalSourceMassProd[ip];
      totalMeanRate[ip] = sinkMeanRate[ip] + sourceMeanRate[ip];
    }
  }
};


/**
 * @brief Verification that the source flux statistics are correct for the current timestep
 * @param expectedMasses the expected mass values per phase
 * @param expectedRates the expected rate values per phase
 * @param expectedElementCount the number of expected targeted elements
 * @param stats the timestep stats
 * @param context a context string to provide in any error message.
 */
void checkFluxStats( arraySlice1d< real64 > const & expectedMasses,
                     arraySlice1d< real64 > const & expectedRates,
                     integer const expectedElementCount,
                     SourceFluxStatsAggregator::WrappedStats const & stats,
                     string_view context )
{
  for( int ip = 0; ip < stats.stats().getPhaseCount(); ++ip )
  {
    EXPECT_DOUBLE_EQ( stats.stats().m_producedMass[ip], expectedMasses[ip] ) << GEOS_FMT( "The flux named '{}' did not produce the expected mass ({}, phase = {}).",
                                                                                          stats.getFluxName(), context, ip );
    EXPECT_DOUBLE_EQ( stats.stats().m_productionRate[ip], expectedRates[ip] ) << GEOS_FMT( "The flux named '{}' did not produce at the expected rate ({}, phase = {}).",
                                                                                           stats.getFluxName(), context, ip );
  }
  EXPECT_DOUBLE_EQ( stats.stats().m_elementCount, expectedElementCount ) << GEOS_FMT( "The flux named '{}' did not produce in the expected elements ({}).",
                                                                                      stats.getFluxName(), context );
}

/**
 * @brief Verification that the source flux statistics are correct over the whole simulation
 * @param problem the simulation ProblemManager
 * @param testSet the simulation TestSet
 */
void checkWholeSimFluxStatistics( ProblemManager & problem, TestSet const & testSet )
{
  DomainPartition & domain = problem.getDomainPartition();
  SourceFluxStatsAggregator & wholeSimStats = problem.getGroupByPath< SourceFluxStatsAggregator >( testSet.inputs.wholeSimFluxStatsPath );
  wholeSimStats.forMeshLevelStatsWrapper( domain,
                                          [&] ( MeshLevel & meshLevel,
                                                SourceFluxStatsAggregator::WrappedStats & meshLevelStats )
  {
    wholeSimStats.forAllFluxStatsWrappers( meshLevel,
                                           [&] ( MeshLevel &,
                                                 SourceFluxStatsAggregator::WrappedStats & fluxStats )
    {
      if( fluxStats.getFluxName() == testSet.inputs.sourceFluxName )
      {
        checkFluxStats( testSet.totalSourceMassProd,
                        testSet.sourceMeanRate,
                        testSet.inputs.sourceElementsCount,
                        fluxStats, "over whole simulation" );
      }
      else if( fluxStats.getFluxName() == testSet.inputs.sinkFluxName )
      {
        checkFluxStats( testSet.totalSinkMassProd,
                        testSet.sinkMeanRate,
                        testSet.inputs.sinkElementsCount,
                        fluxStats, "over whole simulation" );
      }
      else
      {
        FAIL() << "Unexpected SourceFlux found!";
      }
    } );

    for( int ip = 0; ip < meshLevelStats.stats().getPhaseCount(); ++ip )
    {
      EXPECT_DOUBLE_EQ( meshLevelStats.stats().m_producedMass[ip], testSet.totalMassProd[ip] ) << "The fluxes did not produce the expected total mass (over whole simulation, phase = " << ip << ").";
      EXPECT_DOUBLE_EQ( meshLevelStats.stats().m_productionRate[ip], testSet.totalMeanRate[ip] ) << "The fluxes did not produce at the expected rate (over whole simulation, phase = " << ip << ").";
    }
  } );
}


/**
 * @brief This Task allows to extract and check each timestep stats during the simulation.
 */
class TimeStepChecker : public TaskBase
{
public:
  TimeStepChecker( string const & name, Group * const parent ):
    TaskBase( name, parent )
  {}

  void postProcessInput() override
  {}

  void setTestSet( TestSet const & testSet ) { m_testSet = &testSet; }
  integer getTestedTimeStepCount() { return m_timestepId; }

  static string catalogName() { return "SinglePhaseStatsTimeStepChecker"; }

  virtual bool execute( real64 const time_n,
                        real64 const GEOS_UNUSED_PARAM( dt ),
                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition & domain )
  {
    EXPECT_LT( m_timestepId, m_testSet->timestepCount ) << "The tested time-step count were higher than expected.";
    SourceFluxStatsAggregator & timestepStats = getGroupByPath< SourceFluxStatsAggregator >( m_testSet->inputs.timeStepFluxStatsPath );
    timestepStats.forMeshLevelStatsWrapper( domain,
                                            [&] ( MeshLevel & meshLevel,
                                                  SourceFluxStatsAggregator::WrappedStats & )
    {
      timestepStats.forAllFluxStatsWrappers( meshLevel,
                                             [&] ( MeshLevel &,
                                                   SourceFluxStatsAggregator::WrappedStats & fluxStats )
      {
        if( fluxStats.getFluxName() == m_testSet->inputs.sourceFluxName )
        {
          checkFluxStats( m_testSet->sourceMassProd[m_timestepId],
                          m_testSet->sourceRates[m_timestepId],
                          m_testSet->inputs.sourceElementsCount,
                          fluxStats, GEOS_FMT( "for timestep at t = {} s", time_n ) );
        }
        else if( fluxStats.getFluxName() == m_testSet->inputs.sinkFluxName )
        {
          checkFluxStats( m_testSet->sinkMassProd[m_timestepId],
                          m_testSet->sinkRates[m_timestepId],
                          m_testSet->inputs.sinkElementsCount,
                          fluxStats, GEOS_FMT( "for timestep at t = {} s", time_n ) );
        }
        else
        {
          FAIL() << "Unexpected SourceFlux found!";
        }
      } );
    } );

    ++m_timestepId;

    return false;
  }
private:
  TestSet const * m_testSet = nullptr;
  int m_timestepId = 0;
};
REGISTER_CATALOG_ENTRY( TaskBase, TimeStepChecker, string const &, Group * const )


class FluidStatisticsTest : public ::testing::Test
{
public:

  void writeTableFiles( std::map< string, string > const & files )
  {
    for( auto const & [fileName, content] : files )
    {
      std::ofstream os( fileName );
      ASSERT_TRUE( os.is_open() );
      os << content;
      os.close();

      m_tableFileNames.push_back( fileName );
    }
  }

  void TearDown() override
  {
    // removing temp table files
    for( string const & fileName : m_tableFileNames )
    {
      ASSERT_TRUE( std::remove( fileName.c_str() ) == 0 );
    }
    m_tableFileNames.clear();
  }

private:
  std::vector< string > m_tableFileNames;
};


//////////////////////////////// SinglePhase Flux Statistics Test ////////////////////////////////
namespace SinglePhaseFluxStatisticsTest
{


TestSet getTestSet()
{
  TestInputs testInputs;

  testInputs.xmlInput =
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
         xMin="{ -0.01, -0.01, -0.01 }"
         xMax="{ 2.01, 1.01, 1.01 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{ 4.99, 8.99, -0.01 }"
         xMax="{ 10.01, 10.01, 1.01 }" />
  </Geometry>

  <!-- We are adding 500s to the whole sim time to force the wholeSimStatsEvent to be executed -->
  <Events maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />
    <PeriodicEvent name="timestepStatsEvent"
                   timeFrequency="500.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/timeStepFluxStats" />
    <PeriodicEvent name="timestepsCheckEvent"
                   timeFrequency="500.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/timeStepChecker" />
    <PeriodicEvent name="wholeSimStatsEvent"
                   timeFrequency="5000.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/wholeSimFluxStats" />
  </Events>

  <Tasks>
    <SourceFluxStatistics name="timeStepFluxStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver" />
    <SourceFluxStatistics name="wholeSimFluxStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver" />

    <SinglePhaseStatsTimeStepChecker name="timeStepChecker" />
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


  testInputs.sourceFluxName = "sourceFlux";
  testInputs.sinkFluxName = "sinkFlux";
  testInputs.timeStepCheckerPath = "/Tasks/timeStepChecker";
  testInputs.timeStepFluxStatsPath = "/Tasks/timeStepFluxStats";
  testInputs.wholeSimFluxStatsPath = "/Tasks/wholeSimFluxStats";

  testInputs.dt = 500.0;
  testInputs.sourceElementsCount = 2;
  testInputs.sinkElementsCount = 5;

  // FluxRate table from 0.0s to 5000.0s
  testInputs.setFluxRates( { { 0.000 },
                             { 0.000 },
                             { 0.767 },
                             { 0.894 },
                             { 0.561 },
                             { 0.234 },
                             { 0.194 },
                             { 0.178 },
                             { 0.162 },
                             { 0.059 },
                             { 0.000 } } );

  // sink is 3x source production
  testInputs.sourceRateFactor = -1.0;
  testInputs.sinkRateFactor = 3.0;

  return TestSet( testInputs );
}

TEST_F( FluidStatisticsTest, checkSinglePhaseFluxStatistics )
{
  TestSet const testSet = getTestSet();

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();

  setupProblemFromXML( problem, testSet.inputs.xmlInput.data() );

  TimeStepChecker & timeStepChecker = problem.getGroupByPath< TimeStepChecker >( testSet.inputs.timeStepCheckerPath );
  timeStepChecker.setTestSet( testSet );

  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  EXPECT_EQ( timeStepChecker.getTestedTimeStepCount(), testSet.timestepCount ) << "The tested time-step were different than expected.";

  checkWholeSimFluxStatistics( problem, testSet );
}


} /* namespace SinglePhaseFluxStatisticsTest */


//////////////////////////////// Multiphase Flux Statistics Test ////////////////////////////////
namespace MultiPhaseFluxStatisticsTest
{


TestSet getTestSet()
{
  TestInputs testInputs;

  testInputs.xmlInput =
    R"xml(
<Problem>

  <Solvers>
    <CompositionalMultiphaseFVM name="testSolver"
                                discretization="fluidTPFA"
                                targetRegions="{ reservoir }"
                                temperature="366.483" >
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
                          phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
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
    <!-- We are injecting CO2 (negative production values) -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="0"
                scale="-1"
                functionName="FluxRate"
                setNames="{ sourceBox }" />
    <!-- We are depleting water -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="1"
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
         xMax="{ 1.01, 1.01, 1.01 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{ 8.99, 8.99, -0.01 }"
         xMax="{ 10.01, 10.01, 1.01 }" />
  </Geometry>

  <!-- We are adding 500s to the whole sim time to force the wholeSimStatsEvent to be executed -->
  <Events maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />
    <PeriodicEvent name="timestepStatsEvent"
                   timeFrequency="500.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/timeStepFluxStats" />
    <PeriodicEvent name="timestepsCheckEvent"
                   timeFrequency="500.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/timeStepChecker" />
    <PeriodicEvent name="wholeSimStatsEvent"
                   timeFrequency="5000.0"
                   targetExactTimestep="1"
                   targetExactStartStop="1"
                   beginTime="0"
                   target="/Tasks/wholeSimFluxStats" />
  </Events>

  <Tasks>
    <SourceFluxStatistics name="timeStepFluxStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver"
                          logLevel="1" />
    <SourceFluxStatistics name="wholeSimFluxStats"
                          fluxNames="{ all }"
                          flowSolverName="testSolver"
                          logLevel="1" />

    <SinglePhaseStatsTimeStepChecker name="timeStepChecker" />
  </Tasks>

  <Functions>
    <!-- CO2 injection rates -->
    <TableFunction
      name="FluxInjectionRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.767,  0.561,  0.194,  0.102,  0.059,  0.000,  0.000,  0.000,  0.000,  0.000 }"
    />
    <!-- water depletion rates -->
    <TableFunction
      name="FluxProductionRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.003,  0.062,  0.121,  0.427,  0.502,  0.199,  0.117,  0.088,  0.059,  0.000 }"
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
                   coordinates="{     0.0000,       0.0175,      0.0350,      0.0525,      0.0700,     0.0875,     0.1050,     0.1225,     0.1400,     0.1595,     0.1790,     0.1945,     0.2100,     0.2275,     0.2450,     0.2625,     0.2800,     0.2975,     0.3150,    0.3325,    0.3500,    0.3675,    0.3850,    0.4025,    0.4200,    0.4375,    0.4550,    0.4725,    0.4900,    0.5075,    0.5250,    0.5425,    0.5600,    0.5775,    0.5950,    0.6125,    0.6300,    0.6475,    0.6650,    0.6825,    0.7000 }"
                   values="{     0.000000000, 0.0008885248, 0.002483741, 0.004583224, 0.007135315, 0.01012132, 0.01353719, 0.01738728, 0.02168159, 0.02701850, 0.03295183, 0.03808925, 0.04363513, 0.05042783, 0.05779578, 0.06577020, 0.07438478, 0.08367565, 0.09368138, 0.1044429, 0.1160032, 0.1284076, 0.1417029, 0.1559376, 0.1711607, 0.1874214, 0.2047679, 0.2232459, 0.2428968, 0.2637550, 0.2858446, 0.3091747, 0.3337331, 0.3594782, 0.3863263, 0.4141347, 0.4426735, 0.4715782, 0.5002513, 0.5275887, 0.5500000 }" />
  </Functions>

</Problem>
)xml";

  testInputs.tableFiles["pvtgas.txt"] = "DensityFun SpanWagnerCO2Density 1.0e5 5.0e7 1e5 285.15 395.15 2\n"
                                        "ViscosityFun FenghourCO2Viscosity 1.0e5 5.0e7 1e5 285.15 395.15 2\n";

  testInputs.tableFiles["pvtliquid.txt"] = "DensityFun EzrokhiBrineDensity 0.1033 -2.2991e-5 -2.3658e-6\n"
                                           "ViscosityFun EzrokhiBrineViscosity 0 0 0\n";

  testInputs.tableFiles["co2flash.txt"] = "FlashModel CO2Solubility 1.0e5 4e7 1e5 285.15 395.15 2 0\n";


  testInputs.sourceFluxName = "sourceFlux";
  testInputs.sinkFluxName = "sinkFlux";
  testInputs.timeStepCheckerPath = "/Tasks/timeStepChecker";
  testInputs.timeStepFluxStatsPath = "/Tasks/timeStepFluxStats";
  testInputs.wholeSimFluxStatsPath = "/Tasks/wholeSimFluxStats";

  testInputs.dt = 500.0;
  testInputs.sourceElementsCount = 1;
  testInputs.sinkElementsCount = 1;
  // FluxRate table from 0.0s to 5000.0s
  testInputs.setFluxRates( {
      { 0.000, 0.000 },
      { 0.000, 0.000 },
      { 0.767, 0.000 },
      { 0.561, 0.000 },
      { 0.194, 0.121 },
      { 0.102, 0.427 },
      { 0.059, 0.502 },
      { 0.000, 0.199 },
      { 0.000, 0.117 },
      { 0.000, 0.088 },
      { 0.000, 0.059 },
      { 0.000, 0.000 } } );

  testInputs.sourceRateFactor = -1.0;
  testInputs.sinkRateFactor = 1.0;

  return TestSet( testInputs );
}


TEST_F( FluidStatisticsTest, checkMultiPhaseFluxStatistics )
{
  TestSet const testSet = getTestSet();
  writeTableFiles( testSet.inputs.tableFiles );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();

  setupProblemFromXML( problem, testSet.inputs.xmlInput.data() );

  //!\\ TODO : récupération du timestepChecker (à ajouter dans le xml)

  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  // EXPECT_EQ( timeStepChecker.getTestedTimeStepCount(), testSet.timestepCount ) << "The tested time-step were different than expected.";

  checkWholeSimFluxStatistics( problem, testSet );
}


}   /* namespace MultiPhaseFluxStatisticsTest */


//////////////////////////////// Main ////////////////////////////////


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
