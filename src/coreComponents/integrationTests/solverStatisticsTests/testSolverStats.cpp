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
#include "codingUtilities/Parsing.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"

#include "physicsSolvers/SolverStatistics.hpp"
#include <filesystem>

#include <gtest/gtest.h>

using namespace geos;

bool compareWithTolerance( const std::string & valueStr, double expected, double tolerance )
{
  double value = std::stod( valueStr );
  return std::fabs( value - expected ) <= tolerance;
}

void assertIterationValuesEquals( IterationsStatistics const & stats )
{
  auto const getValue = [&stats]( char const * const key ) -> integer const &
  { return stats.getReference< integer >( key ); };

  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numTimeStepsString() ), 20 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numTimeStepCutsString() ), 0 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numSuccessfulConfigIterationsString() ), 0 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numSuccessfulNonlinearIterationsString() ), 20 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numSuccessfulLinearIterationsString() ), 20 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numDiscardedConfigIterationsString() ), 0 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numDiscardedNonlinearIterationsString() ), 0 );
  EXPECT_EQ( getValue( IterationsStatistics::viewKeyStruct::numDiscardedLinearIterationsString() ), 0 );
}

static const string solverLogOutput =
  R"xml(
    <?xml version="1.0" ?>
    <Problem xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:noNamespaceSchemaLocation="../../src/coreComponents/schema/schema.xsd">
    <Solvers
        gravityVector="{ 0.0, 0.0, -9.81 }">
        <SinglePhaseFVM
        name="SinglePhaseFlow"
        discretization="singlePhaseTPFA"
        targetRegions="{ Channel }"
        writeStatistics="iteration" >
        <NonlinearSolverParameters
            newtonTol="1.0e-6"
            newtonMaxIter="8"/>
        <LinearSolverParameters
            directParallel="0"/>
        </SinglePhaseFVM>
    </Solvers>
    )xml";

static const string solverCSVOutput =
  R"xml(
    <?xml version="1.0" ?>
    <Problem xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:noNamespaceSchemaLocation="../../src/coreComponents/schema/schema.xsd">
        <Solvers
            gravityVector="{ 0.0, 0.0, -9.81 }">
            <SinglePhaseFVM
            name="SinglePhaseFlow"
            discretization="singlePhaseTPFA"
            targetRegions="{ Channel }"
            writeStatistics="convergence" >
            <NonlinearSolverParameters
                newtonTol="1.0e-6"
                newtonMaxIter="8"/>
            <LinearSolverParameters
                directParallel="0"/>
            </SinglePhaseFVM>
        </Solvers>
    )xml";

static const string pattern =
  R"xml(
    <Mesh>
        <InternalMesh
        name="mesh1"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 5, 10 }"
        yCoords="{ 0, 5, 10 }"
        zCoords="{ 0, 2.5, 5, 7.5, 10 }"
        nx="{ 5, 5 }"
        ny="{ 5, 5 }"
        nz="{ 3, 3, 3, 3 }"
        cellBlockNames="{ cb-0_0_0, cb-1_0_0, cb-0_1_0, cb-1_1_0,
                            cb-0_0_1, cb-1_0_1, cb-0_1_1, cb-1_1_1,
                            cb-0_0_2, cb-1_0_2, cb-0_1_2, cb-1_1_2,
                            cb-0_0_3, cb-1_0_3, cb-0_1_3, cb-1_1_3 }"/>
    </Mesh>

    <Geometry>
        <Box
        name="source"
        xMin="{ 7.99, -0.01, 9.99 }"
        xMax="{ 10.01, 2.01, 10.01 }"/>

        
    </Geometry>

    <Events
        maxTime="1e5">

        <PeriodicEvent
        name="solverApplications"
        forceDt="5e3"
        target="/Solvers/SinglePhaseFlow"/>
        
    </Events>

    <ElementRegions>
        <CellElementRegion
        name="Channel"
        cellBlocks="{ cb-1_0_0, cb-0_0_0, cb-0_0_1, cb-0_1_1, cb-0_1_2, cb-1_1_2, cb-1_1_3, cb-1_0_3 }"
        materialList="{ water, rock }"/>

        <CellElementRegion
        name="Barrier"
        cellBlocks="{ cb-0_1_0, cb-1_1_0, cb-1_1_1, cb-1_0_1, cb-1_0_2, cb-0_0_2, cb-0_0_3, cb-0_1_3 }"
        materialList="{ }"/>
    </ElementRegions>

    <NumericalMethods>
        <FiniteVolume>
        <TwoPointFluxApproximation
            name="singlePhaseTPFA"/>
        </FiniteVolume>
    </NumericalMethods>

    <Constitutive>
        <CompressibleSinglePhaseFluid
        name="water"
        defaultDensity="1000"
        defaultViscosity="0.001"
        referencePressure="0.0"
        compressibility="5e-10"
        viscosibility="0.0"/>

        <CompressibleSolidCarmanKozenyPermeability
        name="rock"
        solidModelName="nullSolid"
        porosityModelName="rockPorosity"
        permeabilityModelName="channelPerm"/>

        <NullModel
        name="nullSolid"/>

        <PressurePorosity
        name="rockPorosity"
        defaultReferencePorosity="0.1"
        referencePressure="0.0"
        compressibility="1.0e-9"/>

        <CarmanKozenyPermeability
        name="channelPerm"
        particleDiameter="5e-6"
        sphericity="0.9"/>

    </Constitutive>

    </Problem>
)xml";

CommandLineOptions g_commandLineOptions;

TEST( testSolverStats, testLog )
{
  // setup
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  std::ostringstream xmlInput;
  xmlInput << solverLogOutput << pattern;
  problem.parseInputString( xmlInput.str() );
  problem.problemSetup();
  problem.applyInitialConditions();
  problem.runSimulation();

  PhysicsSolverBase & solver = problem.getGroupByPath< PhysicsSolverBase >( string( "/Solvers/SinglePhaseFlow" ) );
  assertIterationValuesEquals( solver.getIterationStats() );

}

TEST( testSolverStats, testOutputFiles )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  std::ostringstream xmlInput;
  xmlInput << solverCSVOutput << pattern;
  problem.parseInputString( xmlInput.str() );
  problem.problemSetup();
  problem.applyInitialConditions();
  problem.runSimulation();

  PhysicsSolverBase & solver = problem.getGroupByPath< PhysicsSolverBase >( string( "/Solvers/SinglePhaseFlow" ) );
  ConvergenceStatistics const & convergenceStat = solver.getConvergenceStats();
  IterationsStatistics const & iterationStat = solver.getIterationStats();

  auto loadCsvLines = []( string const & filename, stdVector< string > & lines ) {

    if( !std::filesystem::exists( filename ))
    {
      GEOS_ERROR( GEOS_FMT( "Error: File '{}' does not exist!", filename ) );
      GEOS_ERROR( GEOS_FMT( "Current directory: {}", std::filesystem::current_path().string() ) );
      return false;
    }

    std::ifstream is( filename );

    if( !is.is_open())
    {
      GEOS_ERROR( GEOS_FMT( "Error: Cannot open file '{}'", filename ) );

      if( is.fail())
      {
        GEOS_ERROR( "Reason: Opening failed (permissions or file doesn't exist)" );
      }
      if( is.bad())
      {
        GEOS_ERROR( "Reason: Serious stream error" );
      }

      auto perms = std::filesystem::status( filename ).permissions();
      if((perms & std::filesystem::perms::owner_read) == std::filesystem::perms::none )
      {
        GEOS_ERROR( "File does not have read permissions" );
      }

      return false;
    }

    std::string line;
    size_t lineCount = 0;

    while( std::getline( is, line ))
    {
      lines.emplace_back( line );
      lineCount++;
    }

    if( lineCount == 0 )
    {
      GEOS_WARNING( "Warning: File appears to be empty" );
    }
    else
    {
      GEOS_WARNING( "Number of lines read: " );
    }

    is.close();
    return true;
  };

  stdVector< string > csvLines;
  loadCsvLines( iterationStat.getFilename(), csvLines );

  EXPECT_EQ( csvLines[0],
             "Number of time steps,Number of time step cuts,"
             "Successful configuration,Successful nonlinear,Successful linear,"
             "Discarded configuration,Discarded nonlinear,Discarded linear,"
             "Setup time,Solve time" );

  assertIterationValuesEquals( iterationStat );

  stdVector< string > csvLines2;
  loadCsvLines( convergenceStat.getFilename(), csvLines2 );

  EXPECT_EQ( csvLines2[0], "Cycle number,time_n (s),dt (s),iteration,R,Rflow" );

  ASSERT_TRUE( std::remove( convergenceStat.getFilename().c_str() ) == 0 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  basicCleanup();
  return result;
}
