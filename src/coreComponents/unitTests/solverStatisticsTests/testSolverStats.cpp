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

#include <gtest/gtest.h>

using namespace geos;

static const string basicXml =
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
        writeSolverIterationsCSV="1"
        writeSolvingConvergenceCSV="1" >
        <NonlinearSolverParameters
            newtonTol="1.0e-6"
            newtonMaxIter="8"/>
        <LinearSolverParameters
            directParallel="0"/>
        </SinglePhaseFVM>
    </Solvers>

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

TEST( testSolverStats, testOutputFiles )
{
  // setup
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  problem.parseInputString( basicXml );
  problem.problemSetup();
  problem.applyInitialConditions();
  problem.runSimulation();

  PhysicsSolverBase const & solver = problem.getGroupByPath< PhysicsSolverBase >( string( "/Solvers/SinglePhaseFlow" ) );

  SolverStatistics const & solverStat = solver.getSolverStatistics();

  EXPECT_EQ( solverStat.m_iterationsStats.m_numTimeSteps, 20 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numTimeStepCuts, 0 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numSuccessfulOuterLoopIterations, 0 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numSuccessfulNonlinearIterations, 20 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numSuccessfulLinearIterations, 20 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numDiscardedOuterLoopIterations, 0 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numDiscardedNonlinearIterations, 0 );
  EXPECT_EQ( solverStat.m_iterationsStats.m_numDiscardedLinearIterations, 0 );

  auto loadCsvLines = []( string const & filename, std::vector< string > & lines ) {
    std::ifstream is( filename );
    EXPECT_TRUE( is.good());
    string line;
    while( std::getline( is, line ))
    {
      lines.emplace_back( line );
    }
    is.close();
  };

  std::vector< string > csvLines;
  loadCsvLines( "convergence/SinglePhaseFlow_iterations.csv", csvLines );

  EXPECT_EQ( csvLines[0], "Time-steps,Time steps cut,Successful outer loop,Successful nonlinear,Successful linear,Discarded outer loop,Discarded nonlinear,Discarded linear" );
  EXPECT_EQ( csvLines[20], "20,0,0,20,20,0,0,0" );

  std::vector< string > csvLines2;
  loadCsvLines( "convergence/SinglePhaseFlow_convergence.csv", csvLines2 );

  EXPECT_EQ( csvLines2[0], "Time-steps,Newton Iter,RFlow" );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  basicCleanup();
  return result;
}
