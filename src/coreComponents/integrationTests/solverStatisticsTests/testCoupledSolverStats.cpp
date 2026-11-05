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

#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "common/MpiWrapper.hpp"
#include "common/Path.hpp"

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

//////////////////////////////// Utilities ////////////////////////////////

string replaceAll( string str, string_view token, string_view value )
{
  for( size_t pos = 0; ( pos = str.find( token, pos ) ) != string::npos; pos += value.size() )
    str.replace( pos, token.size(), value );
  return str;
}

struct CSVStats
{
  bool exists = false;
  integer headerCount = 0;
  integer dataRowCount = 0;
};

CSVStats inspectCSV( string const & path, string_view headerToken )
{
  CSVStats stats;
  std::ifstream file( path );
  if( !file.is_open() )
    return stats;
  stats.exists = true;
  string line;
  while( std::getline( file, line ) )
  {
    if( line.find( headerToken ) != string::npos )
      stats.headerCount++;
    else if( !line.empty() )
      stats.dataRowCount++;
  }
  return stats;
}

string statsFilePath( string_view solverName, string_view suffix )
{
  return joinPath( OutputBase::getOutputDirectory(),
                   "convergence",
                   GEOS_FMT( "{}_{}.csv", solverName, suffix ) );
}

void expectSingleHeaderCSV( string const & path, string_view headerToken, integer const minDataRows )
{
  CSVStats const stats = inspectCSV( path, headerToken );
  EXPECT_TRUE( stats.exists ) << path;
  EXPECT_EQ( stats.headerCount, 1 ) << path;
  EXPECT_GE( stats.dataRowCount, minDataRows ) << path;
}

//////////////////////////////// Input deck ////////////////////////////////

// Mandel-type consolidation problem (PoroElastic_Mandel_smoke), self-contained:
// initial-displacement table files removed, constant step load instead of loadFunction.
static char const * xmlTemplate =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, 0.0 }">
      <SinglePhasePoromechanics
        name="poroSolve@SUFFIX@"
        solidSolverName="lagsolve@SUFFIX@"
        flowSolverName="flowSolve@SUFFIX@"
        writeStatistics="@CSV_MODE@"
        targetRegions="{ Domain }">
        <NonlinearSolverParameters
          couplingType="@COUPLING_TYPE@"
          lineSearchAction="@LINE_SEARCH_ACTION@"
          newtonMaxIter="10"
          newtonTol="1.0e-3"
          maxTimeStepCuts="2"/>
        <LinearSolverParameters directParallel="0"/>
      </SinglePhasePoromechanics>

      <SolidMechanicsLagrangianFEM
        name="lagsolve@SUFFIX@"
        discretization="FE1"
        writeStatistics="@CSV_MODE@"
        targetRegions="{ Domain }">
        <NonlinearSolverParameters newtonTol="1.0e-4" newtonMaxIter="40"/>
        <LinearSolverParameters directParallel="0"/>
      </SolidMechanicsLagrangianFEM>

      <SinglePhaseFVM
        name="flowSolve@SUFFIX@"
        discretization="singlePhaseTPFA"
        writeStatistics="@CSV_MODE@"
        targetRegions="{ Domain }">
        <NonlinearSolverParameters newtonTol="1.0e-4" newtonMaxIter="40"/>
        <LinearSolverParameters directParallel="0"/>
      </SinglePhaseFVM>
    </Solvers>

    <Mesh>
      <InternalMesh
        name="mesh1" elementTypes="{ C3D8 }"
        xCoords="{ 0.0, 1.0 }" yCoords="{ 0.0, 0.1 }" zCoords="{ 0.0, 1.0 }"
        nx="{ 2 }" ny="{ 1 }" nz="{ 2 }"
        cellBlockNames="{ cb1 }"/>
    </Mesh>

    <NumericalMethods>
      <FiniteElements>
        <FiniteElementSpace name="FE1" order="1"/>
      </FiniteElements>
      <FiniteVolume>
        <TwoPointFluxApproximation name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>

    <ElementRegions>
      <CellElementRegion name="Domain" cellBlocks="{ * }" materialList="{ shale, water }"/>
    </ElementRegions>

    <Constitutive>
      <PorousElasticIsotropic name="shale" solidModelName="shaleSolid"
        porosityModelName="shalePorosity" permeabilityModelName="shalePerm"/>
      <ElasticIsotropic name="shaleSolid" defaultDensity="0"
        defaultBulkModulus="6.6667e7" defaultShearModulus="4.0e7"/>
      <BiotPorosity name="shalePorosity" defaultGrainBulkModulus="1.0e27"
        defaultReferencePorosity="0.375"/>
      <ConstantPermeability name="shalePerm"
        permeabilityComponents="{ 1.0e-12, 0.0, 1.0e-12 }"/>
      <CompressibleSinglePhaseFluid name="water" defaultDensity="1000"
        defaultViscosity="0.001" referencePressure="0.0" referenceDensity="1"
        compressibility="4.4e-10" referenceViscosity="0.001" viscosibility="0.0"/>
    </Constitutive>

    <FieldSpecifications>
      <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }"
        objectPath="ElementRegions/Domain/cb1" fieldName="pressure" scale="0.0"/>
      <FieldSpecification name="xnegconstraint" objectPath="nodeManager"
        fieldName="totalDisplacement" component="0" scale="0.0" setNames="{ xneg }"/>
      <FieldSpecification name="yconstraint" objectPath="nodeManager"
        fieldName="totalDisplacement" component="1" scale="0.0" setNames="{ yneg, ypos }"/>
      <FieldSpecification name="zconstraint" objectPath="nodeManager"
        fieldName="totalDisplacement" component="2" scale="0.0" setNames="{ zneg }"/>
      <FieldSpecification name="normalDisplacement" objectPath="nodeManager"
        fieldName="totalDisplacement" component="2" scale="-1.0e-4" setNames="{ zpos }"/>
      <FieldSpecification name="boundaryPressure" objectPath="faceManager"
        fieldName="pressure" scale="0.0" setNames="{ xpos }"/>
    </FieldSpecifications>

    <Events maxTime="0.05">
      <PeriodicEvent name="solverApplication" forceDt="0.025"
        target="/Solvers/poroSolve@SUFFIX@"/>
    </Events>
  </Problem>
  )xml";

//////////////////////////////// Tests ////////////////////////////////

struct Params
{
  string couplingType;     // "FullyImplicit" or "Sequential"
  string csvMode;          // "all", "none", ...
  string suffix;           // unique per variant so output files do not collide
};

// Run a single coupled-statistics scenario and verify its CSV output.
static void runScenario( Params const & p )
{
  SCOPED_TRACE( GEOS_FMT( "{}_{}", p.couplingType, p.csvMode ) );

  // the serial and parallel variants of this test may run concurrently: keep their output files apart
  string const suffix = GEOS_FMT( "{}_np{}", p.suffix, MpiWrapper::commSize() );
  stdVector< string > const solverNames = { "poroSolve" + suffix, "flowSolve" + suffix, "lagsolve" + suffix };

  // remove files left over by a previous run, so that the checks below only see this run's output
  if( MpiWrapper::commRank() == 0 )
  {
    for( string const & solverName : solverNames )
    {
      std::filesystem::remove( statsFilePath( solverName, "iterations" ) );
      std::filesystem::remove( statsFilePath( solverName, "convergence" ) );
    }
  }
  MpiWrapper::barrier();

  // Sequential coupling forbids line search; use "None" there.
  string const lineSearchAction = ( p.couplingType == "Sequential" ) ? "None" : "Attempt";

  string xmlInput = replaceAll( xmlTemplate, "@COUPLING_TYPE@", p.couplingType );
  xmlInput = replaceAll( xmlInput, "@CSV_MODE@", p.csvMode );
  xmlInput = replaceAll( xmlInput, "@SUFFIX@", suffix );
  xmlInput = replaceAll( xmlInput, "@LINE_SEARCH_ACTION@", lineSearchAction );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  setupProblemFromXML( problem, xmlInput.data() );

  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  MpiWrapper::barrier();
  if( MpiWrapper::commRank() != 0 )
    return;

  string_view const iterationHeaderToken = "Number of time steps";
  string_view const convergenceHeaderToken = "Cycle number";
  string const & coupledName = solverNames[0];

  // 1) Coupled solver CSVs: exist, exactly ONE header (a single rank writes), and data rows.
  expectSingleHeaderCSV( statsFilePath( coupledName, "iterations" ), iterationHeaderToken, 2 ); // 2 forced steps
  // Without subcycling, the Sequential scheme does a single pass and computes no outer-loop residual.
  if( p.couplingType == "FullyImplicit" )
  {
    expectSingleHeaderCSV( statsFilePath( coupledName, "convergence" ), convergenceHeaderToken, 2 );
  }

  // 2) Sub-solver CSVs: in Sequential mode the sub-solvers iterate on their own,
  //    so their iteration and convergence files must contain data.
  if( p.couplingType == "Sequential" )
  {
    for( string const & solverName : { solverNames[1], solverNames[2] } )
    {
      expectSingleHeaderCSV( statsFilePath( solverName, "iterations" ), iterationHeaderToken, 2 );
      expectSingleHeaderCSV( statsFilePath( solverName, "convergence" ), convergenceHeaderToken, 2 );
    }
  }
  // In FullyImplicit mode sub-solvers take no independent time steps,
  // so writeIterationStatsToTable() early-returns: no iteration file expected.
  else
  {
    for( string const & solverName : { solverNames[1], solverNames[2] } )
    {
      EXPECT_FALSE( inspectCSV( statsFilePath( solverName, "iterations" ), iterationHeaderToken ).exists ) << solverName;
    }
  }
}

TEST( CoupledSolverStatisticsCSVTest, checkCSVOutput )
{
  stdVector< Params > const scenarios = {
    { "FullyImplicit", "all", "FIM" },
    { "Sequential", "all", "Seq" },
  };

  for( Params const & p : scenarios )
  {
    runScenario( p );
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
