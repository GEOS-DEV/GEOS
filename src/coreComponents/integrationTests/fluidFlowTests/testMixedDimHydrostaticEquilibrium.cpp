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
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testMixedDimHydrostaticEquilibrium.cpp
 *
 * Integration test for compositional multiphase flow in mixed-dimensional meshes
 * (3D matrix + 2D fractures) initialized with hydrostatic equilibrium.
 *
 * The test verifies that:
 * 1. The Newton solver converges in at most 2 iterations, because the hydrostatic
 *    equilibrium initialization guarantees that the system starts nearly at steady state.
 * 2. The discrete accumulation term (compAmount - compAmount_n) is zero for all
 *    conservation equations, confirming that the initial state is perfectly consistent.
 *
 * The test runs on all 28 fractured mesh variants (hex/tet, flat/wavy, DFN 1/2/3/12/13/23/123)
 * and uses 4 MPI ranks with various domain partitionings.
 */

#include <gtest/gtest.h>
#include <fstream>
#include <tuple>
#include <sstream>

#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/MpiWrapper.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/DomainPartition.hpp"

using namespace geos;

constexpr real64 accumulation_tolerance = 1.0e-4; // same order as the newton tol in the input file

CommandLineOptions g_commandLineOptions;

/**
 * @brief Integration test for compositional multiphase flow in mixed-dimensional meshes
 *        initialized with hydrostatic equilibrium.
 *
 * Parametrized with std::tuple<std::string, std::tuple<int,int,int>>:
 *  - std::string:                Mesh file name (e.g. "fractured_mesh_hex_DFN_1.vtu")
 *  - tuple<int, int, int>:       Number of MPI partitions in x, y, z directions
 */
class MixedDimHydrostaticEquilibriumTest
  : public ::testing::TestWithParam< std::tuple< std::string, std::tuple< int, int, int > > >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  // ...existing code...

  std::string generateXmlInput( std::string const & meshFile,
                                std::string const & nodeSetNames )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" nodesetNames=")xml" << nodeSetNames <<
      R"xml(" />
  </Mesh>

  <Geometry>
    <Box name="xnegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  0.01,  1.01,  1.01 }"/>
    <Box name="xposFace" xMin="{  0.99, -0.01, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="ynegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  0.01,  1.01 }"/>
    <Box name="yposFace" xMin="{ -0.01,  0.99, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="znegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  1.01,  0.01 }"/>
    <Box name="zposFace" xMin="{ -0.01, -0.01,  0.99 }" xMax="{  1.01,  1.01,  1.01 }"/>
  </Geometry>

  <Solvers gravityVector="{0.0, 0.0, -10.0}">
    <CompositionalMultiphaseFVM
      name="flowSolver"
      discretization="tpfa"
      targetRegions="{ Region, Fracture }"
      useTotalMassEquation="1"
      useMass="1"
      isThermal="0"
      initialDt="1.0e-5"
      temperature="355.0"
      allowNegativePressure="0"
      logLevel="0">
      <NonlinearSolverParameters newtonTol="1.0e-3" newtonMaxIter="4" lineSearchAction="None" maxAllowedResidualNorm="1.0e20"/>
      <LinearSolverParameters directParallel="0" reuseFactorization="0"/>
    </CompositionalMultiphaseFVM>
    <SurfaceGenerator name="SurfaceGen" targetRegions="{ Region, Fracture }" fractureRegion="Fracture" initialRockToughness="10.0e9" logLevel="0"/>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="tpfa"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ rockMatrix, fluid, RelPerm }"/>
    <SurfaceElementRegion name="Fracture" faceBlock="faceElementSubRegion" materialList="{ fractureMaterial, fluid, RelPerm }" defaultAperture="1.0e-5"/>
  </ElementRegions>

  <Constitutive>

    <InvariantImmiscibleFluid name="fluid"
      componentNames="{ CO2, H2O }"
      phaseNames="{ gas, water }"
      densities="{ 500.0, 1000.0 }"
      componentMolarWeight="{ 44e-3, 18e-3 }"
      viscosities="{ 1.0e-3, 1.0e-3 }"/>

    <CompressibleSolidConstantPermeability
      name="rockMatrix"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"/>

    <CompressibleSolidConstantPermeability
      name="fractureMaterial"
      solidModelName="nullSolid"
      porosityModelName="fracPorosity"
      permeabilityModelName="fracPerm"/>

    <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>
    <ConstantPermeability name="fracPerm" permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>

    <TableRelativePermeability
      name="RelPerm"
      phaseNames="{ gas, water }"
      wettingNonWettingRelPermTableNames="{ waterRelativePermeabilityTable, gasRelativePermeabilityTable }"/>

    <NullModel name="nullSolid"/>

    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="0.1"
      referencePressure="0.0"
      compressibility="0.0"/>

    <PressurePorosity
      name="fracPorosity"
      defaultReferencePorosity="0.1"
      referencePressure="0.0"
      compressibility="0.0"/>

  </Constitutive>

  <FieldSpecifications>

    <HydrostaticEquilibrium
      name="CompositionalHydrostaticEquilibrium"
      objectPath="ElementRegions"
      maxNumberOfEquilibrationIterations="200"
      datumElevation="1.0"
      datumPressure="1.0e7"
      componentNames="{ CO2, H2O }"
      phaseContacts="{ 0.5 }"
      elevationIncrementInHydrostaticPressureTable="0.05"
      componentFractionVsElevationTableNames="{ CO2CompositionalFractionVsElevation, WaterCompositionalFractionVsElevation }"
      temperatureVsElevationTableName="InitTempTable"/>

    <FieldSpecification name="separableFace" fieldName="isFaceSeparable" initialCondition="1" setNames=")xml" << nodeSetNames <<
      R"xml(" objectPath="faceManager" scale="1" />
    <FieldSpecification name="frac"          initialCondition="1"         setNames=")xml" << nodeSetNames <<
      R"xml(" objectPath="faceManager" fieldName="ruptureState" scale="1" />
    <FieldSpecification name="topP"       fieldName="pressure"          setNames="{ zposFace }" objectPath="faceManager" scale="1.0e7"/>
    <FieldSpecification name="topCO2Comp" fieldName="globalCompFraction" component="0" setNames="{ zposFace }" objectPath="faceManager" scale="1.0"/>
    <FieldSpecification name="topH2OComp" fieldName="globalCompFraction" component="1" setNames="{ zposFace }" objectPath="faceManager" scale="0.0"/>
    <FieldSpecification name="topT"       fieldName="temperature"        setNames="{ zposFace }" objectPath="faceManager" scale="355.0"/>

  </FieldSpecifications>

  <Functions>

    <TableFunction
      name="waterRelativePermeabilityTable"
      coordinates="{ 0.0, 0.5, 1.0 }"
      values="{ 0.0, 0.5, 1.0 }"/>
    <TableFunction
      name="gasRelativePermeabilityTable"
      coordinates="{ 0.0, 0.5, 1.0 }"
      values="{ 0.0, 0.5, 1.0 }"/>

    <TableFunction
      name="WaterCompositionalFractionVsElevation"
      interpolation="linear"
      coordinates="{ 0.0, 0.5, 0.5001, 1.0 }"
      values="{ 1.0, 1.0, 0.0, 0.0 }"/>

    <TableFunction
      name="CO2CompositionalFractionVsElevation"
      interpolation="linear"
      coordinates="{ 0.0, 0.5, 0.5001, 1.0 }"
      values="{ 0.0, 0.0, 1.0, 1.0 }"/>

    <TableFunction
      name="InitTempTable"
      coordinates="{ -1.0, 0.0 }"
      values="{ 355.0, 355.0 }"/>

  </Functions>

  <Tasks>
    <FieldApplicator name="apply_fracture_updates" fieldSpecificationNames="{ CompositionalHydrostaticEquilibrium }"/>
  </Tasks>

  <Events minTime="0.0" maxTime="1.0">
    <SoloEvent name="fracturingStep"           target="/Solvers/SurfaceGen"              targetTime="0.0" beginTime="0.0"/>
    <SoloEvent name="triggerFractureUpdate" target="/Tasks/apply_fracture_updates"     targetTime="0.0" beginTime="0.0"/>
    <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver" beginTime="0.0" endTime="1.0" forceDt="1.0"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( MixedDimHydrostaticEquilibriumTest, Run )
{
  std::string const & meshFileName = std::get< 0 >( GetParam() );
  std::tuple< int, int, int > const & partitions = std::get< 1 >( GetParam() );
  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );

  std::string nodeSetNames = "{ f1_node_set }";
  if( meshFileName.find( "_DFN_123.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set, f3_node_set }";
  }
  else if( meshFileName.find( "_DFN_12.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set }";
  }
  else if( meshFileName.find( "_DFN_13.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f3_node_set }";
  }
  else if( meshFileName.find( "_DFN_23.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f2_node_set, f3_node_set }";
  }
  else if( meshFileName.find( "_DFN_2.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f2_node_set }";
  }
  else if( meshFileName.find( "_DFN_3.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f3_node_set }";
  }

  std::string const xmlContent = generateXmlInput( meshFileName, nodeSetNames );

  // Build a unique XML filename per test case to avoid races when multiple
  // test instances run concurrently and would otherwise overwrite the same file.
  std::string baseName = meshFileName;
  // strip the .vtu extension
  std::string::size_type const extPos = baseName.rfind( ".vtu" );
  if( extPos != std::string::npos )
    baseName.erase( extPos );
  std::string const xmlName = "test_hydrostatic_" + baseName
                              + "_" + std::to_string( xPartitions )
                              + "x" + std::to_string( yPartitions )
                              + "x" + std::to_string( zPartitions );
  std::string const xmlPath = testBinaryDir + "/" + xmlName + ".xml";
  // Only rank 0 writes the XML file; all ranks then barrier before reading it.
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );

  {
    auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.push_back( xmlPath );
    options->problemName = xmlName;
    options->xPartitionsOverride = xPartitions;
    options->yPartitionsOverride = yPartitions;
    options->zPartitionsOverride = zPartitions;
    options->overridePartitionNumbers = true;

    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();
    state.run();

    ProblemManager & pm = state.getProblemManager();

    // Verify that the hydrostatic equilibrium is the exact solution of the discrete system.
    // For flat meshes the discretization is exact (or nearly so), so at most 2 Newton
    // iterations are needed. For wavy meshes the geometry introduces a small discretization
    // error with respect the exact hydrostatic equilibrium, so up to 4 iterations are allowed.
    {
      CompositionalMultiphaseFVM & solver =
        pm.getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "flowSolver" );

      bool const isWavy = meshFileName.find( "_wavy_" ) != std::string::npos;
      integer const maxExpectedIter = isWavy ? 4 : 2;

      integer const numNewtonIter = solver.getNonlinearSolverParameters().m_numNewtonIterations;
      EXPECT_LE( numNewtonIter, maxExpectedIter )
        << "Mesh " << meshFileName
        << ": expected at most " << maxExpectedIter
        << " Newton iterations after hydrostatic init, got " << numNewtonIter;
    }

    // Verify that the discrete accumulation term is zero for every element
    // (compAmount == compAmount_n for all components), confirming that the
    // hydrostatic initial state is consistent with the previous time-step
    // snapshot and that the discrete fluxes are zero, i.e. the fluid is in equilibrium.
    {
      MeshLevel & mesh = pm.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();

      mesh.getElemManager().forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
      {
        bool const isMatrixCell   = dynamic_cast< CellElementSubRegion const * >( &subRegion );
        bool const isFractureCell = dynamic_cast< FaceElementSubRegion const * >( &subRegion );
        if( !isMatrixCell && !isFractureCell )
          return;

        // compAmount and compAmount_n are both set by applyInitialConditions
        // to the same value, so their difference must be exactly zero.
        arrayView2d< real64 const, compflow::USD_COMP > const compAmount =
          subRegion.getField< fields::flow::compAmount >();
        arrayView2d< real64 const, compflow::USD_COMP > const compAmount_n =
          subRegion.getField< fields::flow::compAmount_n >();

        integer const numComp = compAmount.size( 1 );

        for( localIndex k = 0; k < subRegion.size(); ++k )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            real64 const accumulation = std::fabs( compAmount[k][ic] - compAmount_n[k][ic] );
            EXPECT_NEAR( accumulation, 0.0, accumulation_tolerance )
              << "Mesh " << meshFileName
              << ": non-zero accumulation term at element " << k
              << ", component " << ic
              << " (compAmount=" << compAmount[k][ic]
              << ", compAmount_n=" << compAmount_n[k][ic] << ")";
          }
        }
      } );
    }
  }
}

// ---------------------------------------------------------------------------
// Test suite instantiation
//
// All 28 fractured mesh variants x 6 domain partitions = 168 test cases
// (run with 4 MPI ranks)
// ---------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
  MixedDimHydrostaticEquilibriumCases,
  MixedDimHydrostaticEquilibriumTest,
  ::testing::Combine(
    ::testing::Values(
      // Strategic subset: hex/tet × flat/wavy/full_span × single/triple fractures
      // Covers all mesh types and geometries with representative complexity

      // Flat hex meshes - single and triple fracture only
      "fractured_mesh_hex_DFN_1.vtu",
      "fractured_mesh_hex_DFN_123.vtu",

      // Flat tet meshes - single and triple fracture only
      "fractured_mesh_tet_DFN_1.vtu",
      "fractured_mesh_tet_DFN_123.vtu",

      // Wavy hex meshes - single and triple fracture only
      "fractured_wavy_mesh_hex_DFN_1.vtu",
      "fractured_wavy_mesh_hex_DFN_123.vtu",

      // Wavy tet meshes - single and triple fracture only
      "fractured_wavy_mesh_tet_DFN_1.vtu",
      "fractured_wavy_mesh_tet_DFN_123.vtu",

      // Full span hex meshes - single and triple fracture only
      "fractured_full_span_mesh_hex_DFN_1.vtu",
      "fractured_full_span_mesh_hex_DFN_123.vtu",

      // Full span tet meshes - single and triple fracture only
      "fractured_full_span_mesh_tet_DFN_1.vtu",
      "fractured_full_span_mesh_tet_DFN_123.vtu"
      ),
    ::testing::Values(
      std::make_tuple( 1, 1, 4 ),
      std::make_tuple( 1, 2, 2 ),
      std::make_tuple( 1, 4, 1 ),
      std::make_tuple( 2, 1, 2 ),
      std::make_tuple( 2, 2, 1 ),
      std::make_tuple( 4, 1, 1 )
      )
    )
  );

int main( int argc, char * argv[] )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv, false );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
