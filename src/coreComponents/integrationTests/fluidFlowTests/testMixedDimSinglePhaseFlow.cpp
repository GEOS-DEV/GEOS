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

#include <gtest/gtest.h>
#include <fstream>
#include <tuple>
#include <sstream>
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "codingUtilities/Utilities.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/DomainPartition.hpp"

using namespace geos;

constexpr real64 relative_tolerance = 1.0e-6; // exact up to the order of defaultAperture="1.0e-7"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Integration test for single-phase flow in mixed-dimensional meshes (3D matrix + 2D fractures).
 *
 * This test verifies the correct initialization and solution of pressure fields in a system
 * coupled with fracture elements. It checks both:
 * 1. Initial condition application (when runSolver = false)
 * 2. Steady-state solution with a pressure gradient (when runSolver = true)
 *
 * Parametrized with std::tuple<std::string, bool, std::tuple<int, int, int>>:
 * - std::string: Mesh file name (e.g., "fractured_mesh_hex_DFN_1.vtu").
 * - bool: Whether to run the flow solver (true) or just check initialization (false).
 * - tuple<int, int, int>: Number of partitions in the x, y, z directions.
 */
class MixedDimSinglePhaseFlowTest : public ::testing::TestWithParam< std::tuple< std::string, bool, std::tuple< int, int, int > > >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  std::string generateXmlInput( std::string const & meshFile, std::string const & nodeSetNames, bool const runSolver )
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
    </Geometry>

  <Solvers gravityVector="{0.0, 0.0, 0.0}">
    <SinglePhaseFVM
      name="flowSolver"
      targetRegions="{ Region, Fracture }"
      discretization="tpfa"
      logLevel="1">
      <NonlinearSolverParameters newtonTol="1.0e-6" newtonMaxIter="20"/>
      <LinearSolverParameters directParallel="0"/>    
    </SinglePhaseFVM>
    <SurfaceGenerator name="SurfaceGen" targetRegions="{ Region, Fracture }" fractureRegion="Fracture" initialRockToughness="10.0e9" logLevel="1"/>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="tpfa"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ rockMatrix, water }"/>
    <SurfaceElementRegion name="Fracture" faceBlock="faceElementSubRegion" materialList="{ fractureMaterial, water }" defaultAperture="1.0e-7"/>
  </ElementRegions>

  <Constitutive>
  
    <CompressibleSinglePhaseFluid
      name="water"
      defaultDensity="1.0" 
      defaultViscosity="1.0"
      referenceDensity="1.0" 
      referenceViscosity="1.0"
      referencePressure="0.0" 
      compressibility="0.0" 
      viscosibility="0.0"/>

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

    <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>
    <ConstantPermeability name="fracPerm" permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>      

    <NullModel
      name="nullSolid"/>

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
    <FieldSpecification name="separableFace" fieldName="isFaceSeparable" initialCondition="1" setNames=")xml" << nodeSetNames <<
      R"xml(" objectPath="faceManager" scale="1" />
    <FieldSpecification name="frac" initialCondition="1" setNames=")xml" << nodeSetNames <<
      R"xml(" objectPath="faceManager" fieldName="ruptureState" scale="1" />
    <FieldSpecification name="initialP" fieldName="pressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions" scale="1.5"/>
    <FieldSpecification name="initialPf" fieldName="pressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions/Fracture/faceElementSubRegion" scale="2.0"/>
    <FieldSpecification name="inletP" fieldName="pressure" setNames="{ xnegFace }" objectPath="faceManager" scale="2.0"/>
    <FieldSpecification name="outletP" fieldName="pressure" setNames="{ xposFace }" objectPath="faceManager" scale="1.0"/>
  </FieldSpecifications>
  <Tasks>
    <DynamicFieldSpecification name="apply_fracture_updates" fieldSpecificationNames="{initialPf}"/>
  </Tasks> 
  <Outputs>
    <VTK name="vtkOutputM" fieldNames="{ pressure }" outputRegionType="cell" onlyPlotSpecifiedFieldNames="1"/>
    <VTK name="vtkOutputF" fieldNames="{ pressure }" outputRegionType="surface" onlyPlotSpecifiedFieldNames="1"/>
  </Outputs>
  <Events minTime="-1.0" maxTime="1.0">
    <SoloEvent name="preFracture" target="/Solvers/SurfaceGen" targetTime="-1.0" beginTime="-1.0" />
    <SoloEvent name="TriggerFractureUpdate" target="/Tasks/apply_fracture_updates" targetTime="-1.0" beginTime="-1.0"/>    
    <PeriodicEvent name="outputsM" timeFrequency="1.0" target="/Outputs/vtkOutputM"/>
    <PeriodicEvent name="outputsF" timeFrequency="1.0" target="/Outputs/vtkOutputF"/>
)xml";
    if( runSolver )
    {
      oss << R"xml(  <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver" forceDt="1.0"/>
)xml";
    }
    oss << R"xml(  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( MixedDimSinglePhaseFlowTest, Run )
{
  std::string const & meshFileName = std::get< 0 >( GetParam() );
  bool const runSolver = std::get< 1 >( GetParam() );
  std::tuple< int, int, int > const & partitions = std::get< 2 >( GetParam() );
  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );

  std::string xmlPath = testBinaryDir + "/test_mixed_dim_single_phase_flow.xml";

  std::string nodeSetNames = "{ f1_node_set }";
  if( meshFileName.find( "_DFN_12.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set }";
  }
  else if( meshFileName.find( "_DFN_123.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set, f3_node_set }";
  }

  std::string xmlContent = generateXmlInput( meshFileName, nodeSetNames, runSolver );

  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }

  // Scoped state to ensure cleanup
  {
    auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.push_back( xmlPath );
    options->problemName = "test_mixed_dim_single_phase_flow";

    options->xPartitionsOverride = xPartitions;
    options->yPartitionsOverride = yPartitions;
    options->zPartitionsOverride = zPartitions;
    options->overridePartitionNumbers = true;

    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();

    state.run();

    {
      ProblemManager & pm = state.getProblemManager();
      MeshLevel & mesh = pm.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();

      mesh.getElemManager().forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
      {
        bool const isMatrixCell = dynamic_cast< CellElementSubRegion const * >( &subRegion );
        bool const isFractureCell = dynamic_cast< FaceElementSubRegion const * >( &subRegion );
        if( !isMatrixCell && !isFractureCell )
        {
          return;
        }

        arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
        arrayView2d< real64 const > const center = subRegion.getElementCenter();

        for( localIndex k = 0; k < subRegion.size(); ++k )
        {
          real64 numericalPressure = pressure[k];
          real64 exactPressure = 0.0;
          if( runSolver )
          {
            real64 const x = center[k][0];
            exactPressure = 2.0 * ( 1.0 - x ) + 1.0 * x;
          }
          else
          {
            exactPressure = isMatrixCell ? 1.5 : 2.0;
          }
          real64 const relativeError = std::fabs( numericalPressure - exactPressure ) / exactPressure;
          EXPECT_NEAR( relativeError, 0.0, relative_tolerance ) << "Element " << k << " inexact pressure.";
        }
      } );
    }
  }
}

/**
 * @brief Serial execution test cases (1 MPI rank).
 *
 * Combinations:
 * - Meshes: 1 fracture, 2 fractures, 3 fractures
 * - Solver: Run solver (true) vs Initialization only (false)
 * - Partitioning: 1x1x1 (Serial)
 */
INSTANTIATE_TEST_SUITE_P(
  MixedDimFlowCases,
  MixedDimSinglePhaseFlowTest,
  ::testing::Combine(
    ::testing::Values(
      "fractured_mesh_hex_DFN_1.vtu",
      "fractured_mesh_hex_DFN_12.vtu",
      "fractured_mesh_hex_DFN_123.vtu"
      ),
    ::testing::Bool(),
    ::testing::Values( std::make_tuple( 1, 1, 1 ) )
    )
  );

/**
 * @brief Parallel execution test cases (4 MPI ranks).
 *
 * Combinations:
 * - Meshes: 1 fracture, 2 fractures, 3 fractures
 * - Solver: Run solver (true) vs Initialization only (false)
 * - Partitioning: 1x1x4, 1x2x2, 1x4x1, 2x1x2, 2x2x1, 4x1x1 (4 ranks)
 */
INSTANTIATE_TEST_SUITE_P(
  MixedDimPartitionedFlowCases,
  MixedDimSinglePhaseFlowTest,
  ::testing::Combine(
    ::testing::Values(
      "fractured_mesh_hex_DFN_1.vtu",
      "fractured_mesh_hex_DFN_12.vtu",
      "fractured_mesh_hex_DFN_123.vtu"
      ),
    ::testing::Bool(),
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
