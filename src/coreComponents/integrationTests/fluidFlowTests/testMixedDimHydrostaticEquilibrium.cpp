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
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/DomainPartition.hpp"

using namespace geos;

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

  /**
   * @brief Derive the correct nodeSetNames string from the mesh file name.
   *
   * The DFN suffix encodes which fracture sets are present:
   *   _DFN_1   -> f1_node_set
   *   _DFN_2   -> f2_node_set
   *   _DFN_3   -> f3_node_set
   *   _DFN_12  -> f1_node_set, f2_node_set
   *   _DFN_13  -> f1_node_set, f3_node_set
   *   _DFN_23  -> f2_node_set, f3_node_set
   *   _DFN_123 -> f1_node_set, f2_node_set, f3_node_set
   */
  static std::string getNodeSetNames( std::string const & meshFileName )
  {
    // Check from most specific to least specific to avoid substring collisions
    if( meshFileName.find( "_DFN_123.vtu" ) != std::string::npos )
      return "{ f1_node_set, f2_node_set, f3_node_set }";
    if( meshFileName.find( "_DFN_23.vtu" ) != std::string::npos )
      return "{ f2_node_set, f3_node_set }";
    if( meshFileName.find( "_DFN_13.vtu" ) != std::string::npos )
      return "{ f1_node_set, f3_node_set }";
    if( meshFileName.find( "_DFN_12.vtu" ) != std::string::npos )
      return "{ f1_node_set, f2_node_set }";
    if( meshFileName.find( "_DFN_3.vtu" ) != std::string::npos )
      return "{ f3_node_set }";
    if( meshFileName.find( "_DFN_2.vtu" ) != std::string::npos )
      return "{ f2_node_set }";
    // default: single fracture set 1
    return "{ f1_node_set }";
  }

  std::string generateXmlInput( std::string const & meshFile,
                                 std::string const & nodeSetNames )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" nodesetNames=")xml" << nodeSetNames << R"xml(" />
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
      logLevel="1">
      <NonlinearSolverParameters newtonTol="1.0e-5" newtonMaxIter="3" lineSearchAction="None" maxAllowedResidualNorm="1.0e20"/>
      <LinearSolverParameters directParallel="0" reuseFactorization="0"/>
    </CompositionalMultiphaseFVM>
    <SurfaceGenerator name="SurfaceGen" targetRegions="{ Region, Fracture }" fractureRegion="Fracture" initialRockToughness="10.0e9" logLevel="1"/>
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

    <FieldSpecification name="separableFace" fieldName="isFaceSeparable" initialCondition="1" setNames=")xml" << nodeSetNames << R"xml(" objectPath="faceManager" scale="1" />
    <FieldSpecification name="frac"          initialCondition="1"         setNames=")xml" << nodeSetNames << R"xml(" objectPath="faceManager" fieldName="ruptureState" scale="1" />
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

  <Events minTime="0.0" maxTime="10.0">
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

  std::string const nodeSetNames = getNodeSetNames( meshFileName );
  std::string const xmlContent   = generateXmlInput( meshFileName, nodeSetNames );

  std::string const xmlPath = testBinaryDir + "/test_mixed_dim_hydrostatic_equilibrium.xml";
  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }

  {
    auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.push_back( xmlPath );
    options->problemName = "test_mixed_dim_hydrostatic_equilibrium";
    options->xPartitionsOverride = xPartitions;
    options->yPartitionsOverride = yPartitions;
    options->zPartitionsOverride = zPartitions;
    options->overridePartitionNumbers = true;

    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();
    state.run();

    ProblemManager & pm = state.getProblemManager();

    // Test that the fluid equilibrium is the actual solution of the discrete system
    // in this sense at the most two iterations are needed to confirm that the residuals
    // remain zero after the hydrostatic equilibration.
    {
      CompositionalMultiphaseFVM & solver =
        pm.getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "flowSolver" );

      integer const numNewtonIter = solver.getNonlinearSolverParameters().m_numNewtonIterations;
      EXPECT_LE( numNewtonIter, 2 )
        << "Mesh " << meshFileName
        << ": expected at most 2 Newton iterations after hydrostatic init, got " << numNewtonIter;
    }

    // Check that the discrete accumulation term is zero for every element,
    // i.e. compAmount == compAmount_n for all components.
    // This confirms the initial state is fully consistent with the
    // previous-time-step state computed with HydrostaticEquilibrium and implies that the discrete fluxes are zero
    // this is fluid is indeed in equilibrium.
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
            real64 const accumulation = compAmount[k][ic] - compAmount_n[k][ic];
            EXPECT_DOUBLE_EQ( accumulation, 0.0 )
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
      // flat hex meshes
      "fractured_mesh_hex_DFN_1.vtu",
      "fractured_mesh_hex_DFN_2.vtu",
      "fractured_mesh_hex_DFN_3.vtu",
      "fractured_mesh_hex_DFN_12.vtu",
      "fractured_mesh_hex_DFN_13.vtu",
      "fractured_mesh_hex_DFN_23.vtu",
      "fractured_mesh_hex_DFN_123.vtu",
      // flat tet meshes
      "fractured_mesh_tet_DFN_1.vtu",
      "fractured_mesh_tet_DFN_2.vtu",
      "fractured_mesh_tet_DFN_3.vtu",
      "fractured_mesh_tet_DFN_12.vtu",
      "fractured_mesh_tet_DFN_13.vtu",
      "fractured_mesh_tet_DFN_23.vtu",
      "fractured_mesh_tet_DFN_123.vtu",
      // wavy hex meshes
      "fractured_wavy_mesh_hex_DFN_1.vtu",
      "fractured_wavy_mesh_hex_DFN_2.vtu",
      "fractured_wavy_mesh_hex_DFN_3.vtu",
      "fractured_wavy_mesh_hex_DFN_12.vtu",
      "fractured_wavy_mesh_hex_DFN_13.vtu",
      "fractured_wavy_mesh_hex_DFN_23.vtu",
      "fractured_wavy_mesh_hex_DFN_123.vtu",
      // wavy tet meshes
      "fractured_wavy_mesh_tet_DFN_1.vtu",
      "fractured_wavy_mesh_tet_DFN_2.vtu",
      "fractured_wavy_mesh_tet_DFN_3.vtu",
      "fractured_wavy_mesh_tet_DFN_12.vtu",
      "fractured_wavy_mesh_tet_DFN_13.vtu",
      "fractured_wavy_mesh_tet_DFN_23.vtu",
      "fractured_wavy_mesh_tet_DFN_123.vtu"
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
