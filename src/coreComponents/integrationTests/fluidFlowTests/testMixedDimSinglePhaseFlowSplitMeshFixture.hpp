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
 * @file testMixedDimSinglePhaseFlowSplitMeshFixture.hpp
 *
 * Shared test fixture and TEST_P body for single-phase flow in mixed-dimensional
 * pre-split meshes (VTM format: 3D matrix + 2D fractures defined via faceBlocks).
 * Both the serial and MPI test executables include this header and only differ in
 * the INSTANTIATE_TEST_SUITE_P parameters.
 *
 * Key differences from the WF1 (SurfaceGenerator) fixture:
 *  - VTKMesh uses faceBlocks="{ Fault }" instead of nodesetNames
 *  - No SurfaceGenerator solver, no preFracture event
 *  - SurfaceElementRegion uses faceBlock="Fault"
 *  - Fracture pressure initialized directly via objectPath="ElementRegions/Fracture/Fault"
 */

#ifndef GEOS_TEST_MIXEDDIMSINGLEPHASEFLOW_SPLITMESH_FIXTURE_HPP
#define GEOS_TEST_MIXEDDIMSINGLEPHASEFLOW_SPLITMESH_FIXTURE_HPP

#include <gtest/gtest.h>
#include <fstream>
#include <tuple>
#include <sstream>
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/MpiWrapper.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/DomainPartition.hpp"

using namespace geos;

constexpr real64 relative_tolerance_split = 1.0e-6;

extern CommandLineOptions g_commandLineOptions;

/**
 * @brief Integration test for single-phase flow in pre-split mixed-dimensional meshes.
 *
 * Mirrors MixedDimSinglePhaseFlowTest but uses VTM composite meshes where fractures
 * are defined via faceBlocks instead of the SurfaceGenerator workflow.
 *
 * Parametrized with std::tuple<std::string, bool, std::tuple<int, int, int>>:
 * - std::string: VTM mesh file name (e.g., "fractured_mesh_hex_DFN_123.vtm").
 * - bool: Whether to run the flow solver (true) or just check initialization (false).
 * - tuple<int, int, int>: Number of partitions in the x, y, z directions.
 */
class MixedDimSinglePhaseFlowSplitMeshTest
  : public ::testing::TestWithParam< std::tuple< std::string, bool, std::tuple< int, int, int > > >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  std::string generateXmlInput( std::string const & meshFile, bool const runSolver )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" useGlobalIds="1" file=")xml" << meshFile << R"xml(" faceBlocks="{ Fault }"/>
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
      logLevel="0">
      <NonlinearSolverParameters newtonTol="1.0e-6" newtonMaxIter="20"/>
      <LinearSolverParameters directParallel="0"/>
    </SinglePhaseFVM>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="tpfa"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ rockMatrix, water }"/>
    <SurfaceElementRegion name="Fracture" faceBlock="Fault" materialList="{ fractureMaterial, water }" defaultAperture="1.0e-7"/>
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
    <FieldSpecification name="initialP"  fieldName="pressure" initialCondition="1"
                        setNames="{ all }" objectPath="ElementRegions" scale="1.5"/>
    <FieldSpecification name="initialPf" fieldName="pressure" initialCondition="1"
                        setNames="{ all }" objectPath="ElementRegions/Fracture/Fault" scale="2.0"/>
    <FieldSpecification name="inletP"  fieldName="pressure" setNames="{ xnegFace }"
                        objectPath="faceManager" scale="2.0"/>
    <FieldSpecification name="outletP" fieldName="pressure" setNames="{ xposFace }"
                        objectPath="faceManager" scale="1.0"/>
  </FieldSpecifications>

  <Events minTime="0.0" maxTime="1.0">
)xml";
    if( runSolver )
    {
      oss << R"xml(    <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver" forceDt="1.0"/>
)xml";
    }
    oss << R"xml(  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( MixedDimSinglePhaseFlowSplitMeshTest, Run )
{
  std::string const & meshFileName = std::get< 0 >( GetParam() );
  bool const runSolver = std::get< 1 >( GetParam() );
  std::tuple< int, int, int > const & partitions = std::get< 2 >( GetParam() );
  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );

  std::string xmlPath = testBinaryDir + "/test_mixed_dim_single_phase_flow_split_mesh.xml";

  std::string xmlContent = generateXmlInput( meshFileName, runSolver );

  // Only rank 0 writes the XML file; barrier ensures all ranks see it before proceeding.
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );

  // Scoped state to ensure cleanup
  {
    auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.push_back( xmlPath );
    options->problemName = "test_mixed_dim_single_phase_flow_split_mesh";

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
        bool const isMatrixCell   = dynamic_cast< CellElementSubRegion const * >( &subRegion ) != nullptr;
        bool const isFractureCell = dynamic_cast< FaceElementSubRegion const * >( &subRegion ) != nullptr;
        if( !isMatrixCell && !isFractureCell )
        {
          return;
        }

        arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
        arrayView2d< real64 const > const center   = subRegion.getElementCenter();

        localIndex const numElems = subRegion.size();

        RAJA::ReduceMax< parallelDeviceReduce, real64 > maxRelError( 0.0 );

        if( runSolver )
        {
          forAll< geos::parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            real64 const x = center[k][0];
            real64 const exactPressure = 2.0 * ( 1.0 - x ) + 1.0 * x;
            real64 const relErr = LvArray::math::abs( pressure[k] - exactPressure ) / exactPressure;
            maxRelError.max( relErr );
          } );
        }
        else
        {
          real64 const exactPressure = isMatrixCell ? 1.5 : 2.0;
          forAll< geos::parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            real64 const relErr = LvArray::math::abs( pressure[k] - exactPressure ) / exactPressure;
            maxRelError.max( relErr );
          } );
        }

        EXPECT_NEAR( maxRelError.get(), 0.0, relative_tolerance_split )
          << ( isMatrixCell ? "Matrix" : "Fracture" ) << " subregion has inexact pressure.";
      } );
    }
  }
}

#endif /* GEOS_TEST_MIXEDDIMSINGLEPHASEFLOW_SPLITMESH_FIXTURE_HPP */
