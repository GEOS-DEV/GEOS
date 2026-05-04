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
 * @file testMixedDimHydrostaticEquilibriumSplitMeshFixture.hpp
 *
 * Shared test fixture and TEST_P body for compositional multiphase flow in
 * mixed-dimensional meshes initialized with hydrostatic equilibrium,
 * using pre-split (VTM) meshes where fractures are defined via faceBlocks
 * instead of the SurfaceGenerator workflow (WF2).
 *
 * Key XML differences vs WF1 (testMixedDimHydrostaticEquilibriumFixture.hpp):
 *  - VTKMesh uses faceBlocks="{ Fault }" instead of nodesetNames + isFaceSeparable / ruptureState
 *  - SurfaceElementRegion faceBlock="Fault" instead of "faceElementSubRegion"
 *  - No SurfaceGenerator solver, no FieldApplicator fracturing step
 *  - HydrostaticEquilibrium FieldApplicator is still used to apply initial conditions
 *    to the already-existing fracture region
 *  - minTime="0.0" (no pre-fracture step at t=-1)
 *
 * The test verifies the same two properties as WF1:
 * 1. The Newton solver converges in at most 2 (flat) or 4 (wavy) iterations.
 * 2. The discrete accumulation term (compAmount - compAmount_n) is zero everywhere.
 */

#ifndef GEOS_TEST_MIXEDDIMHYDROSTATICEQUILIBRIUMSPLITMESH_FIXTURE_HPP
#define GEOS_TEST_MIXEDDIMHYDROSTATICEQUILIBRIUMSPLITMESH_FIXTURE_HPP

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

constexpr real64 accumulation_tolerance_split = 1.0e-4;

extern CommandLineOptions g_commandLineOptions;

/**
 * @brief Integration test for compositional multiphase flow in pre-split mixed-dimensional
 *        meshes initialized with hydrostatic equilibrium (WF2 - external fracture workflow).
 *
 * Parametrized with std::tuple<std::string, std::tuple<int,int,int>>:
 *  - std::string:           VTM mesh file name (e.g. "fractured_mesh_tet_DFN_1.vtm")
 *  - tuple<int, int, int>:  Number of MPI partitions in x, y, z directions
 */
class MixedDimHydrostaticEquilibriumSplitMeshTest
  : public ::testing::TestWithParam< std::tuple< std::string, std::tuple< int, int, int > > >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  std::string generateXmlInput( std::string const & meshFile )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" faceBlocks="{ Fault }" useGlobalIds="1" checkEulerCharacteristic="0"/>
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
      <NonlinearSolverParameters newtonTol="1.0e-1" newtonMaxIter="4" lineSearchAction="None" maxAllowedResidualNorm="1.0e20"/>
      <LinearSolverParameters directParallel="0" reuseFactorization="0"/>
    </CompositionalMultiphaseFVM>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="tpfa"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ rockMatrix, fluid, RelPerm }"/>
    <SurfaceElementRegion name="Fracture" faceBlock="Fault" materialList="{ fractureMaterial, fluid, RelPerm }" defaultAperture="1.0e-5"/>
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
    <SoloEvent name="triggerFractureUpdate" target="/Tasks/apply_fracture_updates" targetTime="0.0" beginTime="0.0"/>
    <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver" beginTime="0.0" endTime="1.0" forceDt="1.0"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( MixedDimHydrostaticEquilibriumSplitMeshTest, Run )
{
  std::string const & meshFileName = std::get< 0 >( GetParam() );
  std::tuple< int, int, int > const & partitions = std::get< 1 >( GetParam() );
  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );

  std::string const xmlContent = generateXmlInput( meshFileName );

  // Build a unique XML filename per test case to avoid races when multiple
  // test instances run concurrently.
  std::string baseName = meshFileName;
  std::string::size_type const extPos = baseName.rfind( ".vtm" );
  if( extPos != std::string::npos )
    baseName.erase( extPos );
  std::string const xmlName = "test_hydrostatic_split_" + baseName
                              + "_" + std::to_string( xPartitions )
                              + "x" + std::to_string( yPartitions )
                              + "x" + std::to_string( zPartitions );
  std::string const xmlPath = testBinaryDir + "/" + xmlName + ".xml";

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

    // Verify Newton convergence — same thresholds as WF1.
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

    // Verify that the discrete accumulation term is zero for every element.
    {
      MeshLevel & mesh = pm.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();

      mesh.getElemManager().forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
      {
        bool const isMatrixCell   = dynamic_cast< CellElementSubRegion const * >( &subRegion );
        bool const isFractureCell = dynamic_cast< FaceElementSubRegion const * >( &subRegion );
        if( !isMatrixCell && !isFractureCell )
          return;

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
            EXPECT_NEAR( accumulation, 0.0, accumulation_tolerance_split )
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

#endif /* GEOS_TEST_MIXEDDIMHYDROSTATICEQUILIBRIUMSPLITMESH_FIXTURE_HPP */
