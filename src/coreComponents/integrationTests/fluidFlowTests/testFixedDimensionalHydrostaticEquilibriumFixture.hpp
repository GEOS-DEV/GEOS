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
 * @file testFixedDimensionalHydrostaticEquilibriumFixture.hpp
 *
 * Shared test fixture and TEST_P body for compositional multiphase flow
 * on fixed-dimensional (InternalMesh) meshes initialized with hydrostatic
 * equilibrium.  Both the serial and (optional) MPI test executables include
 * this header and only differ in the INSTANTIATE_TEST_SUITE_P parameters.
 *
 * Two cases are covered:
 *  - Two-phase   (water / gas,   2 components)
 *  - Three-phase (oil / water / gas, 3 components)
 *
 * The test verifies that:
 * 1. The Newton solver converges in at most @p maxExpectedNewtonIter iterations,
 *    because the hydrostatic equilibrium initialization guarantees that the system
 *    starts nearly at steady state.
 * 2. The discrete accumulation term (compAmount - compAmount_n) is zero for all
 *    conservation equations, confirming that the initial state is perfectly consistent.
 */

#ifndef GEOS_TEST_FIXEDDIMHYDROSTATICEQUILIBRIUM_FIXTURE_HPP
#define GEOS_TEST_FIXEDDIMHYDROSTATICEQUILIBRIUM_FIXTURE_HPP

#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <string>

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

constexpr real64 accumulation_tolerance = 1.0e-5; // same order as the newton tol in the input file

extern CommandLineOptions g_commandLineOptions;

/**
 * @brief Integration test for compositional multiphase flow on fixed-dimensional
 *        meshes (InternalMesh) initialized with hydrostatic equilibrium.
 *
 * Parametrized with std::string that selects the case:
 *  - "two_phases"   : 2-component, 2-phase  (water, gas)
 *  - "three_phases" : 3-component, 3-phase  (oil, water, gas)
 */
class FixedDimensionalHydrostaticEquilibriumTest
  : public ::testing::TestWithParam< std::string >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  /**
   * @brief Generate the two-phase XML input string.
   */
  static std::string generateTwoPhaseXml()
  {
    return
      R"xml(<?xml version="1.0" ?>
<Problem>
  <Solvers gravityVector="{0.0, 0.0, -9.81}">
    <CompositionalMultiphaseFVM name="FlowSolver" logLevel="0" discretization="fluidTPFA" temperature="355.0"
      targetRelativePressureChangeInTimeStep="1" targetPhaseVolFractionChangeInTimeStep="1"
      targetRegions="{reservoir}" useTotalMassEquation="1" useMass="1">
      <NonlinearSolverParameters newtonMaxIter="40" newtonTol="1.0e-5" maxTimeStepCuts="8" />
      <LinearSolverParameters directParallel="0" reuseFactorization="0"/>
    </CompositionalMultiphaseFVM>
  </Solvers>

  <Mesh>
    <InternalMesh name="mesh1" elementTypes="{C3D8}" nx="{3}" ny="{3}" nz="{9}"
      xCoords="{0, 1}" yCoords="{0, 1}" zCoords="{0, 3}" cellBlockNames="{cb1}"/>
  </Mesh>

  <Geometry>
    <Box name="zposFace" xMin="{-1e9, -1e9, 2.6}" xMax="{1e9, 1e9, 3.01}"/>
  </Geometry>

  <NumericalMethods>
    <FiniteVolume><TwoPointFluxApproximation name="fluidTPFA"/></FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="reservoir" cellBlocks="{*}" materialList="{fluid, rock, relPerm_matrix}"/>
  </ElementRegions>

  <Constitutive>
    <InvariantImmiscibleFluid name="fluid" componentNames="{ C0, C1 }" phaseNames="{water, gas}" densities="{1000.0, 10.0}"
      componentMolarWeight="{18e-3, 16e-3}" viscosities="{ 1.0e-3, 1.0e-4 }"/>
    <BrooksCoreyRelativePermeability name="relPerm_matrix" phaseNames="{water, gas}"
      phaseMinVolumeFraction="{0.13, 0.06}" phaseRelPermExponent="{1.5, 1.5}" phaseRelPermMaxValue="{0.9, 0.9}"/>
    <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>
    <NullModel name="nullSolid"/>
    <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.01" referencePressure="0.0" compressibility="1e-9"/>
    <ConstantPermeability name="rockPerm" permeabilityComponents="{1e-14, 1e-14, 1e-14}"/>
  </Constitutive>

  <Functions>
    <TableFunction name="WaterCompositionalFractionVsElevation" coordinates="{0, 1, 1.01, 2, 2.01, 3}" values="{1, 1, 0, 0, 0, 0}"/>
    <TableFunction name="GasCompositionalFractionVsElevation" coordinates="{0, 1, 1.01, 2, 2.01, 3}" values="{0, 0, 1, 1, 1, 1}"/>
    <TableFunction name="InitTempTable" coordinates="{-1.0, 3.0}" values="{355.0, 355.0}"/>
  </Functions>

  <FieldSpecifications>
    <HydrostaticEquilibrium name="CompositionalHydrostaticEquilibrium" objectPath="ElementRegions"
      datumElevation="3.0" datumPressure="1e7" componentNames="{C0, C1}" phaseContacts="{1.0}"
      temperatureVsElevationTableName="InitTempTable" elevationIncrementInHydrostaticPressureTable="0.01"
      componentFractionVsElevationTableNames="{WaterCompositionalFractionVsElevation, GasCompositionalFractionVsElevation}"/>
    <FieldSpecification name="topP" fieldName="pressure" setNames="{zposFace}" objectPath="faceManager" scale="1e7"/>
    <FieldSpecification name="topWaterComp" fieldName="globalCompFraction" component="0" setNames="{zposFace}" objectPath="faceManager" scale="0.0"/>
    <FieldSpecification name="topGasComp" fieldName="globalCompFraction" component="1" setNames="{zposFace}" objectPath="faceManager" scale="1.0"/>
    <FieldSpecification name="topT" fieldName="temperature" setNames="{zposFace}" objectPath="faceManager" scale="355.0"/>
  </FieldSpecifications>

  <Events minTime="0.0" maxTime="1.0">
    <PeriodicEvent name="solverApplications0" forceDt="1.0" endTime="1.0" target="/Solvers/FlowSolver"/>
  </Events>
</Problem>
)xml";
  }

  /**
   * @brief Generate the three-phase XML input string.
   */
  static std::string generateThreePhaseXml()
  {
    return
      R"xml(<?xml version="1.0" ?>
<Problem>
  <Solvers gravityVector="{0.0, 0.0, -9.81}">
    <CompositionalMultiphaseFVM name="FlowSolver" logLevel="0" discretization="fluidTPFA" temperature="355.0"
      targetRelativePressureChangeInTimeStep="1" targetPhaseVolFractionChangeInTimeStep="1"
      targetRegions="{reservoir}" useTotalMassEquation="1" useMass="1">
      <NonlinearSolverParameters newtonMaxIter="40" newtonTol="1.0e-5" maxTimeStepCuts="8" />
      <LinearSolverParameters directParallel="0" reuseFactorization="0"/>
    </CompositionalMultiphaseFVM>
  </Solvers>

  <Mesh>
    <InternalMesh name="mesh1" elementTypes="{C3D8}" nx="{3}" ny="{3}" nz="{9}"
      xCoords="{0, 1}" yCoords="{0, 1}" zCoords="{0, 3}" cellBlockNames="{cb1}"/>
  </Mesh>

  <Geometry>
    <Box name="zposFace" xMin="{-1e9, -1e9, 2.6}" xMax="{1e9, 1e9, 3.01}"/>
  </Geometry>

  <NumericalMethods>
    <FiniteVolume><TwoPointFluxApproximation name="fluidTPFA"/></FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="reservoir" cellBlocks="{*}" materialList="{fluid, rock, relPerm_matrix}"/>
  </ElementRegions>

  <Constitutive>
    <InvariantImmiscibleFluid name="fluid" componentNames="{ C0, C1, C2 }" phaseNames="{oil, water, gas}" densities="{600, 1000.0, 10.0}"
      componentMolarWeight="{114e-3, 18e-3, 16e-3}" viscosities="{ 1.0e-2, 1.0e-3, 1.0e-4 }"/>
    <BrooksCoreyRelativePermeability name="relPerm_matrix" phaseNames="{oil, water, gas}"
      phaseMinVolumeFraction="{0.08, 0.13, 0.06}" phaseRelPermExponent="{1.5, 1.5, 1.5}" phaseRelPermMaxValue="{0.9, 0.9, 0.9}"/>
    <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>
    <NullModel name="nullSolid"/>
    <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.01" referencePressure="0.0" compressibility="1e-9"/>
    <ConstantPermeability name="rockPerm" permeabilityComponents="{1e-14, 1e-14, 1e-14}"/>
  </Constitutive>

  <Functions>
    <TableFunction name="OilCompositionalFractionVsElevation" coordinates="{0, 1, 1.01, 2, 2.01, 3}" values="{0, 0, 1, 1, 0, 0}"/>
    <TableFunction name="WaterCompositionalFractionVsElevation" coordinates="{0, 1, 1.01, 2, 2.01, 3}" values="{1, 1, 0, 0, 0, 0}"/>
    <TableFunction name="GasCompositionalFractionVsElevation" coordinates="{0, 1, 1.01, 2, 2.01, 3}" values="{0, 0, 0, 0, 1, 1}"/>
    <TableFunction name="InitTempTable" coordinates="{-1.0, 3.0}" values="{355.0, 355.0}"/>
  </Functions>

  <FieldSpecifications>
    <HydrostaticEquilibrium name="CompositionalHydrostaticEquilibrium" objectPath="ElementRegions"
      datumElevation="3.0" datumPressure="1e7" componentNames="{C0, C1, C2}" phaseContacts="{1.0, 2.0}"
      temperatureVsElevationTableName="InitTempTable" elevationIncrementInHydrostaticPressureTable="0.01"
      componentFractionVsElevationTableNames="{OilCompositionalFractionVsElevation, WaterCompositionalFractionVsElevation, GasCompositionalFractionVsElevation}"/>
    <FieldSpecification name="topP" fieldName="pressure" setNames="{zposFace}" objectPath="faceManager" scale="1e7"/>
    <FieldSpecification name="topOilComp" fieldName="globalCompFraction" component="0" setNames="{zposFace}" objectPath="faceManager" scale="0.0"/>
    <FieldSpecification name="topWaterComp" fieldName="globalCompFraction" component="1" setNames="{zposFace}" objectPath="faceManager" scale="0.0"/>
    <FieldSpecification name="topGasComp" fieldName="globalCompFraction" component="2" setNames="{zposFace}" objectPath="faceManager" scale="1.0"/>
    <FieldSpecification name="topT" fieldName="temperature" setNames="{zposFace}" objectPath="faceManager" scale="355.0"/>
  </FieldSpecifications>

  <Events minTime="0.0" maxTime="1.0">
    <PeriodicEvent name="solverApplications0" forceDt="1.0" endTime="1.0" target="/Solvers/FlowSolver"/>
  </Events>
</Problem>
)xml";
  }

  std::string testBinaryDir;
};

TEST_P( FixedDimensionalHydrostaticEquilibriumTest, Run )
{
  std::string const & caseName = GetParam();

  // Select the XML content based on the case name.
  std::string xmlContent;
  if( caseName == "two_phases" )
  {
    xmlContent = generateTwoPhaseXml();
  }
  else if( caseName == "three_phases" )
  {
    xmlContent = generateThreePhaseXml();
  }
  else
  {
    FAIL() << "Unknown case name: " << caseName;
  }

  // Build a unique XML filename per test case.
  std::string const xmlName = "test_fixed_dim_hydrostatic_" + caseName;
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

    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();
    state.run();

    ProblemManager & pm = state.getProblemManager();

    // -----------------------------------------------------------------------
    // Check 1: Newton convergence
    // On a regular InternalMesh the hydrostatic discretization should be exact
    // (or nearly so), so only a small number of Newton iterations are expected.
    // -----------------------------------------------------------------------
    {
      CompositionalMultiphaseFVM & solver =
        pm.getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "FlowSolver" );

      integer const maxExpectedIter = 2;

      integer const numNewtonIter = solver.getNonlinearSolverParameters().m_numNewtonIterations;
      EXPECT_LE( numNewtonIter, maxExpectedIter )
        << "Case " << caseName
        << ": expected at most " << maxExpectedIter
        << " Newton iterations after hydrostatic init, got " << numNewtonIter;
    }

    // -----------------------------------------------------------------------
    // Check 2: Zero accumulation
    // compAmount == compAmount_n for every element and every component,
    // confirming that the initial state is a true discrete equilibrium.
    // -----------------------------------------------------------------------
    {
      MeshLevel & mesh = pm.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();

      mesh.getElemManager().forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
      {
        bool const isMatrixCell = dynamic_cast< CellElementSubRegion const * >( &subRegion );
        if( !isMatrixCell )
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
            real64 const accumulation = compAmount[k][ic] - compAmount_n[k][ic];
            EXPECT_NEAR( accumulation, 0.0, accumulation_tolerance )
              << "Case " << caseName
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

#endif /* GEOS_TEST_FIXEDDIMHYDROSTATICEQUILIBRIUM_FIXTURE_HPP */
