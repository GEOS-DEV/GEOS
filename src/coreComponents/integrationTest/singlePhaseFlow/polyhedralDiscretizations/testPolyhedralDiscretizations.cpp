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
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Define the XML input for the test
char const * xmlInput =
  R"xml(
  <Problem>

  <Mesh>
      <VTKMesh
        name="mesh"
        logLevel="5"  
        partitionRefinement="0"
        useGlobalIds="0"
        file="polyhedral_extruded.vtk"/>
    </Mesh>

    <Geometry>
        <Box
        name="westBC"
        xMin="{ -0.001, 0.0, 0.0}"
        xMax="{ +0.001, 1.0, 1.0}"/>
        <Box
        name="eastBC"
        xMin="{ +0.999, 0.0, 0.0}"
        xMax="{ +1.001, 1.0, 1.0}"/>
    </Geometry>

    <ElementRegions>
      <CellElementRegion
        name="Domain"
        cellBlocks="{ * }"
        materialList="{rock, fluid }"/>
    </ElementRegions>

    <Solvers gravityVector="{ 0.0, 0.0, 0.0}"> </Solvers>

    <Constitutive>

      <CompressibleSinglePhaseFluid
        name="fluid"
        defaultDensity="1000"
        defaultViscosity="0.001"
        referencePressure="0.0"
        compressibility="0.0"
        viscosibility="0.0"/>

      <CompressibleSolidConstantPermeability
        name="rock"
        solidModelName="nullSolid"
        porosityModelName="rockPorosity"
        permeabilityModelName="rockPerm"/>

      <NullModel
        name="nullSolid"/>

      <PressurePorosity
        name="rockPorosity"
        defaultReferencePorosity="0.1"
        referencePressure="0.0"
        compressibility="0.0"/>

      <ConstantPermeability
        name="rockPerm"
        permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>

    </Constitutive>

    <FieldSpecifications>

      <FieldSpecification
        name="initialPressure"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/Domain"
        fieldName="pressure"
        scale="1.0e7"/>    
      <FieldSpecification
        name="west_pressure"
        setNames="{ westBC }"
        objectPath="faceManager"
        fieldName="pressure"
        scale="2.0e7" />
      <FieldSpecification
        name="east_pressure"
        setNames="{ eastBC }"
        objectPath="faceManager"
        fieldName="pressure"
        scale="1.0e7" />      

    </FieldSpecifications>

    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation
          name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>

    <Solvers>
      <SinglePhaseFVM
        name="SinglePhaseFlow"
        logLevel="1"
        discretization="singlePhaseTPFA"
        targetRegions="{ Domain }">
        <NonlinearSolverParameters
          newtonTol="1.0e-5"
          newtonMaxIter="2"/>
        <LinearSolverParameters
          directParallel="0"/>
      </SinglePhaseFVM>
    </Solvers>

    <Events
      minTime="0.0"
      maxTime="86400">
      <PeriodicEvent
        name="outputs"
        timeFrequency="86400"
        target="/Outputs/vtkConsistencyTPFA"/>
      <PeriodicEvent
        name="solverApplications"
        endTime="86400"
        maxEventDt="86400"
        target="/Solvers/SinglePhaseFlow"/>
    </Events>

    <Outputs>
      <VTK
        name="vtkConsistencyTPFA"/>
    </Outputs>

  </Problem>
  )xml";

class TPFAIntegrationTest : public ::testing::Test {

public:
  TPFAIntegrationTest() : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override {
    // Setup problem from XML input
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  GeosxState state;
};

TEST_F(TPFAIntegrationTest, PressureFieldL2Error) {
  ProblemManager & problemManager = state.getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();
  // Add test logic here to validate the pressure field or other properties
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  g_commandLineOptions = *geos::basicSetup(argc, argv);
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
