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

static const string header =
  R"xml(
    <?xml version="1.0" ?>

<Problem xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:noNamespaceSchemaLocation="../../src/coreComponents/schema/schema.xsd">
  <Solvers
    gravityVector="{ 0.0, 0.0, -9.81 }">
    <SinglePhaseFVM
      name="SinglePhaseFlow"
      discretization="singlePhaseTPFA"
      logLevel="0"
      targetRegions="{ Channel  }">
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
      nx="{ 2, 2 }"
      ny="{ 2, 2 }"
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
        )xml";

static const string footer =
  R"xml(
    <Box
      name="aquifer"
      xMin="{ -0.01, 3.99, 3.99 }"
      xMax="{  0.01, 5.01, 6.01 }"/> 
    
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
      cellBlocks="{ cb-1_0_0, cb-1_0_3, cb-0_0_0, cb-0_0_1, cb-0_1_1, cb-0_1_2, cb-1_1_2, cb-1_1_3 }"
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

  <FieldSpecifications>
    <FieldSpecification
      name="initialPressure_channel"
      initialCondition="1"
      objectPath="ElementRegions/Channel"
      setNames="{ all }"
      fieldName="pressure"
      scale="0.0"/>

    <FieldSpecification
      name="sinkTerm"
      objectPath="faceManager"
      fieldName="pressure"
      scale="-5e6"
      setNames="{ sink }"/>

    <Aquifer
      name="aquiferBC"
      aquiferPorosity="2e-1"
      aquiferPermeability="3e-13"
      aquiferInitialPressure="6e6"
      aquiferWaterViscosity="0.00089"
      aquiferWaterDensity="962.81"
      aquiferTotalCompressibility="1e-10"
      aquiferElevation="5"
      aquiferThickness="18"
      aquiferInnerRadius="2000"
      aquiferAngle="20"
      setNames="{ aquifer }"/>

  </FieldSpecifications>

  <Tasks>
  </Tasks>
  
  <Outputs>
  </Outputs>
</Problem>

)xml";

CommandLineOptions g_commandLineOptions;

void setupAndPlayWrongBC( string const & wrongXml )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  std::ostringstream xmlInput;
  xmlInput << header << wrongXml << footer;
  problem.parseInputString( xmlInput.str() );
  problem.problemSetup();
  problem.applyInitialConditions();
}

TEST( testIncorrectFieldSpecification, testWrongSetNames )
{
  string const xmlTemplate = R"xml(
    <Box
        name="sink"
        xMin="{ 20, 20, -5 }"
        xMax="{ 21, 21, -6 }"/>
    )xml";

  bool trowHappened = false;
  try
  {
    setupAndPlayWrongBC( xmlTemplate );
  }

  catch( std::exception const & e )
  {
    string const & expectedMsg1 = "does not capture anything in the mesh";
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg1 ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" << e.what() <<
                          "  expected sequence :\n" << expectedMsg1 );
    trowHappened = true;
  }

  ASSERT_TRUE( trowHappened );
}

TEST( testIncorrectFieldSpecification, testSetNames )
{
  string const xmlTemplate = R"xml(
    <Box
        name="sink"
        xMin="{ 2.99, -0.01, -0.01 }"
        xMax="{ 8.99, 5.01, 1.01 }"/>
    )xml";

  bool trowHappened = false;
  try
  {
    setupAndPlayWrongBC( xmlTemplate );
  }
  catch( std::exception const & e )
  {
    string const & expectedMsg1 = "- Does not capture: faceManager";
    string const & expectedMsg2 = "- Instead, captures: edgeManager, nodeManager";
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg1 ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" << e.what() <<
                          "  expected sequence :\n" << expectedMsg1 );
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg2 ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" << e.what() <<
                          "  expected sequence :\n" << expectedMsg2 );
    trowHappened = true;
  }

  ASSERT_TRUE( trowHappened );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  basicCleanup();
  return result; 
}
