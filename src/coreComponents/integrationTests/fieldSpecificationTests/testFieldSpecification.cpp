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

static const string part1 =
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
      targetRegions="{ Channel1  }">
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

    <Box
      name="aquifer"
      xMin="{ -0.01, 3.99, 3.99 }"
      xMax="{  0.01, 5.01, 6.01 }"/> 
)xml";

static const string part2 =
  R"xml(
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
      name="Channel1"
      cellBlocks="{ cb-1_0_0, cb-0_0_0, cb-0_0_1, cb-0_1_1, cb-0_1_2, cb-1_1_2, cb-1_1_3 }"
      materialList="{ water, rock }"/>

    <CellElementRegion
      name="Channel2"
      cellBlocks="{ cb-1_0_3 }"
      materialList="{  }"/>

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
      setNames="{ all }"
      objectPath="ElementRegions/Channel1"
      errorSetMode="silent"
      fieldName="pressure"
      scale="0.0"/>
    )xml";

static const string part3 =
  R"xml(
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

void setupAndPlayWrongBC( string const & xml )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  problem.parseInputString( xml );
  problem.problemSetup();
  problem.applyInitialConditions();
}

std::vector< string > splitStringByDelimiter( string const & s, string const & delimiter )
{
  string cpyTokens = s;
  std::vector< string > tokens;
  size_t pos = 0;
  string token;
  while((pos = cpyTokens.find( delimiter )) != string::npos )
  {
    token = cpyTokens.substr( 0, pos );
    tokens.push_back( token );
    cpyTokens.erase( 0, pos + delimiter.length());
  }
  tokens.push_back( cpyTokens );

  return tokens;
}

TEST( testIncorrectFieldSpecification, testWrongFieldNames )
{
  static const string xmlWrongFieldNames =
    R"xml(
    <FieldSpecification
      name="sourceTerm"
      objectPath="ElementRegions/Channel1"
      errorSetMode="silent"
      fieldName="pdfd"
      scale="5e6"
      setNames="{ source }"/>
  )xml";

  static const string boxComponent =
    R"xml(
    <Box
      name="sink"
      xMin="{ 7.99, -0.01, -0.01 }"
      xMax="{ 10.01, 2.01, 0.01 }"/>
  )xml";

  static constexpr auto expectedMsg1 =  "Available fields in ElementRegions/Channel1 are:";
  static constexpr auto expectedMsg2 =  "{ channelPerm_dPerm_dPressure, channelPerm_permeability, deltaPressure, elementCenter, elementVolume, "
                                        "ghostRank, localToGlobalMap, mass, pressure, rockPorosity_initialPorosity, rockPorosity_porosity, "
                                        "rockPorosity_referencePorosity, temperature, water_dDensity, water_dEnthalpy, water_dInternalEnergy, "
                                        "water_dViscosity, water_density, water_enthalpy, water_internalEnergy, water_viscosity }";
  std::ostringstream xmlInput;
  xmlInput << part1 << boxComponent <<  part2 << xmlWrongFieldNames << part3;

  bool trowHappened = false;
  try
  {
    setupAndPlayWrongBC( xmlInput.str() );
  }
  catch( std::exception const & e )
  {
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg1 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg1 ) );
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg2 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg2 ) );
    trowHappened = true;
  }
  ASSERT_TRUE( trowHappened );
}

TEST( testIncorrectFieldSpecification, testRightFieldNames )
{

  static const string boxComponent =
    R"xml(
    <Box
      name="sink"
      xMin="{ 7.99, -0.01, -0.01 }"
      xMax="{ 10.01, 2.01, 0.01 }"/>
  )xml";
  static constexpr auto tokens =  "channelPerm_dPerm_dPressure, channelPerm_permeability, deltaPressure, elementCenter, elementVolume, "
                                  "ghostRank, localToGlobalMap, mass, pressure, rockPorosity_initialPorosity, rockPorosity_porosity, "
                                  "rockPorosity_referencePorosity, temperature, water_dDensity, water_dEnthalpy, water_dInternalEnergy, "
                                  "water_dViscosity, water_density, water_enthalpy, water_internalEnergy, water_viscosity";
  std::vector< string > const splitToken = splitStringByDelimiter( tokens, "," );
  for( auto const & token : splitToken )
  {
    string const xmlTemplate = R"xml(
    <FieldSpecification
      name="sourceTerm"
      objectPath="ElementRegions/Channel1"
      errorSetMode="silent"
      fieldName=")xml" + string( stringutilities::trimSpaces( token )) + R"xml("
      scale="5e6"
      setNames="{ source }"/>
    )xml";

    std::ostringstream xmlInput;
    xmlInput << part1 << boxComponent << part2 << xmlTemplate << part3;

    EXPECT_NO_THROW( setupAndPlayWrongBC( xmlInput.str() ));
  }
}

TEST( testIncorrectFieldSpecification, testSetNameTargetNothing )
{
  string const boxComponent = R"xml(
    <Box
        name="sink"
        xMin="{ 20, 20, -5 }"
        xMax="{ 21, 21, -6 }"/>
    )xml";

  std::ostringstream xmlInput;
  xmlInput << part1 << boxComponent << part2 << part3;

  bool trowHappened = false;
  try
  {
    setupAndPlayWrongBC( xmlInput.str() );
  }

  catch( std::exception const & e )
  {
    string const & expectedMsg1 = "does not capture anything in the mesh";
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg1 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg1 ) );
    trowHappened = true;
  }

  ASSERT_TRUE( trowHappened );
}

TEST( testIncorrectFieldSpecification, testSetName )
{

  static string const boxComponent = R"xml(
    <Box
        name="sink2"
        xMin="{ 2.99, -0.01, -0.01 }"
        xMax="{ 8.99, 5.01, 1.01 }"/>
    )xml";

  static string const fsComponent =
    R"xml(
    <FieldSpecification
      name="bcTemp"
      objectPath="ElementRegions/Channel2"
      fieldName="temperature"
      scale="-5e6"
      setNames="{ sink2 }"/>
      )xml";

  std::ostringstream xmlInput;
  xmlInput << part1 << boxComponent<< part2 << fsComponent << part3;

  bool trowHappened = false;
  try
  {
    setupAndPlayWrongBC( xmlInput.str() );
  }
  catch( std::exception const & e )
  {
    string const & expectedMsg1 = "- Does not capture: ElementRegions/Channel2";
    string const & expectedMsg2 = "- Instead, captures: nodeManager, edgeManager, faceManager";
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg1 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg1 ) );
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg2 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg2 ) );
    trowHappened = true;
  }

  ASSERT_TRUE( trowHappened );
}

TEST( testIncorrectFieldSpecification, testMultiSetNames )
{
  static string const boxComponent1 = R"xml(
   <Box
      name="sink"
      xMin="{ 0.99, -0.01, -0.01 }"
      xMax="{ 10.01, 10.01, 2.01 }"/>
    )xml";

  static string const fsComponent =
    R"xml(
    <FieldSpecification
      name="bcTemp"
      initialCondition="1"
      objectPath="ElementRegions/Channel1"
      fieldName="pressure"
      scale="0.0"
      setNames="{ sink, xneg }"/>
      )xml";

  std::ostringstream xmlInput;
  xmlInput << part1 <<boxComponent1<< part2 << fsComponent << part3;

  bool trowHappened = false;
  try
  {
    setupAndPlayWrongBC( xmlInput.str() );
  }
  catch( std::exception const & e )
  {
    string const & expectedMsg1 = "Set 'xneg':";
    string const & expectedMsg2 = "- Instead, captures: nodeManager, edgeManager, faceManager";
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg1 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg1 ) );
    GEOS_ERROR_IF_EQ_MSG( string( e.what() ).find( expectedMsg2 ), string::npos,
                          GEOS_FMT( "The error message was not containing the expected sequence.\n"
                                    "  Error message :\n{}"
                                    "  expected sequence :\n{}",
                                    e.what(),
                                    expectedMsg2 ) );
    trowHappened = true;
  }

  ASSERT_TRUE( trowHappened );
}



int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();
  return result;
}
