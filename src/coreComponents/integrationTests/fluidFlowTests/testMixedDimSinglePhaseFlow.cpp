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
#include "codingUtilities/Utilities.hpp"

using namespace geos;

CommandLineOptions g_commandLineOptions;

class MixedDimSinglePhaseFlowTest : public ::testing::TestWithParam< std::string >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  std::string generateXmlInput( std::string const & meshFile, std::string const & nodeSetNames )
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
    </Geometry>

  <Solvers>
    <SinglePhaseFVM
      name="flowSolver"
      targetRegions="{ Region, Fracture }"
      discretization="tpfa"
      logLevel="1">
      <NonlinearSolverParameters newtonTol="1.0e-6" newtonMaxIter="20"/>
      <LinearSolverParameters solverType="gmres" preconditionerType="amg" krylovTol="1.0e-10"/>
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
    <SurfaceElementRegion name="Fracture" faceBlock="faceElementSubRegion" materialList="{ fractureMaterial, water }" defaultAperture="1.0e-5"/>
  </ElementRegions>

  <Constitutive>
  
    <CompressibleSinglePhaseFluid
      name="water"
      defaultDensity="1000"
      defaultViscosity="0.001"
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

    <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0e-15, 1.0e-15, 1.0e-15 }"/>
    <ConstantPermeability name="fracPerm" permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>      

    <NullModel
      name="nullSolid"/>

    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="0.05"
      referencePressure="0.0"
      compressibility="0.0"/>

    <PressurePorosity
      name="fracPorosity"
      defaultReferencePorosity="1.0"
      referencePressure="0.0"
      compressibility="0.0"/>      

  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification name="separableFace" fieldName="isFaceSeparable" initialCondition="1" setNames=")xml" << nodeSetNames << R"xml(" objectPath="faceManager" scale="1" />
    <FieldSpecification name="frac" initialCondition="1" setNames=")xml" << nodeSetNames << R"xml(" objectPath="faceManager" fieldName="ruptureState" scale="1" />
    <FieldSpecification name="initialP" fieldName="pressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions" scale="15.0e6"/>
    <FieldSpecification name="initialPf" fieldName="pressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions/Fracture/faceElementSubRegion" scale="20.0e6"/>
    <FieldSpecification name="inletP" fieldName="pressure" setNames="{ xnegFace }" objectPath="faceManager" scale="20.0e6"/>
    <FieldSpecification name="outletP" fieldName="pressure" setNames="{ xposFace }" objectPath="faceManager" scale="10.0e6"/>
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
    <PeriodicEvent name="solverApplications" target="/Solvers/flowSolver" forceDt="1.0"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( MixedDimSinglePhaseFlowTest, Run )
{
  std::string const & meshFileName = GetParam();
  std::string xmlPath = testBinaryDir + "/test_mixed_dim_single_phase_flow.xml";

  std::string nodeSetNames = "{ f1_node_set }";
  if( meshFileName.find( "_DFN_2.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f2_node_set }";
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
  else if( meshFileName.find( "_DFN_123.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set, f3_node_set }";
  }

  std::string xmlContent = generateXmlInput( meshFileName, nodeSetNames );

  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }

  auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
  options->inputFileNames.push_back( xmlPath );
  options->problemName = "test_mixed_dim_single_phase_flow";

  // Scoped state to ensure cleanup
  {
    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();
    state.run();
  }
}

INSTANTIATE_TEST_SUITE_P(
  MixedDimFlowCases,
  MixedDimSinglePhaseFlowTest,
  ::testing::Values(
    "fractured_mesh_tet_DFN_12.vtu"
  )
);

int main( int argc, char * argv[] )
{
  MPI_Init( &argc, &argv );
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv, false );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  MPI_Finalize();
  return result;
}
