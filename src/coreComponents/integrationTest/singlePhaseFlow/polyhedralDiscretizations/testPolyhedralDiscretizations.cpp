/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 * ------------------------------------------------------------------------------------------------------------
 */

#include "gtest/gtest.h"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "finiteVolume/SinglePhaseFVM.hpp"
#include "finiteVolume/TwoPointFluxApproximation.hpp"
#include "constitutive/CompressibleSinglePhaseFluid.hpp"
#include "constitutive/CompressibleSolidConstantPermeability.hpp"
#include "constitutive/NullModel.hpp"
#include "constitutive/PressurePorosity.hpp"
#include "constitutive/ConstantPermeability.hpp"

using namespace geos;

constexpr double maxCoordInX = 1.0;
constexpr double maxCoordInY = 1.0;
constexpr double maxCoordInZ = 1.0;
constexpr localIndex numElemsInX = 10;
constexpr localIndex numElemsInY = 10;
constexpr localIndex numElemsInZ = 10;

class TPFAIntegrationTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Setup mesh and problem from XML string (unit cube, TPFA, single phase)
    string const inputStream = GEOS_FMT(
      "<Problem>"
      "  <Mesh>"
      "    <InternalMesh name=\"mesh1\" elementTypes=\"{{C3D8}}\" xCoords=\"{{0,1}}\" yCoords=\"{{0,1}}\" zCoords=\"{{0,1}}\" nx=\"{{10}}\" ny=\"{{10}}\" nz=\"{{10}}\" cellBlockNames=\"{{cb1}}\"/>"
      "  </Mesh>"
      "  <ElementRegions>"
      "    <CellElementRegion name=\"Domain\" cellBlocks=\"{{cb1}}\" materialList=\"{{rock,fluid}}\"/>"
      "  </ElementRegions>"
      "  <NumericalMethods>"
      "    <FiniteVolume>"
      "      <TwoPointFluxApproximation name=\"singlePhaseTPFA\"/>"
      "    </FiniteVolume>"
      "  </NumericalMethods>"
      "  <Solvers>"
      "    <SinglePhaseFVM name=\"SinglePhaseFlow\" logLevel=\"1\" discretization=\"singlePhaseTPFA\" targetRegions=\"{Domain}\">"
      "      <NonlinearSolverParameters newtonTol=\"1.0e-5\" newtonMaxIter=\"2\"/>"
      "      <LinearSolverParameters directParallel=\"0\"/>"
      "    </SinglePhaseFVM>"
      "  </Solvers>"
      "  <Constitutive>"
      "    <CompressibleSinglePhaseFluid name=\"fluid\" defaultDensity=\"1000\" defaultViscosity=\"0.001\" referencePressure=\"0.0\" compressibility=\"0.0\" viscosibility=\"0.0\"/>"
      "    <CompressibleSolidConstantPermeability name=\"rock\" solidModelName=\"nullSolid\" porosityModelName=\"rockPorosity\" permeabilityModelName=\"rockPerm\"/>"
      "    <NullModel name=\"nullSolid\"/>"
      "    <PressurePorosity name=\"rockPorosity\" defaultReferencePorosity=\"0.1\" referencePressure=\"0.0\" compressibility=\"0.0\"/>"
      "    <ConstantPermeability name=\"rockPerm\" permeabilityComponents=\"{1.0e-13,1.0e-13,1.0e-13}\"/>"
      "  </Constitutive>"
      "  <FieldSpecifications>"
      "    <FieldSpecification name=\"initialPressure\" initialCondition=\"1\" setNames=\"{all}\" objectPath=\"ElementRegions/Domain\" fieldName=\"pressure\" scale=\"1.0e7\"/>"
      "    <FieldSpecification name=\"west_pressure\" setNames=\"{westBC}\" objectPath=\"faceManager\" fieldName=\"pressure\" scale=\"2.0e7\"/>"
      "    <FieldSpecification name=\"east_pressure\" setNames=\"{eastBC}\" objectPath=\"faceManager\" fieldName=\"pressure\" scale=\"1.0e7\"/>"
      "  </FieldSpecifications>"
      "  <Geometry>"
      "    <Box name=\"westBC\" xMin=\"{-0.001,0.0,0.0}\" xMax=\"{+0.001,1.0,1.0}\"/>"
      "    <Box name=\"eastBC\" xMin=\"{+0.999,0.0,0.0}\" xMax=\"{+1.001,1.0,1.0}\"/>"
      "  </Geometry>"
      "</Problem>");

    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult = xmlDocument.loadString(inputStream);
    ASSERT_TRUE(xmlResult);
    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild(dataRepository::keys::ProblemManager);
    ProblemManager & problemManager = getGlobalState().getProblemManager();
    problemManager.processInputFileRecursive(xmlDocument, xmlProblemNode);
    DomainPartition & domain = problemManager.getDomainPartition();
    MeshManager & meshManager = problemManager.getGroup<MeshManager>(problemManager.groupKeys.meshManager);
    meshManager.generateMeshLevels(domain);
    ElementRegionManager & elementManager = domain.getMeshBody(0).getBaseDiscretization().getElemManager();
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child(elementManager.getName().c_str());
    elementManager.processInputFileRecursive(xmlDocument, topLevelNode);
    elementManager.postInputInitializationRecursive();
    problemManager.problemSetup();
    problemManager.applyInitialConditions();
  }
};

TEST_F(TPFAIntegrationTest, PressureFieldL2Error) {
  ProblemManager & problemManager = getGlobalState().getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();
  MeshLevel & mesh = domain.getMeshBody(0).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion(0).getSubRegion<CellElementSubRegion>(0);
  localIndex numElems = subRegion.size();
  arrayView1d<real64 const> pressure = subRegion.getField<real64>("pressure");
  arrayView2d<real64 const> centers = subRegion.getElementCenter();
  // Reference solution: linear from 2e7 to 1e7 in x
  std::vector<real64> reference(numElems);
  for(localIndex elemID = 0; elemID < numElems; ++elemID) {
    reference[elemID] = 2.0e7 - centers[elemID][0] * 1.0e7;
  }
  // Compute L2 error
  double l2err = 0.0;
  for(localIndex elemID = 0; elemID < numElems; ++elemID) {
    l2err += std::pow(pressure[elemID] - reference[elemID], 2);
  }
  l2err = std::sqrt(l2err / numElems);
  std::cout << "TPFA L2 error: " << l2err << std::endl;
  EXPECT_GT(l2err, 1e-6) << "TPFA solution should be inexact";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  GeosxState state(geos::basicSetup(argc, argv));
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
