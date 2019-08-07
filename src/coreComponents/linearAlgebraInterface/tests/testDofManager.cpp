/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file testDofManager.cpp
 * @brief This test file is part of the ctest suite and tests the DofManager functionalities.
 */

#include "gtest/gtest.h"

#ifdef __clang__
#define __null nullptr
#endif

#include <numeric>
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/initialization.hpp"
#include "common/TimingMacros.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "createConnLocPattern.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

#include "DofManager.hpp"

using namespace geosx;

namespace
{
int global_argc;
char** global_argv;
}

class DofManagerTest : public ::testing::Test
{
public:
  /**
   * @brief Set the timer function.
   *
   * @param [out] time double actual time.
   */
  void setTimer( double &time )
  {
    time = MPI_Wtime();
  }

  /**
   * @brief Get elapsed time.
   *
   * @param [inout] time double
   * - on input: actual time from setTimer;
   * - on output: elapsed time.
   */
  void getElapsedTime( double &time )
  {
    time = MPI_Wtime() - time;
  }

protected:

  /**
   * @brief Provide a very simple xml input, with an internally generated mesh.
   */
  static void SetUpTestCase()
  {
    problemManager = new ProblemManager("Problem", nullptr);

    string const inputStream =
    "<?xml version=\"1.0\" ?>"
    "<Problem xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">"
    "  <Mesh>"
    "    <InternalMesh name=\"mesh1\""
    "                  elementTypes=\"{C3D8}\""
    "                  xCoords=\"{0, 1, 2, 3, 4}\""
    "                  yCoords=\"{0, 1}\""
    "                  zCoords=\"{0, 1}\""
    "                  nx=\"{4, 4, 4, 4}\""
    "                  ny=\"{4}\""
    "                  nz=\"{5}\""
    "                  cellBlockNames=\"{block1, block2, block3, block4}\"/>"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <ElementRegion name=\"region1\" cellBlocks=\"{block1}\" materialList=\"{dummy_material}\" />"
    "    <ElementRegion name=\"region2\" cellBlocks=\"{block2}\" materialList=\"{dummy_material}\" />"
    "    <ElementRegion name=\"region3\" cellBlocks=\"{block3}\" materialList=\"{dummy_material}\" />"
    "    <ElementRegion name=\"region4\" cellBlocks=\"{block4}\" materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";

    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
    if (!xmlResult)
    {
      GEOS_LOG_RANK_0("XML parsed with errors!");
      GEOS_LOG_RANK_0("Error description: " << xmlResult.description());
      GEOS_LOG_RANK_0("Error offset: " << xmlResult.offset);
    }

    int mpiSize = CommunicationTools::MPI_Size( MPI_COMM_GEOSX );
    dataRepository::ManagedGroup * commandLine =
      problemManager->GetGroup<dataRepository::ManagedGroup>( problemManager->groupKeys.commandLine );
    commandLine->RegisterViewWrapper<integer>( problemManager->viewKeys.xPartitionsOverride.Key() )->
      setApplyDefaultValue(mpiSize);

    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( "Problem" );
    problemManager->InitializePythonInterpreter();
    problemManager->ProcessInputFileRecursive( xmlProblemNode );

    // Open mesh levels
    DomainPartition * domain  = problemManager->getDomainPartition();
    MeshManager * meshManager = problemManager->GetGroup<MeshManager>( problemManager->groupKeys.meshManager );
    meshManager->GenerateMeshLevels(domain);

    ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager->getName().c_str() );
    elementManager->ProcessInputFileRecursive( topLevelNode );
    elementManager->PostProcessInputRecursive();

    problemManager->ProblemSetup();
  }

  /**
   * @brief Destructor.
   */
  static void TearDownTestCase()
  {
    delete problemManager;
    problemManager = nullptr;
  }

  static ProblemManager * problemManager;
};

ProblemManager * DofManagerTest::problemManager = nullptr; //!< the main problemManager.

/**
 * @name Ctest tests.
 * @brief Run testing functions using different DofManager functionalities.
 */

/**
 * @function TestFEM_partial.
 * @brief Define a vectorial nodal field (displacement) connected through elements (FEM) and check the
 * global number of DoFs. The field is defined on 3 regions out of 4.
 */
TEST_F(DofManagerTest, TestFEM_partial)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  string_array displacementRegion;
  displacementRegion.push_back( "region1" );
  displacementRegion.push_back( "region3" );
  displacementRegion.push_back( "region4" );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3,
                       displacementRegion );

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "displacement", "displacement" );

  // Total number of nodes, sum of regions 1 and 3+4
  constexpr globalIndex nNodes = 9*5*6 + 5*5*6;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes );
  EXPECT_EQ( pattern.globalCols(), 3*nNodes );
}

/**
 * @function TestFEM_all.
 * @brief Define a vectorial nodal field (displacement) connected through elements (FEM) and check the
 * global number of DoFs. The field is defined on all regions.
 */
TEST_F(DofManagerTest, TestFEM_all)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3 );

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "displacement", "displacement" );

  // Total number of nodes, sum of all regions
  constexpr globalIndex nNodes = 17*5*6;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes );
  EXPECT_EQ( pattern.globalCols(), 3*nNodes );
}

/**
 * @function TestFVM_partial.
 * @brief Define a scalar element field (pressure) connected through faces (FVM) and check the
 * global number of DoFs. The field is defined on 3 regions out of 4.
 */
TEST_F(DofManagerTest, TestFVM_partial)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  string_array pressureRegion;
  pressureRegion.push_back( "region1" );
  pressureRegion.push_back( "region2" );
  pressureRegion.push_back( "region3" );

  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );

  // Total number of cells
  constexpr globalIndex nCells = 12*4*5;

  EXPECT_EQ( pattern.globalRows(), nCells );
  EXPECT_EQ( pattern.globalCols(), nCells );
}

/**
 * @function TestFVM_all.
 * @brief Define a scalar element field (pressure) connected through faces (FVM) and check the
 * global number of DoFs. The field is defined on all regions.
 */
TEST_F(DofManagerTest, TestFVM_all)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );

  // Total number of cells
  constexpr globalIndex nCells = 16*4*5;

  EXPECT_EQ( pattern.globalRows(), nCells );
  EXPECT_EQ( pattern.globalCols(), nCells );
}

/**
 * @function TestCoupling
 * @brief Define two fields, one is displacement (FEM) and the other pressure (FVM), and a coupling
 * among them. Check the sizes of the coupling block pattern.
 */
TEST_F(DofManagerTest, TestCoupling)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  string_array displacementRegion;
  displacementRegion.push_back( "region1" );
  displacementRegion.push_back( "region3" );
  displacementRegion.push_back( "region4" );

  string_array pressureRegion;
  pressureRegion.push_back( "region1" );
  pressureRegion.push_back( "region2" );
  pressureRegion.push_back( "region3" );

  string_array couplingRegion;
  couplingRegion.push_back( "region3" );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3,
                       displacementRegion );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.addCoupling( "displacement", "pressure", DofManager::Connectivity::Elem, couplingRegion, false );

  ParallelMatrix pattern, emptyPattern;

  dofManager.setSparsityPattern( pattern, "displacement", "pressure" );
  dofManager.setSparsityPattern( emptyPattern, "pressure", "displacement" );

  // Total number of nodes
  constexpr globalIndex nNodes = 9*5*6 + 5*5*6;

  // Total number of cells
  constexpr globalIndex nCells = 12*4*5;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes );
  EXPECT_EQ( pattern.globalCols(), nCells );
  EXPECT_EQ( emptyPattern.globalRows(), nCells );
  EXPECT_EQ( emptyPattern.globalCols(), 3*nNodes );
  EXPECT_DOUBLE_EQ( emptyPattern.normFrobenius(), 0. );
}

/**
 * @function TestUserDefinedPattern
 * @brief Create an external pattern and pass it to the DofManager.
 */
TEST_F(DofManagerTest, TestUserDefinedPattern)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  ParallelMatrix connLocInput;
  // Create a TPFA Finite Volume stencil
  createConnLocPattern( domain, 0, 0, 1, connLocInput );

  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );
  dofManager.addField( "user-defined", connLocInput, DofManager::Connectivity::Face );

  ParallelMatrix pattern, userPattern;

  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );
  dofManager.setSparsityPattern( userPattern, "user-defined", "user-defined" );

  pattern.set( 1. );
  userPattern.set( 1. );

  EXPECT_EQ( pattern.globalRows(), userPattern.globalRows() );
  EXPECT_EQ( pattern.globalCols(), userPattern.globalCols() );
  EXPECT_DOUBLE_EQ( pattern.normFrobenius(), userPattern.normFrobenius() );
}

/**
 * @function TestFEM_FVM
 * @brief Define two fields, one is displacement (FEM) and the other pressure (FVM), and check the
 * size of the global pattern.
 */
TEST_F(DofManagerTest, TestFEM_FVM)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3 );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern );

  // Total number of nodes
  constexpr globalIndex nNodes = 17*5*6;

  // Total number of cells
  constexpr globalIndex nCells = 16*4*5;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes+nCells );
  EXPECT_EQ( pattern.globalCols(), 3*nNodes+nCells );
}

/**
 * @function TestMassMatrix
 * @brief Define a self-connected field, like a mass matrix.
 */
TEST_F(DofManagerTest, TestMassMatrix)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "massmatrix", DofManager::Location::Elem, DofManager::Connectivity::Elem );

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern );

  // Total number of cells
  constexpr globalIndex nCells = 16*4*5;

  pattern.set( 1. );

  EXPECT_EQ( pattern.globalRows(), nCells );
  EXPECT_EQ( pattern.globalCols(), nCells );
  EXPECT_DOUBLE_EQ( pattern.normFrobenius(), sqrt(nCells) );
}

/**
 * @function TestIndices
 * @brief Check the getIndices functions for FEM and FVM.
 */
TEST_F(DofManagerTest, TestIndices)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3 );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );

  constexpr globalIndex indicesDefaultDisp[24] = { 72, 73, 74, 75, 76, 77, 108, 109, 110,
    111, 112, 113, 78, 79, 80, 81, 82, 83, 114, 115, 116, 117, 118, 119 };
  constexpr globalIndex indicesDefaultPres[2] = { 1530, 1550 };

  globalIndex_array indicesDisp, indicesPres;
  dofManager.getIndices( indicesDisp, DofManager::Connectivity::Elem, 0, 0, 10, "displacement" );
  dofManager.getIndices( indicesPres, DofManager::Connectivity::Face, 3, "pressure" );

  int mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );
  if( mpiRank==0 )
  {
    for( localIndex i=0; i<std::min(indicesDisp.size(), integer_conversion<localIndex>(24)); ++i )
    {
      EXPECT_EQ( indicesDisp[i], indicesDefaultDisp[i] );
    }
    for( localIndex i=0; i<std::min(indicesPres.size(), integer_conversion<localIndex>(2)); ++i )
    {
      EXPECT_EQ( indicesPres[i], indicesDefaultPres[i] );
    }
    EXPECT_EQ( indicesDisp.size(), 24 );
    EXPECT_EQ( indicesPres.size(), 2 );
  }
  else
  {
    // Only rank 0 there will be always present, so it is the only with a real check!
    EXPECT_EQ( 1, 1 );
  }
}

/**
 * @function TestPermutation
 * @brief Check the permutation from the field-based ordering to the MPI rank-based one.
 */
TEST_F(DofManagerTest, TestPermutation)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  string_array displacementRegion;
  displacementRegion.push_back( "region1" );
  displacementRegion.push_back( "region3" );
  displacementRegion.push_back( "region4" );

  string_array pressureRegion;
  pressureRegion.push_back( "region1" );
  pressureRegion.push_back( "region2" );
  pressureRegion.push_back( "region3" );

  string_array couplingRegion;
  couplingRegion.push_back( "region3" );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3,
                       displacementRegion );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.addCoupling( "displacement", "pressure", DofManager::Connectivity::Elem, couplingRegion, false );

  ParallelMatrix pattern;

  dofManager.setSparsityPattern( pattern );

  // Create the permutation collecting all unknowns belonging to each process
  ParallelMatrix permutation;
  dofManager.createPermutation( permutation );

  ParallelMatrix permutedPattern;
  // Apply the permutation
  dofManager.permuteSparsityPattern( pattern, permutation, permutedPattern );

  EXPECT_EQ( pattern.globalRows(), permutedPattern.globalRows() );
  EXPECT_EQ( pattern.globalCols(), permutedPattern.globalCols() );
  EXPECT_DOUBLE_EQ( pattern.normFrobenius(), permutedPattern.normFrobenius() );
}

/**
 * @function TestWithTimes
 * @brief This is not a real test (it is always true!). It collects timings.
 */
TEST_F(DofManagerTest, TestWithTimes)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  const bool printPattern = false;

  string_array displacementRegion;
  displacementRegion.push_back( "region1" );
  displacementRegion.push_back( "region3" );
  displacementRegion.push_back( "region4" );

  string_array pressureRegion;
  pressureRegion.push_back( "region1" );
  pressureRegion.push_back( "region2" );
  pressureRegion.push_back( "region3" );

  string_array testRegion1;
  testRegion1.push_back( "region1" );

  string_array testRegion2;
  testRegion2.push_back( "region2" );
  testRegion2.push_back( "region4" );

  string_array testRegion3;
  testRegion3.push_back( "region3" );

  ParallelMatrix connLocInput;
  createConnLocPattern( domain, 0, 0, 1, connLocInput );

  double timeAddField, timeAddCoupling, timeGetSingleSparsityPattern, timeGetGlobalSparsityPattern, timeGetIndices;

  setTimer( timeAddField );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3,
                       displacementRegion );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.addField( "massmatrix", DofManager::Location::Elem, DofManager::Connectivity::None, 2, testRegion3 );
  dofManager.addField( "user-defined", connLocInput, DofManager::Connectivity::Face );

  getElapsedTime( timeAddField );
  setTimer( timeAddCoupling );

  dofManager.addCoupling( "displacement", "pressure", DofManager::Connectivity::Elem, testRegion3, false );
  dofManager.addCoupling( "massmatrix", "pressure", DofManager::Connectivity::Elem );

  getElapsedTime( timeAddCoupling );
  setTimer( timeGetSingleSparsityPattern );

  ParallelMatrix pattern;

  dofManager.setSparsityPattern( pattern, "displacement", "displacement" );
  if( printPattern )
  {
    pattern.write( "displacement.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );
  if( printPattern )
  {
    pattern.write( "pressure.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "massmatrix", "massmatrix" );
  if( printPattern )
  {
    pattern.write( "massmatrix.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "displacement", "pressure" );
  if( printPattern )
  {
    pattern.write( "coupling1.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "pressure", "displacement" );
  if( printPattern )
  {
    pattern.write( "coupling1_empty.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "pressure", "massmatrix" );
  if( printPattern )
  {
    pattern.write( "coupling2.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "massmatrix", "pressure" );
  if( printPattern )
  {
    pattern.write( "coupling2_transp.mtx" );
  }

  dofManager.setSparsityPattern( pattern, "user-defined", "user-defined" );
  if( printPattern )
  {
    pattern.write( "user-defined.mtx" );
  }

  getElapsedTime( timeGetSingleSparsityPattern );
  setTimer( timeGetGlobalSparsityPattern );

  dofManager.setSparsityPattern( pattern );
  if( printPattern )
  {
    pattern.write( "global.mtx" );
  }

  // Create the permutation collecting all unknowns belonging to each process
  ParallelMatrix permutation;
  dofManager.createPermutation( permutation );

  ParallelMatrix permutedPattern;
  // Apply the permutation
  dofManager.permuteSparsityPattern( pattern, permutation, permutedPattern );
  if( printPattern )
  {
    pattern.write( "permutatedGlobal.mtx" );
  }

  getElapsedTime( timeGetGlobalSparsityPattern );
  setTimer( timeGetIndices );

  int mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );

  globalIndex_array indices;
  dofManager.getIndices( indices, DofManager::Connectivity::Elem, 1, 0, 10, "displacement" );
  if( indices.size() > 0 )
    std::cout << mpiRank << " - " << "displacement" << " " << indices << std::endl;

  dofManager.getIndices( indices, DofManager::Connectivity::Face, 30, "pressure" );
  if( indices.size() > 0 )
    std::cout << mpiRank << " - " << "pressure" << " " << indices << std::endl;

  getElapsedTime( timeGetIndices );

  std::cout << mpiRank << " - numGlobalDofs - " << dofManager.numGlobalDofs() << std::endl;
  std::cout << mpiRank << " - numGlobalDofs(" "displacement" ") - " << dofManager.numGlobalDofs( "displacement" )
            << std::endl;
  std::cout << mpiRank << " - numGlobalDofs(" "pressure" ") - " << dofManager.numGlobalDofs( "pressure" ) << std::endl;
  std::cout << mpiRank << " - numLocalDofs - " << dofManager.numLocalDofs() << std::endl;
  std::cout << mpiRank << " - numLocalDofs(" "displacement" ") - " << dofManager.numLocalDofs( "displacement" )
            << std::endl;
  std::cout << mpiRank << " - numLocalDofs(" "pressure" ") - " << dofManager.numLocalDofs( "pressure" ) << std::endl;

  dofManager.printConnectivityMatrix();

  // Sum up all timings
  int mpiSize = CommunicationTools::MPI_Size( MPI_COMM_GEOSX );
  array1d<double> timesLocal( 5 ), timesSum( 5 );
  timesLocal[0] = timeAddField;
  timesLocal[1] = timeAddCoupling;
  timesLocal[2] = timeGetSingleSparsityPattern;
  timesLocal[3] = timeGetGlobalSparsityPattern;
  timesLocal[4] = timeGetIndices;

  MPI_Allreduce( timesLocal.data(), timesSum.data(), 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX );

  timeAddField = timesSum[0];
  timeAddCoupling = timesSum[1];
  timeGetSingleSparsityPattern = timesSum[2];
  timeGetGlobalSparsityPattern = timesSum[3];
  timeGetIndices = timesSum[4];

  if( mpiRank == 0 )
  {
    double totalTime = std::accumulate( timesSum.begin(), timesSum.end(), 0.0 );
    std::cout << "TIMING" << std::endl;
    std::cout << "addField: " << std::fixed << std::setprecision( 0 ) << std::trunc( timeAddField * 1e3 ) << " [ms] -- "
              << std::fixed
              << std::setprecision( 2 ) << timeAddField / totalTime * 100.0 << "%" << std::endl;
    std::cout << "addCoupling: " << std::fixed << std::setprecision( 0 ) << std::trunc( timeAddCoupling * 1e3 )
              << " [ms] -- "
              << std::fixed << std::setprecision( 2 ) << timeAddCoupling / totalTime * 100.0 << "%"
              << std::endl;
    std::cout << "setSparsityPattern(local): " << std::fixed << std::setprecision( 0 )
              << std::trunc( timeGetSingleSparsityPattern * 1e3 )
              << " [ms] -- " << std::fixed << std::setprecision( 2 )
              << timeGetSingleSparsityPattern / totalTime * 100.0
              << "%" << std::endl;
    std::cout << "setSparsityPattern(global): " << std::fixed << std::setprecision( 0 )
              << std::trunc( timeGetGlobalSparsityPattern * 1e3 )
              << " [ms] -- " << std::fixed << std::setprecision( 2 )
              << timeGetGlobalSparsityPattern / totalTime * 100.0
              << "%" << std::endl;
    std::cout << "getIndices: " << std::fixed << std::setprecision( 0 ) << std::trunc( timeGetIndices * 1e3 )
              << " [ms] -- "
              << std::fixed << std::setprecision( 2 ) << timeGetIndices / totalTime * 100.0 << "%"
              << std::endl;
    std::cout << "TOTAL: " << std::fixed << std::setprecision( 0 ) << std::trunc( totalTime * 1e3 ) << " [ms] -- "
              << std::fixed
              << std::setprecision( 2 ) << 100.0 << "%" << std::endl;
  }

  // Fake check
  SUCCEED();
}

/**
 * @function main
 * @brief Main function to setup the GEOSX environment, read the xml file and run all cases.
 */
int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  // Global call will not work because CXXUtils has already been initialized in problemManager
  //geosx::basicSetup( argc, argv );
  setupMPI(argc, argv);
  setupOpenMP();
  setupMKL();

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>( global_argc )];
  for( int i = 0 ; i < argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  // Global call will not work because CXXUtils will be destructed by problemManager
  //geosx::basicCleanup();
  FieldSpecificationManager::finalize();
  NewFunctionManager::finalize();
  finalizeMPI();

  return result;
}
