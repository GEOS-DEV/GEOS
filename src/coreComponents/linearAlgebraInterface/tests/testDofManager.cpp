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
#include "managers/initialization.hpp"
#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

#include "codingUtilities/UnitTestUtilities.hpp"

#include "DofManager.hpp"

using namespace geosx;
using namespace geosx::testing;

#ifndef GTEST_SKIP
#define GTEST_SKIP() return
#endif

#define SKIP_TEST_IF( COND, REASON ) \
do \
{ \
  if( COND ) \
  { \
    GEOS_WARNING( "This test is currently known to fail when " #COND " because:\n" REASON "\n" \
                  "Therefore, we skip it entirely for this run (may show as PASSED or SKIPPED)" ); \
    GTEST_SKIP(); \
  } \
} while(0)

#define SKIP_TEST_IN_SERIAL( REASON ) \
do \
{ \
  int const mpiSize = CommunicationTools::MPI_Size( MPI_COMM_GEOSX ); \
  SKIP_TEST_IF( mpiSize == 1, REASON ); \
} while(0)

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
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{block1}\" materialList=\"{dummy_material}\" />"
    "    <CellElementRegion name=\"region2\" cellBlocks=\"{block2}\" materialList=\"{dummy_material}\" />"
    "    <CellElementRegion name=\"region3\" cellBlocks=\"{block3}\" materialList=\"{dummy_material}\" />"
    "    <CellElementRegion name=\"region4\" cellBlocks=\"{block4}\" materialList=\"{dummy_material}\" />"
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
    dataRepository::Group * commandLine =
      problemManager->GetGroup<dataRepository::Group>( problemManager->groupKeys.commandLine );
    commandLine->registerWrapper<integer>( problemManager->viewKeys.xPartitionsOverride.Key() )->
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

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  string_array displacementRegion;
  displacementRegion.push_back( "region1" );
  displacementRegion.push_back( "region3" );
  displacementRegion.push_back( "region4" );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3,
                       displacementRegion );
  dofManager.close();

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "displacement", "displacement" );

  // Total number of nodes, sum of regions 1 and 3+4
  constexpr globalIndex nNodes = 9*5*6 + 5*5*6;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes );
  EXPECT_EQ( pattern.globalCols(), 3*nNodes );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestFEM_all.
 * @brief Define a vectorial nodal field (displacement) connected through elements (FEM) and check the
 * global number of DoFs. The field is defined on all regions.
 */
TEST_F(DofManagerTest, TestFEM_all)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3 );
  dofManager.close();

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "displacement", "displacement" );

  // Total number of nodes, sum of all regions
  constexpr globalIndex nNodes = 17*5*6;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes );
  EXPECT_EQ( pattern.globalCols(), 3*nNodes );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestFVM_partial.
 * @brief Define a scalar element field (pressure) connected through faces (FVM) and check the
 * global number of DoFs. The field is defined on 3 regions out of 4.
 */
TEST_F(DofManagerTest, TestFVM_partial)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  string_array pressureRegion;
  pressureRegion.push_back( "region1" );
  pressureRegion.push_back( "region2" );
  pressureRegion.push_back( "region3" );

  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.close();

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );

  // Total number of cells
  constexpr globalIndex nCells = 12*4*5;

  EXPECT_EQ( pattern.globalRows(), nCells );
  EXPECT_EQ( pattern.globalCols(), nCells );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestFVM_all.
 * @brief Define a scalar element field (pressure) connected through faces (FVM) and check the
 * global number of DoFs. The field is defined on all regions.
 */
TEST_F(DofManagerTest, TestFVM_all)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );
  dofManager.close();

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );

  // Total number of cells
  constexpr globalIndex nCells = 16*4*5;

  EXPECT_EQ( pattern.globalRows(), nCells );
  EXPECT_EQ( pattern.globalCols(), nCells );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestCoupling
 * @brief Define two fields, one is displacement (FEM) and the other pressure (FVM), and a coupling
 * among them. Check the sizes of the coupling block pattern.
 */
TEST_F(DofManagerTest, TestCoupling)
{
  SKIP_TEST_IN_SERIAL( "see: https://github.com/trilinos/Trilinos/issues/5663" );

  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
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
  dofManager.close();

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

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @brief Create a TPFA-type sparsity pattern by performing a fake assembly
 * @param domain the domain
 * @param mesh the mesh to use
 * @param regions list of region names to include (if empty, all regions are used)
 * @param numComp number of components per cell
 * @param sparsity the matrix to be populated
 */
void makeSparsityTPFA( DomainPartition * const domain,
                       MeshLevel * const mesh,
                       array1d<string> const & regionsInput,
                       localIndex const numComp,
                       ParallelMatrix * sparsity,
                       ParallelMatrix * connLoc )
{
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  // make a list of regions
  array1d<string> regions = regionsInput;
  if( regions.empty() )
  {
    elemManager->forElementRegions( [&]( ElementRegionBase const * const region )
    {
      regions.push_back( region->getName() );
    } );
  }

  // make an index lookup for regions (to filter when looping over faces)
  set<localIndex> regionIndices;
  for( string const & region : regions )
  {
    regionIndices.insert( elemManager->GetRegion( region )->getIndexInParent() );
  }

  // count the number of locally owned elements in active regions
  localIndex numLocalElems = 0;
  elemManager->forElementSubRegions( regions, [&]( ElementSubRegionBase * const subRegion )
  {
    numLocalElems += subRegion->GetNumberOfLocalIndices();
  } );

  // do a prefix sum to get rank offset
  globalIndex const firstLocalElem = CommunicationTools::PrefixSum<globalIndex>( numLocalElems ).first;

  // register and populate temporary DoF index arrays
  localIndex elemIndex = 0;
  elemManager->forElementSubRegions<CellElementSubRegion>( regions, [&]( CellElementSubRegion * const subRegion )
  {
    array1d<integer> const & elemGhostRank = subRegion->GhostRank();
    array1d<globalIndex> const & dofNumber =
      subRegion->template registerWrapper< array1d<globalIndex> >( "dof_index" )->reference();

    for( localIndex ei = 0; ei < subRegion->size(); ++ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        dofNumber[ei] = firstLocalElem + elemIndex++;
      }
    }
  } );

  // do the same for local/global faces
  localIndex const numLocalFaces = faceManager->GetNumberOfLocalIndices();
  globalIndex const firstLocalFace = CommunicationTools::PrefixSum<globalIndex>( numLocalFaces ).first;

  // prepare a face indexing array
  array1d<integer> const & faceGhostRank = faceManager->GhostRank();
  array1d<globalIndex> & connIndex = faceManager->registerWrapper< array1d<globalIndex> >( "conn_index" )->reference();

  localIndex faceIndex = 0;
  for( localIndex kf = 0; kf < faceManager->size(); ++kf )
  {
    if( faceGhostRank[kf] < 0 )
    {
      connIndex[kf] = firstLocalFace + faceIndex++;
    }
  }

  // sync DoF index arrays across ranks
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( "dof_index" );
  fieldNames["face"].push_back( "conn_index" );
  CommunicationTools::SynchronizeFields( fieldNames, mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // prepare data for assembly loop
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();
  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> > elemDofIndex =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex const> >( "dof_index" );

  array1d<globalIndex> localDofIndex( 2 * numComp );
  array2d<real64> localValues( 2 * numComp, 2 * numComp );
  localValues = 1.0;

  // create/resize the matrix
  if( sparsity != nullptr )
  {
    sparsity->createWithLocalSize( numLocalElems * numComp, numLocalElems * numComp, 7 * numComp, MPI_COMM_GEOSX );
  }
  if( connLoc != nullptr )
  {
    connLoc->createWithLocalSize( numLocalFaces, numLocalElems, 2, MPI_COMM_GEOSX );
  }

  // Loop over faces and assemble TPFA-style "flux" contributions
  for( localIndex kf = 0; kf < faceManager->size(); ++kf )
  {
    if( faceGhostRank[kf] >= 0 ||
        faceToElem.m_toElementRegion[kf][0] < 0 ||
        faceToElem.m_toElementRegion[kf][1] < 0 ||
        ! regionIndices.contains( faceToElem.m_toElementRegion[kf][0] ) ||
        ! regionIndices.contains( faceToElem.m_toElementRegion[kf][1] ) )
    {
      continue;
    }

    for( localIndex ke = 0; ke < 2; ++ke )
    {
      localIndex const er  = faceToElem.m_toElementRegion[kf][ke];
      localIndex const esr = faceToElem.m_toElementSubRegion[kf][ke];
      localIndex const ei  = faceToElem.m_toElementIndex[kf][ke];

      for( localIndex c = 0; c < numComp; ++c )
      {
        localDofIndex[ke * numComp + c] = elemDofIndex[er][esr][ei] * numComp + c;
      }

      if( connLoc != nullptr )
      {
        connLoc->insert( connIndex[kf], elemDofIndex[er][esr][ei], ke+1 );
      }
    }

    if( sparsity != nullptr )
    {
      sparsity->insert( localDofIndex, localDofIndex, localValues );
    }
  }

  if( sparsity != nullptr )
  {
    sparsity->close();
  }
  if( connLoc != nullptr )
  {
    connLoc->close();
  }

  // delete the temporary DoF index field
  elemManager->forElementSubRegions( regions, [&]( ElementSubRegionBase * const subRegion )
  {
    subRegion->deregisterWrapper( "dof_index" );
  } );
  faceManager->deregisterWrapper( "conn_index" );
}

/**
 * @function TestPatternTPFA
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a simple fake assembly loop
 */
TEST_F(DofManagerTest, TestPatternTPFA)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  ParallelMatrix pattern, patternExpected;

  array1d<string> regions;
  regions.push_back( "region1" );
  regions.push_back( "region3" );
  regions.push_back( "region4" );

  makeSparsityTPFA( domain, domain->getMeshBody(0)->getMeshLevel(0), regions, 3, &patternExpected, nullptr );

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, 3, regions );
  dofManager.close();
  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );

  pattern.set( 1. );
  patternExpected.set( 1. );

  EXPECT_EQ( pattern.globalRows(), patternExpected.globalRows() );
  EXPECT_EQ( pattern.globalCols(), patternExpected.globalCols() );
  EXPECT_DOUBLE_EQ( pattern.normFrobenius(), patternExpected.normFrobenius() );

  compareMatrices( pattern, patternExpected, std::numeric_limits<real64>::epsilon() );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestUserDefinedPattern
 * @brief Create an external pattern and pass it to the DofManager.
 */
TEST_F(DofManagerTest, TestUserDefinedPattern)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  array1d<string> regions;
  regions.push_back( "region1" );
  regions.push_back( "region3" );
  regions.push_back( "region4" );

  ParallelMatrix pattern, patternExpected, connLocUserInput;
  makeSparsityTPFA( domain, domain->getMeshBody(0)->getMeshLevel(0), regions, 3, &patternExpected, &connLocUserInput );

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );
  dofManager.addField( "user-defined", connLocUserInput, 3, DofManager::Connectivity::Face );
  dofManager.close();
  dofManager.setSparsityPattern( pattern, "user-defined", "user-defined" );

  pattern.set( 1. );
  patternExpected.set( 1. );

  EXPECT_EQ( pattern.globalRows(), patternExpected.globalRows() );
  EXPECT_EQ( pattern.globalCols(), patternExpected.globalCols() );
  EXPECT_DOUBLE_EQ( pattern.normFrobenius(), patternExpected.normFrobenius() );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestFEM_FVM
 * @brief Define two fields, one is displacement (FEM) and the other pressure (FVM), and check the
 * size of the global pattern.
 */
TEST_F(DofManagerTest, TestFEM_FVM)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3 );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );
  dofManager.close();

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern );

  // Total number of nodes
  constexpr globalIndex nNodes = 17*5*6;

  // Total number of cells
  constexpr globalIndex nCells = 16*4*5;

  EXPECT_EQ( pattern.globalRows(), 3*nNodes+nCells );
  EXPECT_EQ( pattern.globalCols(), 3*nNodes+nCells );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestMassMatrix
 * @brief Define a self-connected field, like a mass matrix.
 */
TEST_F(DofManagerTest, TestMassMatrix)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "massmatrix", DofManager::Location::Elem, DofManager::Connectivity::Elem );
  dofManager.close();

  ParallelMatrix pattern;
  dofManager.setSparsityPattern( pattern );

  // Total number of cells
  constexpr globalIndex nCells = 16*4*5;

  pattern.set( 1. );

  EXPECT_EQ( pattern.globalRows(), nCells );
  EXPECT_EQ( pattern.globalCols(), nCells );
  EXPECT_DOUBLE_EQ( pattern.normFrobenius(), sqrt(nCells) );

  dofManager.clear(); // remove index arrays from mesh before next test
}

/**
 * @function TestIndices
 * @brief Check the getIndices functions for FEM and FVM.
 */
TEST_F(DofManagerTest, TestIndices)
{
  int const mpiSize = CommunicationTools::MPI_Size( MPI_COMM_GEOSX );
  SKIP_TEST_IF( mpiSize > 2, "Test not designed for more than 2 ranks" );

  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3 );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face );
  dofManager.close();

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  localIndex er = 0, esr = 0, ei = 0;
  std::vector<globalIndex> indicesDispExpected;
  std::vector<globalIndex> indicesPresExpected;

  int const mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    er = 0; esr = 0; ei = 10;

    globalIndex const indicesDefaultDisp[24] = { 72, 73, 74, 75, 76, 77, 108, 109, 110, 111, 112, 113,
                                                 78, 79, 80, 81, 82, 83, 114, 115, 116, 117, 118, 119 };
    globalIndex const indicesDefaultPres[1]  = { mpiSize == 1 ? 1540 : 820 };

    indicesDispExpected.assign( indicesDefaultDisp, indicesDefaultDisp + 24 );
    indicesPresExpected.assign( indicesDefaultPres, indicesDefaultPres + 1 );
  }
  else if( mpiRank == 1 )
  {
    er = 2; esr = 0; ei = 10;

    globalIndex const indicesDefaultDisp[24] = { 756, 757, 758, 1006, 1007, 1008, 774, 775, 776, 1024, 1025, 1026,
                                                 759, 760, 761, 1009, 1010, 1011, 777, 778, 779, 1027, 1028, 1029 };
    globalIndex const indicesDefaultPres[1]  = { 1700 };

    indicesDispExpected.assign( indicesDefaultDisp, indicesDefaultDisp + 24 );
    indicesPresExpected.assign( indicesDefaultPres, indicesDefaultPres + 1 );
  }

  if( mpiRank < 2 )

  std::sort( indicesDispExpected.begin(), indicesDispExpected.end() );
  std::sort( indicesPresExpected.begin(), indicesPresExpected.end() );

  globalIndex_array indicesDisp, indicesPres;

  CellElementSubRegion const * const subregion =
    mesh->getElemManager()->GetRegion( er )->GetSubRegion<CellElementSubRegion>( esr );
  CellElementSubRegion::NodeMapType const & elemNodes = subregion->nodeList();

  arrayView1d<globalIndex const> const & dispDofIndex =
    mesh->getNodeManager()->getReference< array1d<globalIndex> >( dofManager.getKey( "displacement" ) );

  arrayView1d<globalIndex const> const & presDofIndex =
    subregion->getReference< array1d<globalIndex> >( dofManager.getKey( "pressure" ) );

  for( localIndex a = 0; a < elemNodes.size(1); ++a )
  {
    for( localIndex d = 0; d < 3; ++d )
    {
      indicesDisp.push_back( dispDofIndex[ elemNodes( ei, a ) ] + d );
    }
  }
  indicesPres.push_back( presDofIndex[ei] );

  std::sort( indicesDisp.begin(), indicesDisp.end() );
  std::sort( indicesPres.begin(), indicesPres.end() );

  EXPECT_EQ( indicesDisp.size(), indicesDispExpected.size() );
  for( localIndex i = 0; i < indicesDisp.size(); ++i )
  {
    EXPECT_EQ( indicesDisp[i], indicesDispExpected[i] );
  }

  EXPECT_EQ( indicesPres.size(), indicesPresExpected.size() );
  for( localIndex i = 0; i < indicesPres.size(); ++i )
  {
    EXPECT_EQ( indicesPres[i], indicesPresExpected[i] );
  }

  dofManager.clear(); // remove index arrays from mesh before next test
}

#define PRINT_PATTERNS 0

/**
 * @function TestWithTimes
 * @brief This is not a real test (it is always true!). It collects timings.
 */
TEST_F(DofManagerTest, TestWithTimes)
{
  SKIP_TEST_IN_SERIAL( "see: https://github.com/trilinos/Trilinos/issues/5663" );

  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager( "test" );
  dofManager.setMesh( domain, 0, 0 );

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

  ParallelMatrix connLocUserInput;
  makeSparsityTPFA( domain, domain->getMeshBody(0)->getMeshLevel(0), pressureRegion, 3, nullptr, &connLocUserInput );

  double timeAddField, timeAddCoupling, timeGetSingleSparsityPattern, timeGetGlobalSparsityPattern;

  setTimer( timeAddField );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 3,
                       displacementRegion );
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.addField( "massmatrix", DofManager::Location::Elem, DofManager::Connectivity::None, 2, testRegion3 );
  dofManager.addField( "user-defined", connLocUserInput, 3, DofManager::Connectivity::Face );

  getElapsedTime( timeAddField );
  setTimer( timeAddCoupling );

  dofManager.addCoupling( "displacement", "pressure", DofManager::Connectivity::Elem, testRegion3, false );
  dofManager.addCoupling( "massmatrix", "pressure", DofManager::Connectivity::Elem );

  getElapsedTime( timeAddCoupling );
  setTimer( timeGetSingleSparsityPattern );

  dofManager.close();

  ParallelMatrix pattern;

  dofManager.setSparsityPattern( pattern, "displacement", "displacement" );
#if PRINT_PATTERNS
    pattern.write( "displacement.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "pressure", "pressure" );
#if PRINT_PATTERNS
    pattern.write( "pressure.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "massmatrix", "massmatrix" );
#if PRINT_PATTERNS
    pattern.write( "massmatrix.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "displacement", "pressure" );
#if PRINT_PATTERNS
    pattern.write( "coupling1.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "pressure", "displacement" );
#if PRINT_PATTERNS
    pattern.write( "coupling1_empty.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "pressure", "massmatrix" );
#if PRINT_PATTERNS
    pattern.write( "coupling2.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "massmatrix", "pressure" );
#if PRINT_PATTERNS
    pattern.write( "coupling2_transp.mtx" );
#endif

  dofManager.setSparsityPattern( pattern, "user-defined", "user-defined" );
#if PRINT_PATTERNS
    pattern.write( "user-defined.mtx" );
#endif

  getElapsedTime( timeGetSingleSparsityPattern );
  setTimer( timeGetGlobalSparsityPattern );

  dofManager.setSparsityPattern( pattern );
#if PRINT_PATTERNS
    pattern.write( "global.mtx" );
#endif

  getElapsedTime( timeGetGlobalSparsityPattern );

  int mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );

  GEOS_LOG_RANK( "numGlobalDofs = " << dofManager.numGlobalDofs() );
  GEOS_LOG_RANK( "numGlobalDofs(displacement) = " << dofManager.numGlobalDofs( "displacement" ) );
  GEOS_LOG_RANK( "numGlobalDofs(pressure) = " << dofManager.numGlobalDofs( "pressure" ) );
  GEOS_LOG_RANK( "numLocalDofs = " << dofManager.numLocalDofs() );
  GEOS_LOG_RANK( "numLocalDofs(displacement) = " << dofManager.numLocalDofs( "displacement" ) );
  GEOS_LOG_RANK( "numLocalDofs(pressure) = " << dofManager.numLocalDofs( "pressure" ) );

  dofManager.printConnectivityMatrix();

  // Sum up all timings
  array1d<double> timesLocal( 5 ), timesSum( 5 );
  timesLocal[0] = timeAddField;
  timesLocal[1] = timeAddCoupling;
  timesLocal[2] = timeGetSingleSparsityPattern;
  timesLocal[3] = timeGetGlobalSparsityPattern;

  MPI_Allreduce( timesLocal.data(), timesSum.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX );

  timeAddField = timesSum[0];
  timeAddCoupling = timesSum[1];
  timeGetSingleSparsityPattern = timesSum[2];
  timeGetGlobalSparsityPattern = timesSum[3];

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
    std::cout << "TOTAL: " << std::fixed << std::setprecision( 0 ) << std::trunc( totalTime * 1e3 ) << " [ms] -- "
              << std::fixed
              << std::setprecision( 2 ) << 100.0 << "%" << std::endl;
  }

  dofManager.clear(); // remove index arrays from mesh before next test

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
  geosx::basicSetup( argc, argv );

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>( global_argc )];
  for( int i = 0 ; i < argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  // Global call will not work because CXXUtils will be destructed by problemManager
  geosx::basicCleanup();

  return result;
}
