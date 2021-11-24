/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testLAIHelperFunctions.cpp
 */

#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mainInterface/GeosxState.hpp"
#include "unitTests/linearAlgebraTests/testDofManagerUtils.hpp"

#include <gtest/gtest.h>

using namespace geosx;

static real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();
static real64 const tolerance  = machinePrecision;//1e-10;

char const * xmlInput =
  "<Problem>"
  "  <Mesh>"
  "    <InternalMesh name=\"mesh1\""
  "                  elementTypes=\"{C3D8}\""
  "                  xCoords=\"{0, 1}\""
  "                  yCoords=\"{0, 1}\""
  "                  zCoords=\"{0, 1}\""
  "                  nx=\"{6}\""
  "                  ny=\"{9}\""
  "                  nz=\"{5}\""
  "                  cellBlockNames=\"{block1}\"/>"
  "  </Mesh>"
  "  <ElementRegions>"
  "    <CellElementRegion name=\"region1\" cellBlocks=\"{block1}\" materialList=\"{dummy_material}\" />"
  "  </ElementRegions>"
  "</Problem>";

template< typename LAI >
class LAIHelperFunctionsTest : public ::testing::Test
{
protected:

  using Base = ::testing::Test;

  LAIHelperFunctionsTest():
    Base(),
    state( std::make_unique< CommandLineOptions >() )
  {
    geosx::testing::setupProblemFromXML( &state.getProblemManager(), xmlInput );
    mesh = &state.getProblemManager().getDomainPartition().getMeshBody( 0 ).getMeshLevel( 0 );
  }

  GeosxState state;
  MeshLevel * mesh;
};

TYPED_TEST_SUITE_P( LAIHelperFunctionsTest );

TYPED_TEST_P( LAIHelperFunctionsTest, nodalVectorPermutation )
{
  using Matrix = typename TypeParam::ParallelMatrix;
  using Vector = typename TypeParam::ParallelVector;

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = meshLevel.getNodeManager();

  arrayView1d< globalIndex const > const nodeLocalToGlobal = nodeManager.localToGlobalMap();

  DofManager dofManager( "test" );
  dofManager.setDomain( meshLevel );

  string_array regions;
  regions.emplace_back( "region1" );

  dofManager.addField( "nodalVariable", DofManager::Location::Node, 3, regions );
  dofManager.addCoupling( "nodalVariable", "nodalVariable", DofManager::Connector::Elem );
  dofManager.reorderByRank();

  localIndex const numLocalDof = 3 * nodeManager.getNumberOfLocalIndices();

  arrayView1d< globalIndex const > const dofNumber = nodeManager.getReference< globalIndex_array >( dofManager.getKey( "nodalVariable" ) );
  arrayView1d< integer const > const isNodeGhost = nodeManager.ghostRank();

  Vector nodalVariable, expectedPermutedVector;
  nodalVariable.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  nodalVariable.set( 0 );
  expectedPermutedVector.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  expectedPermutedVector.set( 0 );

  nodalVariable.open();
  expectedPermutedVector.open();
  for( localIndex a = 0; a < nodeManager.size(); ++a )
  {
    if( isNodeGhost[a] < 0 )
    {
      for( localIndex d = 0; d < 3; ++d )
      {
        real64 const value = nodeLocalToGlobal[a] * 3 + d;
        nodalVariable.add( dofNumber[a] + d, value );
        expectedPermutedVector.add( nodeLocalToGlobal[a] * 3 + d, value );
      }
    }
  }
  nodalVariable.close();
  expectedPermutedVector.close();

  Matrix permutationMatrix;
  LAIHelperFunctions::createPermutationMatrix( nodeManager,
                                               3,
                                               dofManager.getKey( "nodalVariable" ),
                                               permutationMatrix );
  Vector permutedVector = LAIHelperFunctions::permuteVector( nodalVariable, permutationMatrix );

  permutedVector.axpy( -1, expectedPermutedVector );
  real64 const vectorNorm = permutedVector.norm1();
  EXPECT_LT( vectorNorm, tolerance );
}

TYPED_TEST_P( LAIHelperFunctionsTest, cellCenteredVectorPermutation )
{
  using Matrix = typename TypeParam::ParallelMatrix;
  using Vector = typename TypeParam::ParallelVector;

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  DofManager dofManager( "test" );
  dofManager.setDomain( meshLevel );

  string_array regions;
  regions.emplace_back( "region1" );

  dofManager.addField( "cellCentered", DofManager::Location::Elem, 1, regions );
  dofManager.addCoupling( "cellCentered", "cellCentered", DofManager::Connector::Face );
  dofManager.reorderByRank();

  integer const numLocalDof = dofManager.numLocalDofs( "cellCentered" );

  Vector cellCenteredVariable, expectedPermutedVector;
  cellCenteredVariable.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  cellCenteredVariable.set( 0 );
  expectedPermutedVector.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  expectedPermutedVector.set( 0 );

  cellCenteredVariable.open();
  expectedPermutedVector.open();
  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
    arrayView1d< globalIndex const > const dofNumber = elementSubRegion.getReference< array1d< globalIndex > >( dofManager.getKey( "cellCentered" ) );
    arrayView1d< integer const > const isGhost = elementSubRegion.ghostRank();
    arrayView1d< globalIndex const > const localToGlobal = elementSubRegion.localToGlobalMap();

    for( localIndex k = 0; k < numElems; ++k )
    {
      if( dofNumber[k] >= 0 && isGhost[k] < 0 )
      {
        real64 const value = localToGlobal[k];
        cellCenteredVariable.add( dofNumber[k], value );
        expectedPermutedVector.add( localToGlobal[k], value );
      }
    }
  } );

  cellCenteredVariable.close();
  expectedPermutedVector.close();

  Matrix permutationMatrix;
  LAIHelperFunctions::createPermutationMatrix( elemManager,
                                               1,
                                               dofManager.getKey( "cellCentered" ),
                                               permutationMatrix );
  Vector permutedVector = LAIHelperFunctions::permuteVector( cellCenteredVariable, permutationMatrix );

  permutedVector.axpy( -1, expectedPermutedVector );
  real64 const vectorNorm = permutedVector.norm1();
  EXPECT_LT( vectorNorm, tolerance );
}

REGISTER_TYPED_TEST_SUITE_P( LAIHelperFunctionsTest,
                             nodalVectorPermutation,
                             cellCenteredVectorPermutation );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, LAIHelperFunctionsTest, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, LAIHelperFunctionsTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, LAIHelperFunctionsTest, PetscInterface, );
#endif

/**
 * @function main
 * @brief Main function to setup the GEOSX environment, read the xml file and run all cases.
 */
int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
