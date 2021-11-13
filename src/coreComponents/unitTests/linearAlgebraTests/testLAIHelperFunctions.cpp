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
  "    <CellElementRegion name=\"region1\" cellBlocks=\"{block1}\" materialList=\"{dummy}\" />"
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

  string const fieldName = "nodalVariable";
  integer constexpr numComp = 3;

  DofManager dofManager( "test" );
  dofManager.setMesh( meshLevel );
  dofManager.addField( fieldName, DofManager::Location::Node, numComp );
  dofManager.addCoupling( fieldName, fieldName, DofManager::Connector::Elem );
  dofManager.reorderByRank();

  arrayView1d< globalIndex const > const dofNumber = nodeManager.getReference< globalIndex_array >( dofManager.getKey( fieldName ) );
  arrayView1d< integer const > const isNodeGhost = nodeManager.ghostRank();

  Vector nodalVariable;
  nodalVariable.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< real64 > const nodalVariableView = nodalVariable.open();
  forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    if( dofNumber[a] >= 0 && isNodeGhost[a] < 0 )
    {
      for( localIndex c = 0; c < numComp; ++c )
      {
        nodalVariableView[dofNumber[a] + c] = static_cast< real64 >( nodeLocalToGlobal[a] * numComp + c );
      }
    }
  } );
  nodalVariable.close();

  Matrix permutationMatrix;
  LAIHelperFunctions::createPermutationMatrix( nodeManager,
                                               numComp,
                                               dofManager.getKey( fieldName ),
                                               permutationMatrix );
  Vector permutedVector = LAIHelperFunctions::permuteVector( nodalVariable, permutationMatrix );

  // After permutation, the vector should become an increasing sequence of numbers
  arrayView1d< real64 const > permutedValues = permutedVector.values();
  forAll< serialPolicy >( permutedValues.size(), [=]( localIndex const i )
  {
    EXPECT_EQ( permutedValues[i], rankOffset + i );
  } );
}

TYPED_TEST_P( LAIHelperFunctionsTest, cellCenteredVectorPermutation )
{
  using Matrix = typename TypeParam::ParallelMatrix;
  using Vector = typename TypeParam::ParallelVector;

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  string const fieldName = "cellCenteredVariable";
  integer constexpr numComp = 3;

  DofManager dofManager( "test" );
  dofManager.setMesh( meshLevel );
  dofManager.addField( fieldName, DofManager::Location::Elem, numComp );
  dofManager.addCoupling( fieldName, fieldName, DofManager::Connector::Face );
  dofManager.reorderByRank();

  Vector cellCenteredVariable;
  cellCenteredVariable.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< real64 > const cellCenteredVariableView = cellCenteredVariable.open();
  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & elementSubRegion )
  {
    arrayView1d< globalIndex const > const dofNumber = elementSubRegion.getReference< array1d< globalIndex > >( dofManager.getKey( fieldName ) );
    arrayView1d< integer const > const isGhost = elementSubRegion.ghostRank();
    arrayView1d< globalIndex const > const localToGlobal = elementSubRegion.localToGlobalMap();

    forAll< parallelDevicePolicy<> >( elementSubRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      if( dofNumber[k] >= 0 && isGhost[k] < 0 )
      {
        for( localIndex c = 0; c < numComp; ++c )
        {
          cellCenteredVariableView[dofNumber[k] - rankOffset + c] = static_cast< real64 >( localToGlobal[k] * numComp + c );
        }
      }
    } );
  } );
  cellCenteredVariable.close();

  Matrix permutationMatrix;
  LAIHelperFunctions::createPermutationMatrix( elemManager,
                                               numComp,
                                               dofManager.getKey( fieldName ),
                                               permutationMatrix );
  Vector permutedVector = LAIHelperFunctions::permuteVector( cellCenteredVariable, permutationMatrix );

  // After permutation, the vector should become an increasing sequence of numbers
  arrayView1d< real64 const > permutedValues = permutedVector.values();
  forAll< serialPolicy >( permutedValues.size(), [=]( localIndex const i )
  {
    EXPECT_EQ( permutedValues[i], rankOffset + i );
  } );
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
