/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testDofManagerUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP
#define GEOSX_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP

#include <gtest/gtest.h>

namespace geosx
{

namespace testing
{

void setupProblem( ProblemManager * const problemManager, char const * const xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( xmlInput, strlen( xmlInput ) );
  if (!xmlResult)
  {
    GEOSX_LOG_RANK_0("XML parsed with errors!");
    GEOSX_LOG_RANK_0("Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0("Error offset: " << xmlResult.offset);
  }

  int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
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

string_array getRegions( MeshLevel const * const mesh, std::vector<string> const & input )
{
  string_array regions;
  if( !input.empty() )
  {
    regions.insert( 0, input.data(), input.size() );
  }
  else
  {
    mesh->getElemManager()->forElementRegions( [&]( ElementRegionBase const * const region )
    {
      regions.push_back( region->getName() );
    } );
  }
  return regions;
}

template<DofManager::Location LOC>
struct testMeshHelper {};

template<>
struct testMeshHelper<DofManager::Location::Node>
{
  static auto constexpr managerKey = MeshLevel::groupStructKeys::nodeManagerString;
  static auto constexpr elemMapKey = ElementSubRegionBase::viewKeyStruct::nodeListString;
  template<typename SUBREGION> using ElemToObjMap = typename SUBREGION::NodeMapType;
};

template<>
struct testMeshHelper<DofManager::Location::Face>
{
  static auto constexpr managerKey = MeshLevel::groupStructKeys::faceManagerString;
  static auto constexpr elemMapKey = ElementSubRegionBase::viewKeyStruct::faceListString;
  template<typename SUBREGION> using ElemToObjMap = typename SUBREGION::FaceMapType;
};

template<typename PERM>
localIndex size1( InterObjectRelation<array2d<localIndex, PERM>> const & map, localIndex const GEOSX_UNUSED_ARG( i0 ) )
{
  return map.size( 1 );
}

localIndex size1( InterObjectRelation<array1d<array1d<localIndex>>> const & map, localIndex const i0 )
{
return map[i0].size();
}

localIndex size1( InterObjectRelation<ArrayOfArrays<localIndex>> const & map, localIndex const i0 )
{
  return map.sizeOfArray( i0 );
}

template<DofManager::Location LOC>
struct forLocalObjectsImpl
{
  template<typename LAMBDA>
  static void f( MeshLevel const * const mesh,
                 string_array const & regions,
                 LAMBDA lambda )
  {
    using helper = testMeshHelper<LOC>;
    ObjectManagerBase const * const manager = mesh->GetGroup<ObjectManagerBase>( helper::managerKey );

    arrayView1d<integer const> ghostRank =
      manager->getReference< array1d<integer> >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    array1d<bool> visited( ghostRank.size() );
    visited = false;

    mesh->getElemManager()->forElementSubRegions( regions, [&]( auto const * const subRegion )
    {
      using MapType = typename helper::template ElemToObjMap<std::remove_pointer_t<decltype( subRegion )>>;
      MapType const & elemToObjMap =
        subRegion->template getReference<MapType>( helper::elemMapKey );

      for( localIndex k = 0; k < subRegion->size(); ++k )
      {
        for( localIndex a = 0; a < size1( elemToObjMap, k ); ++a )
        {
          localIndex const idx = elemToObjMap[k][a];
          if( ghostRank[idx] < 0 && !std::exchange( visited[idx], true ) )
          {
            lambda( idx );
          }
        }
      }
    } );
  }
};

template<>
struct forLocalObjectsImpl<DofManager::Location::Elem>
{
  template<typename LAMBDA>
  static void f( MeshLevel const * const mesh,
                 string_array const & regions,
                 LAMBDA lambda )
  {
    // make a list of regions
    ElementRegionManager const * const elemManager = mesh->getElemManager();

    elemManager->forElementSubRegionsComplete( regions, [&]( localIndex const er,
                                                             localIndex const esr,
                                                             ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                                             ElementSubRegionBase const * const subRegion )
    {
      arrayView1d<integer const> ghostRank =
        subRegion->getReference< array1d<integer> > ( ObjectManagerBase::viewKeyStruct::ghostRankString );

      for( localIndex ei = 0; ei < subRegion->size(); ++ei )
      {
        if( ghostRank[ei] < 0 )
        {
          lambda( std::array<localIndex, 3> { er, esr, ei} );
        }
      }
    } );
  }
};

template<DofManager::Location LOC, typename LAMBDA>
void forLocalObjects( MeshLevel const * const mesh,
                      array1d<string> const & regions,
                      LAMBDA && lambda )
{;
  forLocalObjectsImpl<LOC>::template f( mesh, regions, std::forward<LAMBDA>( lambda ) );
}

template<DofManager::Location LOC>
localIndex countLocalObjects( MeshLevel const * const mesh, array1d<string> const & regions )
{
  localIndex numLocal = 0;
  forLocalObjects<LOC>( mesh, regions, [&]( auto const & ) { ++numLocal; } );
  return numLocal;
}

/**
 * @brief Create a TPFA-type sparsity pattern.
 * @param domain the domain
 * @param mesh the mesh to use
 * @param regionsInput list of region names to include (if empty, all regions are used)
 * @param numComp number of components per cell
 * @param sparsity the matrix to be populated, must be properly sized.
 */
void makeSparsityTPFA( MeshLevel const * const mesh,
                       string const & dofIndexKey,
                       string_array const & regions,
                       localIndex const numComp,
                       ParallelMatrix & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> > elemDofIndex =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex const> >( dofIndexKey );

  // Make a set of target region indices to check face fluxes.
  set<localIndex> regionSet;
  elemManager->forElementRegions( regions, [&]( ElementRegionBase const * const region )
  {
    regionSet.insert( region->getIndexInParent() );
  } );

  // prepare data for assembly loop
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();
  localIndex const numElem = faceToElem.size( 1 );

  array1d<globalIndex> localDofIndex( numElem * numComp );
  array2d<real64> localValues( numElem * numComp, numElem * numComp );
  localValues = 1.0;

  // Loop over faces and assemble TPFA-style "flux" contributions
  forLocalObjects<DofManager::Location::Face>( mesh, regions, [&]( localIndex const kf )
  {
    localIndex const er0 = faceToElem.m_toElementRegion[kf][0];
    localIndex const er1 = faceToElem.m_toElementRegion[kf][1];

    if( er0 >= 0 && er1 >= 0 && regionSet.contains( er0 ) && regionSet.contains( er1 ) )
    {
      for( localIndex ke = 0; ke < numElem; ++ke )
      {
        localIndex const er  = faceToElem.m_toElementRegion[kf][ke];
        localIndex const esr = faceToElem.m_toElementSubRegion[kf][ke];
        localIndex const ei  = faceToElem.m_toElementIndex[kf][ke];

        for( localIndex c = 0; c < numComp; ++c )
        {
          localDofIndex[ke * numComp + c] = elemDofIndex[er][esr][ei] + c;
        }
      }
      sparsity.insert( localDofIndex, localDofIndex, localValues );
    }
  } );
}

/**
 * @brief Populate a FEM-type sparsity pattern.
 * @param domain the domain
 * @param mesh the mesh to use
 * @param regionsInput list of region names to include (if empty, all regions are used)
 * @param numComp number of components per cell
 * @param sparsity the matrix to be populated, must be properly sized.
 */
void makeSparsityFEM( MeshLevel const * const mesh,
                      string const & dofIndexKey,
                      string_array const & regions,
                      localIndex const numComp,
                      ParallelMatrix & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();

  arrayView1d<globalIndex const> nodeDofIndex =
    nodeManager->getReference< array1d<globalIndex> >( dofIndexKey );

  // perform assembly loop over elements
  elemManager->forElementSubRegions( regions, [&]( auto const * const subRegion )
  {
    using NodeMapType = TYPEOFPTR( subRegion )::NodeMapType;
    NodeMapType const & nodeMap =
      subRegion->template getReference<NodeMapType>( ElementSubRegionBase::viewKeyStruct::nodeListString );

    localIndex const numNode = subRegion->numNodesPerElement();
    array1d<globalIndex> localDofIndex( numNode * numComp );
    array2d<real64> localValues( numNode * numComp, numNode * numComp );
    localValues = 1.0;

    for( localIndex k = 0; k < subRegion->size(); ++k )
    {
      for( localIndex a = 0; a < numNode; ++a )
      {
        for( localIndex c = 0; c < numComp; ++c )
        {
          localDofIndex[a * numComp + c] = nodeDofIndex[nodeMap[k][a]] + c;
        }
      }

      sparsity.insert( localDofIndex, localDofIndex, localValues );
    }
  } );
}

/**
 * @brief Populate a FEM/FVM coupling sparsity.
 * @param domain the domain
 * @param mesh the mesh to use
 * @param regionsInput list of region names to include (if empty, all regions are used)
 * @param numCompNode number of components per cell
 * @param sparsity the matrix to be populated, must be properly sized.
 */
void makeSparsityFEM_FVM( MeshLevel const * const mesh,
                          string const & dofIndexKeyNode,
                          string const & dofIndexKeyElem,
                          string_array const & regions,
                          localIndex const numCompNode,
                          localIndex const numCompElem,
                          ParallelMatrix & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();

  arrayView1d<globalIndex const> nodeDofIndex =
    nodeManager->getReference< array1d<globalIndex> >( dofIndexKeyNode );

  // perform assembly loop over elements
  elemManager->forElementSubRegions( regions, [&]( auto const * const subRegion )
  {
    using NodeMapType = TYPEOFPTR( subRegion )::NodeMapType;
    NodeMapType const & nodeMap =
      subRegion->template getReference<NodeMapType>( ElementSubRegionBase::viewKeyStruct::nodeListString );

    arrayView1d<globalIndex const> elemDofIndex =
      subRegion->template getReference< array1d<globalIndex> >( dofIndexKeyElem );

    localIndex const numNode = subRegion->numNodesPerElement();

    array1d<globalIndex> localNodeDofIndex( numNode * numCompNode );
    array1d<globalIndex> localElemDofIndex( numCompElem );
    array2d<real64> localValues1( numNode * numCompNode, numCompElem );
    localValues1 = 1.0;
    array2d<real64> localValues2( numCompElem, numNode * numCompNode );
    localValues2 = 1.0;

    for( localIndex k = 0; k < subRegion->size(); ++k )
    {
      for( localIndex c = 0; c < numCompElem; ++c )
      {
        localElemDofIndex[c] = elemDofIndex[k] + c;
      }
      for( localIndex a = 0; a < numNode; ++a )
      {
        for( localIndex c = 0; c < numCompNode; ++c )
        {
          localNodeDofIndex[a * numCompNode + c] = nodeDofIndex[nodeMap[k][a]] + c;
        }
      }

      sparsity.insert( localNodeDofIndex, localElemDofIndex, localValues1 );
      sparsity.insert( localElemDofIndex, localNodeDofIndex, localValues2 );
    }
  } );
}

/**
 * @brief Create a mass matrix-type sparsity pattern (diagonal)
 * @param domain the domain
 * @param mesh the mesh to use
 * @param regionsInput list of region names to include (if empty, all regions are used)
 * @param numComp number of components per cell
 * @param sparsity the matrix to be populated
 */
void makeSparsityMass( MeshLevel const * const mesh,
                       string const & dofIndexKey,
                       string_array const & regions,
                       localIndex const numComp,
                       ParallelMatrix & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> > elemDofIndex =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex const> >( dofIndexKey );

  array1d<globalIndex> localDofIndex( numComp );
  array2d<real64> localValues( numComp, numComp );
  localValues = 1.0;

  forLocalObjects<DofManager::Location::Elem>( mesh, regions, [&]( auto const & idx)
  {
    for( localIndex c = 0; c < numComp; ++c )
    {
      localDofIndex[c] = elemDofIndex[idx[0]][idx[1]][idx[2]] + c;
    }
    sparsity.insert( localDofIndex, localDofIndex, localValues );
  } );
}

/**
 * @brief Create a flux sparsity pattern (face-based dofs coupled via elements)
 * @param domain the domain
 * @param mesh the mesh to use
 * @param regionsInput list of region names to include (if empty, all regions are used)
 * @param numComp number of components per cell
 * @param sparsity the matrix to be populated
 */
void makeSparsityFlux( MeshLevel const * const mesh,
                       string const & dofIndexKey,
                       string_array const & regions,
                       localIndex const numComp,
                       ParallelMatrix & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  arrayView1d<globalIndex const> faceDofIndex =
    faceManager->getReference< array1d<globalIndex> > ( dofIndexKey );

  // perform assembly loop over elements
  elemManager->forElementSubRegions( regions, [&]( auto const * const subRegion )
  {
    using FaceMapType = TYPEOFPTR( subRegion )::FaceMapType;
    FaceMapType const & faceMap =
      subRegion->template getReference<FaceMapType>( ElementSubRegionBase::viewKeyStruct::faceListString );

    localIndex const numFace = subRegion->numFacesPerElement();
    array1d<globalIndex> localDofIndex( numFace * numComp );
    array2d<real64> localValues( numFace * numComp, numFace * numComp );
    localValues = 1.0;

    for( localIndex k = 0; k < subRegion->size(); ++k )
    {
      for( localIndex a = 0; a < numFace; ++a )
      {
        for( localIndex c = 0; c < numComp; ++c )
        {
          localDofIndex[a * numComp + c] = faceDofIndex[faceMap[k][a]] + c;
        }
      }

      sparsity.insert( localDofIndex, localDofIndex, localValues );
    }
  } );
}

}

}

#endif //GEOSX_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP
