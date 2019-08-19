/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file FaceManager.cpp
 */

#include "FaceManager.hpp"
#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

/**
 *
 * @return
 */
FaceManager::FaceManager( string const &, ManagedGroup * const parent ):
  ObjectManagerBase("FaceManager",parent)
{
  this->RegisterViewWrapper( viewKeyStruct::nodeListString, &m_nodeList, false );
  this->RegisterViewWrapper( viewKeyStruct::edgeListString, &m_edgeList, false );
//  m_nodeList.SetRelatedObject( parent->getGroup<NodeManager>(MeshLevel::groupStructKeys::nodeManagerString));

  this->RegisterViewWrapper( viewKeyStruct::elementRegionListString,
                             &(m_toElements.m_toElementRegion),
                             false )->
    setApplyDefaultValue(-1);

  this->RegisterViewWrapper( viewKeyStruct::elementSubRegionListString,
                             &(m_toElements.m_toElementSubRegion),
                             false )->
    setApplyDefaultValue(-1);

  this->RegisterViewWrapper( viewKeyStruct::elementListString,
                             &(m_toElements.m_toElementIndex),
                             false )->
    setApplyDefaultValue(-1);

  this->RegisterViewWrapper( viewKeyStruct::faceAreaString, &m_faceArea, false);
  this->RegisterViewWrapper( viewKeyStruct::faceCenterString, &m_faceCenter, false);
  this->RegisterViewWrapper( viewKeyStruct::faceNormalString, &m_faceNormal, false);

  m_toElements.resize(0,2);

  //0-based; note that the following field is ALSO 0
  //for faces that are not external faces, so check isExternal before using
//  this->AddKeylessDataField<localIndex>("externalFaceIndex", true, true);
//
//  this->AddKeylessDataField<R1Tensor>("FaceCenter",true,true);
}

/**
 *
 * @return
 */
FaceManager::~FaceManager()
{}

/**
 * @class FaceBuilder
 * @brief This class stores the data necessary to construct the various face maps.
 */
struct FaceBuilder
{
  /**
   * @brief Constructor.
   * @param [in] n1_ the second smallest node index that comprise the face.
   * @param [in] n2_ the third smallest node index that comprise the face.
   * @param [in] er_ the element region this face came from.
   * @param [in] esr_ the element sub region this face came from.
   * @param [in] k_ the element this face came from.
   * @param [in] elementLocalFaceIndex_ the element local index of this face.
   */
  FaceBuilder( localIndex const n1_,
               localIndex const n2_,
               localIndex const er_,
               localIndex const esr_,
               localIndex const k_,
               localIndex const elementLocalFaceIndex_ ) :
    n1( uint32_t( n1_ ) ),
    n2( uint32_t( n2_ ) ),
    er( uint32_t( er_ ) ),
    esr( uint32_t( esr_ ) ),
    k( uint32_t( k_ ) ),
    elementLocalFaceIndex( uint32_t( elementLocalFaceIndex_ ) )
  {}

  /**
   * @brief Imposes an ordering on FaceBuilders. First compares n1, then n2, then er, esr, and k.
   * @param [in] rhs the FaceBuilder to compare against.
   */
  bool operator<( FaceBuilder const & rhs ) const
  {
    if ( n1 < rhs.n1 ) return true;
    if ( n1 > rhs.n1 ) return false;
    if ( n2 < rhs.n2 ) return true;
    if ( n2 > rhs.n2 ) return false;
    if ( er < rhs.er ) return true;
    if ( er > rhs.er ) return false;
    if ( esr < rhs.esr ) return true;
    if ( esr > rhs.esr ) return false;
    return k < rhs.k;
  }

  /**
   * @brief Return true if the two FaceBuilders share the same second and third smallest nodes.
   * @param [in] rhs the FaceBuilder to compare against.
   */
  bool operator==( FaceBuilder const & rhs ) const
  { return n1 == rhs.n1 && n2 == rhs.n2; }

  uint32_t n1;
  uint32_t n2;
  uint32_t er;
  uint32_t esr;
  uint32_t k;
  uint32_t elementLocalFaceIndex;
};

/**
 * @brief Get the three smallest values.
 * @param [in] values the array to search.
 * @param [out] minValues an array that will hold the three smallest values, from least to greatest.
 */
void findSmallestThreeValues( arrayView1d< localIndex const > const & values, localIndex (&minValues)[3] )
{
  localIndex const n = values.size();
  GEOS_ASSERT_GE( n, 3 );

  // Pick out the first three values
  minValues[0] = values[0];
  minValues[1] = values[1];
  minValues[2] = values[2];

  // and sort them
  if (minValues[0] > minValues[1]) std::swap( minValues[0], minValues[1] );
  if (minValues[1] > minValues[2]) std::swap( minValues[1], minValues[2] );
  if (minValues[0] > minValues[1]) std::swap( minValues[0], minValues[1] );

  for( localIndex i = 3; i < n; ++i )
  {
    // If the current value is greater than the last minimum values skip it.
    if (values[i] > minValues[2]) continue;

    minValues[2] = values[i];
    if (minValues[1] > minValues[2]) std::swap( minValues[1], minValues[2] );
    if (minValues[0] > minValues[1]) std::swap( minValues[0], minValues[1] );
  }
}

/**
 * @brief Populate the facesByLowestNode map. Return the maximum number of nodes associated with any face.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in/out] facesByLowestNode of size numNodes, where each sub array has been preallocated to hold
 *        *enough* space.
 * For each face of each element, this function gets the three lowest nodes in the face {n0, n1, n2}, creates
 * an EdgeBuilder associated with the face from n1 and n2 and then appends the EdgeBuilder to facesByLowestNode[ n0 ].
 * Finally it sorts the contents of each sub-array of facesByLowestNode from least to greatest.
 */
localIndex createFacesByLowestNode( ElementRegionManager const & elementManager,
                                    ArrayOfArraysView< FaceBuilder > const & facesByLowestNode )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex maxFaceNodes = 0;

  // loop over all the regions
  for( typename dataRepository::indexType er = 0; er < elementManager.numRegions(); ++er )
  {
    ElementRegion const & elemRegion = *elementManager.GetRegion(er);

    // loop over all the subregions
    elemRegion.forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr,
                                                                    CellElementSubRegion const * const subRegion )
    {
      localIndex const numFacesPerElement = subRegion->numFacesPerElement();
      localIndex const numElements = subRegion->size();
      array2d<localIndex const> const & elemsToFaces = subRegion->faceList();

      maxFaceNodes = std::max( maxFaceNodes, subRegion->GetMaxNumFaceNodes() );
      
      // Begin the parallel region so that tempNodeList and lowestNodes are thread private.
      PRAGMA_OMP( omp parallel )
      {
        localIndex_array tempNodeList;
        localIndex lowestNodes[3];

        // Loop over all the elements.
        PRAGMA_OMP( omp for )
        for( localIndex k = 0; k < numElements; ++k )
        {
          for( localIndex elementLocalFaceIndex = 0; elementLocalFaceIndex < numFacesPerElement; ++elementLocalFaceIndex )
          {
            subRegion->GetFaceNodes( k, elementLocalFaceIndex, tempNodeList );
            findSmallestThreeValues( tempNodeList, lowestNodes );

            facesByLowestNode.atomicAppendToArray( RAJA::atomic::auto_atomic{}, lowestNodes[0], FaceBuilder( lowestNodes[1], lowestNodes[2], er, esr, k, elementLocalFaceIndex ) );
          }
        }
      }
    });
  }

  // Loop over all the nodes and sort the associated faces.
  forall_in_range< parallelHostPolicy >( 0, numNodes, [&]( localIndex const nodeID )
  {
    FaceBuilder * const faces = facesByLowestNode[ nodeID ];
    std::sort( faces, faces + facesByLowestNode.sizeOfArray( nodeID ) );
  } );

  return maxFaceNodes;
}

/**
 * @brief Return the total number of unique faces and fill in the uniqueFaceOffsets array.
 * @param [in] facesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [out] uniqueFaceOffsets an array of size numNodes + 1. After this function returns node i contains
 * faces with IDs ranging from uniqueFaceOffsets[ i ] to uniqueFaceOffsets[ i + 1 ] - 1.
 */
localIndex calculateTotalNumberOfFaces( ArrayOfArraysView< FaceBuilder const > const & facesByLowestNode,
                                        arrayView1d< localIndex > const & uniqueFaceOffsets )
{
  localIndex const numNodes = facesByLowestNode.size();
  GEOS_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );

  uniqueFaceOffsets[0] = 0;

  // Loop over all the nodes.
  forall_in_range< parallelHostPolicy >( 0, numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );
    
    // If there are no faces associated with this node we can skip it.
    if ( numFaces == 0 ) return;

    localIndex & numUniqueFaces = uniqueFaceOffsets[ nodeID + 1 ];
    numUniqueFaces = 0;

    // Otherwise since facesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries. Since each face can appear at most twice if we find a match we
    // can skip the next entry as well.
    localIndex j = 0;
    for ( ; j < numFaces - 1; ++j )
    {
      ++numUniqueFaces;
      j += facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 );
    }

    numUniqueFaces += j == numFaces - 1;
  } );

  // At this point uniqueFaceOffsets[ i ] holds the number of unique face associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique face offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueFaceOffsets.begin(), uniqueFaceOffsets.end() );
  return uniqueFaceOffsets.back();
}

/**
 * @brief Add an interior face to the element lists, the face to node map, and the element to face map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in] faceID the ID of the face to add.
 * @param [in] fb0 the FaceBuilder associated with the first element of the current face.
 * @param [in] fb1 the FaceBuilder associated with the second element of the current face.
 * @param [in/out] elemRegionList the face to element region map.
 * @param [in/out] elemSubRegionList the face to element subregion map.
 * @param [in/out] elemList the face to element map.
 * @param [in/out] nodeList the face to node map.
 */
void addInteriorFace( ElementRegionManager & elementManager,
                      localIndex const faceID,
                      FaceBuilder const & fb0,
                      FaceBuilder const & fb1,
                      arrayView2d< localIndex > const & elemRegionList,
                      arrayView2d< localIndex > const & elemSubRegionList,
                      arrayView2d< localIndex > const & elemList,
                      array1d< array1d< localIndex > > const & nodeList )
{
  // Handle the first element.
  {
    localIndex const er = fb0.er;
    localIndex const esr = fb0.esr;
    localIndex const k = fb0.k;
    localIndex const elementLocalFaceIndex = fb0.elementLocalFaceIndex;

    // Get the subRegion associated with the element.
    CellElementSubRegion & subRegion = *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
    
    // The first element defines the node ordering for the face.
    subRegion.GetFaceNodes( k, elementLocalFaceIndex, nodeList[ faceID ] );
    
    // Add the face to the element to face map.
    subRegion.faceList()( k, elementLocalFaceIndex ) = faceID;

    // Populate the face to element maps.
    elemRegionList( faceID, 0 ) = er;
    elemSubRegionList( faceID, 0 ) = esr;
    elemList( faceID, 0 ) = k;
  }
  
  // Handle the second element.
  {
    localIndex const er = fb1.er;
    localIndex const esr = fb1.esr;
    localIndex const k = fb1.k;
    localIndex const elementLocalFaceIndex = fb1.elementLocalFaceIndex;

    CellElementSubRegion & subRegion = *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
    subRegion.faceList()( k, elementLocalFaceIndex ) = faceID;

    elemRegionList( faceID, 1 ) = er;
    elemSubRegionList( faceID, 1 ) = esr;
    elemList( faceID, 1 ) = k;
  }
}

/**
 * @brief Add a boundary face to the element lists, the face to node map, and the element to face map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in] faceID the ID of the face to add.
 * @param [in] fb0 the FaceBuilder associated with the first element of the current face.
 * @param [in/out] elemRegionList the face to element region map.
 * @param [in/out] elemSubRegionList the face to element subregion map.
 * @param [in/out] elemList the face to element map.
 * @param [in/out] nodeList the face to node map.
 */
void addBoundaryFace( ElementRegionManager & elementManager,
                      localIndex const faceID,
                      FaceBuilder const & fb,
                      arrayView2d< localIndex > const & elemRegionList,
                      arrayView2d< localIndex > const & elemSubRegionList,
                      arrayView2d< localIndex > const & elemList,
                      array1d< array1d< localIndex > > const & nodeList )
{
  localIndex const er = fb.er;
  localIndex const esr = fb.esr;
  localIndex const k = fb.k;
  localIndex const elementLocalFaceIndex = fb.elementLocalFaceIndex;

  // Get the subRegion associated with the element.
  CellElementSubRegion & subRegion = *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
  
  // Get the nodes associated with the face.
  subRegion.GetFaceNodes( k, elementLocalFaceIndex, nodeList[ faceID ] );

  // Add the face to the element to face map.
  subRegion.faceList()( k, elementLocalFaceIndex ) = faceID;

  // Populate the face to element maps.
  elemRegionList( faceID, 0 ) = er;
  elemSubRegionList( faceID, 0 ) = esr;
  elemList( faceID, 0 ) = k;

  elemRegionList( faceID, 1 ) = -1;
  elemSubRegionList( faceID, 1 ) = -1;
  elemList( faceID, 1 ) = -1;
}

/**
 * @brief Populate the face to element maps, the face to node map, and the element to face map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in] facesByLowestNode and array of size numNodes of arrays of FaceBuilders associated with each node.
 * @param [in] uniqueFaceOffsets an array containing the unique ID of the first face associated with each node.
 * @param [in/out] elemRegionList the face to element region map.
 * @param [in/out] elemSubRegionList the face to element subregion map.
 * @param [in/out] elemList the face to element map.
 * @param [in/out] nodeList the face to node map.
 * @param [in] maxFaceNodes the maximum number of nodes associated with any face.
 */
void populateMaps( ElementRegionManager & elementManager,
                   ArrayOfArraysView< FaceBuilder const > const & facesByLowestNode,
                   arrayView1d< localIndex const > const & uniqueFaceOffsets,
                   arrayView2d< localIndex > const & elemRegionList,
                   arrayView2d< localIndex > const & elemSubRegionList,
                   arrayView2d< localIndex > const & elemList,
                   array1d< array1d< localIndex > > const & nodeList,
                   localIndex const maxFaceNodes )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  GEOS_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );
  GEOS_ERROR_IF_NE( numUniqueFaces, elemRegionList.size( 0 ) );
  GEOS_ERROR_IF_NE( numUniqueFaces, elemSubRegionList.size( 0 ) );
  GEOS_ERROR_IF_NE( numUniqueFaces, elemList.size( 0 ) );
  GEOS_ERROR_IF_NE( numUniqueFaces, nodeList.size( 0 ) );

  // Need to be smarter about this. Should precalculate the number of nodes associated with each face
  // and store it in an array which would be used here.
  GEOSX_MARK_BEGIN("Reserving space in nodeList");
  for ( localIndex faceID = 0; faceID < numUniqueFaces; ++faceID )
  {
    nodeList[ faceID ].reserve( maxFaceNodes );
  }
  GEOSX_MARK_END("Reserving space in nodeList");

  // loop over all the nodes.
  forall_in_range< parallelHostPolicy >( 0, numNodes, [&]( localIndex const nodeID )
  {
    localIndex curFaceID = uniqueFaceOffsets[ nodeID ];
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );
    
    // loop over all the FaceBuilders associated with the node
    localIndex j = 0;
    for ( ; j < numFaces - 1; ++j )
    {
      // If two subsequent FaceBuilders compare equal then they describe an interior face.
      if ( facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 ) )
      {
        addInteriorFace( elementManager, curFaceID, facesByLowestNode( nodeID, j ), facesByLowestNode( nodeID, j + 1 ), elemRegionList, elemSubRegionList, elemList, nodeList );
        ++j;
      }

      // Otherwise it's a boundary face.
      else
      {
        addBoundaryFace( elementManager, curFaceID, facesByLowestNode( nodeID, j ), elemRegionList, elemSubRegionList, elemList, nodeList );
      }

      ++curFaceID;
    }

    if ( j == numFaces - 1 )
    {
      addBoundaryFace( elementManager, curFaceID, facesByLowestNode( nodeID, j ), elemRegionList, elemSubRegionList, elemList, nodeList );
    }
  } );
}


void FaceManager::BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elementManager )
{
  GEOSX_MARK_FUNCTION;

  m_toElements.setElementRegionManager( elementManager );

  constexpr int MAX_FACES_PER_NODE = 20;
  localIndex const numNodes = nodeManager->size();
  
  ArrayOfArrays<FaceBuilder> facesByLowestNode( numNodes, 2 * MAX_FACES_PER_NODE );
  localIndex const maxFaceNodes = createFacesByLowestNode( *elementManager, facesByLowestNode );
  
  array1d< localIndex > uniqueFaceOffsets( numNodes + 1 );
  localIndex const numFaces = calculateTotalNumberOfFaces( facesByLowestNode, uniqueFaceOffsets );
  
  GEOSX_MARK_BEGIN("FaceManager::resize");
  resize( numFaces );
  GEOSX_MARK_END("FaceManager::resize");

  populateMaps( *elementManager,
                facesByLowestNode,
                uniqueFaceOffsets,
                elementRegionList(),
                elementSubRegionList(),
                elementList(),
                nodeList(),
                maxFaceNodes );

  // First create the sets
  auto const & nodeSets = nodeManager->sets()->wrappers();
  for ( localIndex i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    CreateSet( setName );
  }

  // Then loop over them in parallel and fill them in.
  GEOSX_MARK_BEGIN("Set construction");
  forall_in_range<parallelHostPolicy>( 0, nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    const set<localIndex>& targetSet = nodeManager->sets()->getReference<set<localIndex>>( setName );
    ConstructSetFromSetAndMap( targetSet, m_nodeList, setName );
  } );
  GEOSX_MARK_END("Set construction");

  SetDomainBoundaryObjects( nodeManager );

  computeGeometry( nodeManager );
}


void FaceManager::computeGeometry( NodeManager const * const nodeManager )
{
  real64_array & faceArea  = getReference<real64_array>( viewKeyStruct::faceAreaString);
  r1_array & faceNormal = getReference<r1_array>( viewKeyStruct::faceNormalString);
  r1_array & faceCenter = getReference<r1_array>( viewKeyStruct::faceCenterString);
  r1_array const & X = nodeManager->referencePosition();

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  forall_in_range< parallelHostPolicy >( 0, this->size(), [&]( localIndex const faceID )
  {
    faceArea[ faceID ] = computationalGeometry::Centroid_3DPolygon( m_nodeList[ faceID ],
                                                                    X,
                                                                    faceCenter[ faceID ],
                                                                    faceNormal[ faceID ] );
  } );
}

void FaceManager::SetDomainBoundaryObjects( NodeManager * const nodeManager )
{
  // Set value of domainBounaryIndicator to one if it is found to have only one elements that it
  // is connected to.
  integer_array & faceDomainBoundaryIndicator = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);
  faceDomainBoundaryIndicator = 0;

  array2d<localIndex> const & elemRegionList = this->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = this->elementSubRegionList();
  array2d<localIndex> const & elemList = this->elementList();

  forall_in_range<parallelHostPolicy>( 0, size(), [&]( localIndex const kf )
  {
    if( elemRegionList[kf][1] == -1 )
    {
      faceDomainBoundaryIndicator(kf) = 1;
    }
  } );

  integer_array & nodeDomainBoundaryIndicator = nodeManager->getReference<integer_array>(nodeManager->viewKeys.domainBoundaryIndicator);
  nodeDomainBoundaryIndicator = 0;

  OrderedVariableOneToManyRelation const & faceToNodesMap = this->nodeList();

  forall_in_range<parallelHostPolicy>( 0, size(), [&]( localIndex const k )
  {
    if( faceDomainBoundaryIndicator[k] == 1 )
    {
      arrayView1d<localIndex> const& nodelist = faceToNodesMap[k];
      for( localIndex a=0 ; a< nodelist.size() ; ++a )
      {
        nodeDomainBoundaryIndicator[nodelist[a]] = 1;
      }
    }
  } );
}


void FaceManager::SetIsExternal()
{
  integer_array const &
  isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  m_isExternal = 0;
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k]==1 )
    {
      m_isExternal[k] = 1;
    }
  }
}

//void
//FaceManager::
//SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject )
//{
//  array1d< localIndex_array > const & faceToNodes = this->getReference< array1d< localIndex_array > >( viewKeys.nodeList );
//  globalIndex_array const & nodalGlobalIndex = compositionalObject->m_localToGlobalMap;
//  integer_array const & isDomainBoundary = this->getReference<integer_array>(viewKeys.isDomainBoundary);
//
//  mpiBuffer buffer;
//
//  localIndex numFaces;
//  for( localIndex k=0 ; k<size() ; ++k )
//  {
//    if( isDomainBoundary[k] == 1 )
//    {
//    }
//  }
//  for( localIndex k=0 ; k<size() ; ++k )
//  {
//    if( isDomainBoundary[k] == 1 )
//    {
//      CommBufferOps::Pack( buffer, )
//    }
//  }
//}

localIndex FaceManager::getMaxFaceNodes() const
{
  localIndex max_size = 0;
  OrderedVariableOneToManyRelation const & faces_to_nodes = nodeList();
  for(localIndex kf =0 ; kf < size() ; ++kf )
  {
    const localIndex size = faces_to_nodes[kf].size();
    if ( size > max_size )
    {
      max_size = size;
    }
  }

  return max_size;
}

void FaceManager::SortAllFaceNodes( NodeManager const * const nodeManager,
                                    ElementRegionManager const *const elemManager )
{
  GEOSX_MARK_FUNCTION;

  array2d<localIndex> const & elemRegionList = elementRegionList();
  array2d<localIndex> const & elemSubRegionList = elementSubRegionList();
  array2d<localIndex> const & elemList = elementList();
  arrayView1d<R1Tensor const> const & X = nodeManager->referencePosition();

  const indexType max_face_nodes = getMaxFaceNodes();
  GEOS_ERROR_IF( max_face_nodes >= MAX_FACE_NODES, "More nodes on a face than expected!" );

  elemManager->forElementSubRegions<CellElementSubRegion>([X] (CellElementSubRegion const * const subRegion)
  { subRegion->calculateElementCenters(X); });
  
  forall_in_range<parallelHostPolicy>( 0, size(), [&]( localIndex const kf ) -> void
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion( elemRegionList[kf][0] );
    CellElementSubRegion const * const subRegion = elemRegion->GetSubRegion<CellElementSubRegion>( elemSubRegionList[kf][0] );
    R1Tensor const elementCenter = subRegion->getElementCenter()( elemList[kf][0] );
    const localIndex numFaceNodes = nodeList()[kf].size();
    arrayView1d<localIndex> & faceNodes = nodeList()[kf];
    SortFaceNodes( X, elementCenter, faceNodes, numFaceNodes );
  } );
}

void FaceManager::SortFaceNodes( arrayView1d<R1Tensor const> const & X,
                                 R1Tensor const & elementCenter,
                                 arrayView1d<localIndex> const & faceNodes,
                                 localIndex const numFaceNodes )
{
  localIndex const firstNodeIndex = faceNodes[0];

  // get face center (average vertex location)
  R1Tensor fc(0);
  for( localIndex n =0 ; n < numFaceNodes ; ++n)
  {
    fc += X[faceNodes[n]];
  }
  fc /= realT(numFaceNodes);

  R1Tensor ex, ey, ez;
  // Approximate face normal direction (unscaled)

  if (numFaceNodes == 2)  //2D only.
  {
    ex = X[faceNodes[1]];
    ex -= X[faceNodes[0]];
    ey = elementCenter;
    ey -= fc;

    ez.Cross(ex, ey);
    // The element should be on the right hand side of the vector from node 0 to
    // node 1.
    // This ensure that the normal vector of an external face points to outside
    // the element.
    if (ez[2] > 0)
    {
      localIndex itemp = faceNodes[0];
      faceNodes[0] = faceNodes[1];
      faceNodes[1] = itemp;
    }
  }
  else
  {
    ez = fc;
    ez -= elementCenter;

    /// Approximate in-plane axis
    ex = X[faceNodes[0]];
    ex -= fc;

    ex /= ex.L2_Norm();
    ey.Cross(ez, ex);
    ey /= ey.L2_Norm();

    std::pair<realT, localIndex> thetaOrder[MAX_FACE_NODES];

    /// Sort nodes counterclockwise around face center
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      R1Tensor v = X[faceNodes[n]];
      v -= fc;
      thetaOrder[n] = std::pair<realT, localIndex>(atan2(Dot(v,ey),Dot(v,ex)),faceNodes[n]);
    }

    std::sort(thetaOrder, thetaOrder + numFaceNodes);

    // Reorder nodes on face
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    localIndex tempFaceNodes[MAX_FACE_NODES];

    localIndex firstIndexIndex = 0;
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      tempFaceNodes[n] = thetaOrder[n].second;
      if( tempFaceNodes[n] == firstNodeIndex )
      {
        firstIndexIndex = n;
      }
    }

    for( localIndex n=0 ; n < numFaceNodes ; ++n)
    {
      const localIndex index = (firstIndexIndex + n) % numFaceNodes;
      faceNodes[n] = tempFaceNodes[index];
    }
  }
}


void FaceManager::ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                                   std::vector< std::vector< globalIndex > > & globalFaceNodes )
{
  GEOSX_MARK_FUNCTION;
  nodeManager->CheckTypeID( typeid( NodeManager ) );

  localIndex const numFaces = size();

  arrayView1d< arrayView1d< localIndex const > const > const & faceNodes = this->nodeList().toViewConst();
  arrayView1d< integer const > const & isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  globalFaceNodes.resize( numFaces );

  forall_in_range< parallelHostPolicy >( 0, numFaces, [&]( localIndex const & faceID )
  {
    std::vector< globalIndex > & curFaceGlobalNodes = globalFaceNodes[ faceID ];

    if( isDomainBoundary( faceID ) )
    {
      localIndex const numNodes = faceNodes[faceID].size();
      curFaceGlobalNodes.resize( numNodes );

      for ( localIndex a = 0; a < numNodes ; ++a )
      {
        curFaceGlobalNodes[ a ]= nodeManager->m_localToGlobalMap( faceNodes[ faceID ][ a ] );
      }

      std::sort( curFaceGlobalNodes.begin(), curFaceGlobalNodes.end() );
    }
    else
    {
      curFaceGlobalNodes.resize( 0 );
    }
  } );
}



void FaceManager::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementSubRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementListString));
}


localIndex FaceManager::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex FaceManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template<bool DOPACK>
localIndex FaceManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::nodeListString) );

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_nodeList.Base(),
                                         m_unmappedGlobalIndicesInToNodes,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_nodeList.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::edgeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_edgeList.Base(),
                                         m_unmappedGlobalIndicesInToEdges,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_edgeList.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::elementListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         this->m_toElements,
                                         packList,
                                         m_toElements.getElementRegionManager() );


  return packedSize;
}



localIndex FaceManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const overwriteDownMaps )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_nodeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->m_globalToLocalMap,
                                     m_nodeList.RelatedObjectGlobalToLocal() );


  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_edgeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->m_globalToLocalMap,
                                     m_edgeList.RelatedObjectGlobalToLocal() );


  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::elementListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );



  return unPackedSize;
}

void FaceManager::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_nodeList,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_edgeList,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

}


void FaceManager::depopulateUpMaps( std::set<localIndex> const & receivedFaces,
                                    ElementRegionManager const & elemRegionManager )
{
  for( auto const & targetIndex : receivedFaces )
  {
    for( localIndex k=0 ; k<m_toElements.m_toElementRegion.size(1) ; ++k )
    {
      localIndex const elemRegionIndex    = m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex = m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex          = m_toElements.m_toElementIndex[targetIndex][k];

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 && elemIndex!=-1 )
      {
        CellElementSubRegion const * subRegion = elemRegionManager.GetRegion(elemRegionIndex)->
                                                 GetSubRegion<CellElementSubRegion>(elemSubRegionIndex);
        array2d<localIndex> const & downmap = subRegion->faceList();
        bool hasTargetIndex = false;

        for( localIndex a=0 ; a<downmap.size(1) ; ++a )
        {
          localIndex const compositeLocalIndex = downmap[elemIndex][a];
          if( compositeLocalIndex==targetIndex )
          {
            hasTargetIndex=true;
          }
        }
        if( !hasTargetIndex )
        {
          m_toElements.m_toElementRegion[targetIndex][k] = -1;
          m_toElements.m_toElementSubRegion[targetIndex][k] = -1;
          m_toElements.m_toElementIndex[targetIndex][k] = -1;
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceManager, std::string const &, ManagedGroup * const )

}
