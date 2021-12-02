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
 * @file MeshLevel.cpp
 */

#include "MeshLevel.hpp"

#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name,
                      Group * const parent ):
  Group( name, parent ),
  m_nodeManager( groupStructKeys::nodeManagerString, this ),
  m_edgeManager( groupStructKeys::edgeManagerString, this ),
  m_faceManager( groupStructKeys::faceManagerString, this ),
  m_elementManager( groupStructKeys::elemManagerString, this ),
  m_embSurfNodeManager( groupStructKeys::embSurfNodeManagerString, this ),
  m_embSurfEdgeManager( groupStructKeys::embSurfEdgeManagerString, this )

{

  registerGroup( groupStructKeys::nodeManagerString, &m_nodeManager );

  registerGroup( groupStructKeys::edgeManagerString, &m_edgeManager );


  registerGroup< FaceManager >( groupStructKeys::faceManagerString, &m_faceManager );
  m_faceManager.nodeList().setRelatedObject( m_nodeManager );


  registerGroup< ElementRegionManager >( groupStructKeys::elemManagerString, &m_elementManager );

  registerGroup< EdgeManager >( groupStructKeys::embSurfEdgeManagerString, &m_embSurfEdgeManager );

  registerGroup< EmbeddedSurfaceNodeManager >( groupStructKeys::embSurfNodeManagerString, &m_embSurfNodeManager );

  registerWrapper< integer >( viewKeys.meshLevel );
}


MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel const & source,
                      int const order ):
  MeshLevel( name, parent )
{

//  {
//  Group & thisGroup = *this;
//  Group const & sourceGroup = source;
//  thisGroup = sourceGroup;
//  this->getParent( parent );
//  }


  localIndex numNodes = source.m_nodeManager.size();
  // find out how many node there must be on this rank
  m_nodeManager.resize(numNodes);

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager.referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const refPosNew = m_nodeManager.referencePosition().toView();

  {
    Group & nodeSets = m_nodeManager.sets();
    SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
    allNodes.reserve( m_nodeManager.size() );

    for( localIndex a=0; a<m_nodeManager.size(); ++a )
    {
      allNodes.insert( a );
    }

    // fill coordinate locations here??
    forAll<parallelDevicePolicy<>>( m_nodeManager.size(),
                                    [=]( localIndex const a )
    {
      for( localIndex i=0; i<3; ++i )
      {
        refPosNew(a,i) = refPosSource(a,i); // this needs to be another loop with a linear combination of values.
      }
    });
  }



  source.m_elementManager.forElementRegions<CellElementRegion>([&]( CellElementRegion const & sourceRegion )
  {
    CellElementRegion & region = *(dynamic_cast<CellElementRegion *>( m_elementManager.createChild( sourceRegion.getCatalogName(),
                                                                                                        sourceRegion.getName() ) ) );

    region.addCellBlockNames( sourceRegion.getCellBlockNames() );

    sourceRegion.forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion const & sourceSubRegion )
    {

      CellElementSubRegion & newSubRegion = region.getSubRegions().registerGroup< CellElementSubRegion >( sourceSubRegion.getName() );
      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      newSubRegion.resize( sourceSubRegion.size() );


      Group & sets = newSubRegion.sets();
      SortedArray< localIndex > & allElems  = sets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();


      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      array2d< localIndex, cells::NODE_MAP_PERMUTATION > elemsToNodesNew = sourceSubRegion.nodeList();

      localIndex const numNodesPerElem = 8;// change to pow( 2+( order>1 ? order-1 : 0 ), 3 );
      elemsToNodesNew.resize( elemsToNodesNew.size(0), numNodesPerElem );

      for( localIndex k=0; k<newSubRegion.size(); ++k )
      {
        allElems.insert( k );
        for( localIndex a=0; a<numNodesPerElem; ++a )
        {
          elemsToNodesNew(k,a) = elemsToNodesSource(k,a); // need the logic to map to the nodes here
        }
      }



    });
  });
}



MeshLevel::~MeshLevel()
{}

void MeshLevel::initializePostInitialConditionsPostSubGroups()
{
  m_elementManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    subRegion.calculateElementGeometricQuantities( m_nodeManager, m_faceManager );
  } );
}


void MeshLevel::generateAdjacencyLists( arrayView1d< localIndex const > const & seedNodeList,
                                        localIndex_array & nodeAdjacencyList,
                                        localIndex_array & edgeAdjacencyList,
                                        localIndex_array & faceAdjacencyList,
                                        ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > & elementAdjacencyList,
                                        integer const depth )
{
  NodeManager & nodeManager = getNodeManager();

  ArrayOfArraysView< localIndex const > const & nodeToElementRegionList = nodeManager.elementRegionList().toViewConst();

  ArrayOfArraysView< localIndex const > const & nodeToElementSubRegionList = nodeManager.elementSubRegionList().toViewConst();

  ArrayOfArraysView< localIndex const > const & nodeToElementList = nodeManager.elementList().toViewConst();


  FaceManager & faceManager = this->getFaceManager();
  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  ElementRegionManager & elemManager = this->getElemManager();

  std::set< localIndex > nodeAdjacencySet;
  std::set< localIndex > edgeAdjacencySet;
  std::set< localIndex > faceAdjacencySet;
  std::vector< std::vector< std::set< localIndex > > > elementAdjacencySet( elemManager.numRegions() );

  for( localIndex a=0; a<elemManager.numRegions(); ++a )
  {
    elementAdjacencySet[a].resize( elemManager.getRegion( a ).numSubRegions() );
  }

  nodeAdjacencySet.insert( seedNodeList.begin(), seedNodeList.end() );

  for( integer d=0; d<depth; ++d )
  {
    for( localIndex const nodeIndex : nodeAdjacencySet )
    {
      for( localIndex b=0; b<nodeToElementRegionList.sizeOfArray( nodeIndex ); ++b )
      {
        localIndex const regionIndex = nodeToElementRegionList[nodeIndex][b];
        localIndex const subRegionIndex = nodeToElementSubRegionList[nodeIndex][b];
        localIndex const elementIndex = nodeToElementList[nodeIndex][b];
        elementAdjacencySet[regionIndex][subRegionIndex].insert( elementIndex );
      }
    }

    for( typename dataRepository::indexType kReg=0; kReg<elemManager.numRegions(); ++kReg )
    {
      ElementRegionBase const & elemRegion = elemManager.getRegion( kReg );

      elemRegion.forElementSubRegionsIndex< CellElementSubRegion,
                                            WellElementSubRegion >( [&]( localIndex const kSubReg,
                                                                         auto const & subRegion )
      {
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList();
        arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
        for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
        {
          for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
          {
            nodeAdjacencySet.insert( elemsToNodes[elementIndex][a] );
          }

          for( localIndex a=0; a<elemsToFaces.size( 1 ); ++a )
          {
            faceAdjacencySet.insert( elemsToFaces[elementIndex][a] );

            localIndex const faceID = elemsToFaces[elementIndex][a];
            localIndex const numEdges = faceToEdges.sizeOfArray( faceID );
            for( localIndex b=0; b<numEdges; ++b )
            {
              edgeAdjacencySet.insert( faceToEdges( faceID, b ));
            }

          }

        }
      } );
    }
  }

  nodeAdjacencyList.resize( LvArray::integerConversion< localIndex >( nodeAdjacencySet.size()));
  std::copy( nodeAdjacencySet.begin(), nodeAdjacencySet.end(), nodeAdjacencyList.begin() );

  edgeAdjacencyList.resize( LvArray::integerConversion< localIndex >( edgeAdjacencySet.size()));
  std::copy( edgeAdjacencySet.begin(), edgeAdjacencySet.end(), edgeAdjacencyList.begin() );

  faceAdjacencyList.resize( LvArray::integerConversion< localIndex >( faceAdjacencySet.size()));
  std::copy( faceAdjacencySet.begin(), faceAdjacencySet.end(), faceAdjacencyList.begin() );

  for( localIndex kReg=0; kReg<elemManager.numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = elemManager.getRegion( kReg );

    for( localIndex kSubReg = 0; kSubReg < elemRegion.numSubRegions(); ++kSubReg )
    {
      elementAdjacencyList[kReg][kSubReg].get().clear();
      elementAdjacencyList[kReg][kSubReg].get().resize( LvArray::integerConversion< localIndex >( elementAdjacencySet[kReg][kSubReg].size()) );
      std::copy( elementAdjacencySet[kReg][kSubReg].begin(),
                 elementAdjacencySet[kReg][kSubReg].end(),
                 elementAdjacencyList[kReg][kSubReg].get().begin() );

    }
  }

}

} /* namespace geosx */
