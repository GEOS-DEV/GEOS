// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * MeshLevel.cpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#include "MeshLevel.hpp"
#include "NodeManager.hpp"
//#include "EdgeManager.hpp"
#include "FaceManager.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name,
                      ManagedGroup * const parent ):
  ManagedGroup(name,parent),
  m_nodeManager( groupStructKeys::nodeManagerString,this),
  m_edgeManager( groupStructKeys::edgeManagerString,this),
  m_faceManager( groupStructKeys::faceManagerString,this),
  m_elementManager( groupStructKeys::elemManagerString,this)
{

  RegisterGroup( groupStructKeys::nodeManagerString, &m_nodeManager, false );

  RegisterGroup( groupStructKeys::edgeManagerString, &m_edgeManager, false );


  RegisterGroup<FaceManager>( groupStructKeys::faceManagerString, &m_faceManager, false );
  m_faceManager.nodeList().SetRelatedObject( &m_nodeManager );


  RegisterGroup<ElementRegionManager>( groupStructKeys::elemManagerString, &m_elementManager, false );


  RegisterViewWrapper<integer>( viewKeys.meshLevel );
}

MeshLevel::~MeshLevel()
{
  // TODO Auto-generated destructor stub
}



void MeshLevel::InitializePostSubGroups( ManagedGroup * const )
{

  m_elementManager.forCellBlocks([&]( CellBlockSubRegion * subRegion ) -> void
  {
    subRegion->nodeList().SetRelatedObject(&m_nodeManager);
    subRegion->faceList().SetRelatedObject(&m_faceManager);
  });

}


void MeshLevel::GenerateAdjacencyLists( localIndex_array & seedNodeList,
                                        localIndex_array & nodeAdjacencyList,
                                        localIndex_array & edgeAdjacencyList,
                                        localIndex_array & faceAdjacencyList,
                                        ElementRegionManager::ElementViewAccessor<localIndex_array>& elementAdjacencyList,
                                        integer const depth )
{
  NodeManager * const nodeManager = getNodeManager();

  array<lSet> const & nodeToElementRegionList = nodeManager->elementRegionList();

  array<lSet> const & nodeToElementSubRegionList = nodeManager->elementSubRegionList();

  array<lSet> const & nodeToElementList = nodeManager->elementList();


  FaceManager * const faceManager = this->getFaceManager();
  ElementRegionManager * const elemManager = this->getElemManager();

  localIndex_set nodeAdjacencySet;
  localIndex_set faceAdjacencySet;
  array< array< localIndex_set > > elementAdjacencySet;
  elementAdjacencySet.resize( elemManager->numRegions() );

  for( localIndex a=0 ; a<elemManager->numRegions() ; ++a )
  {
    elementAdjacencySet[a].resize( elemManager->GetRegion(a)->numSubRegions() );
  }

  nodeAdjacencySet.insert( seedNodeList.begin(), seedNodeList.end() );

  for( integer d=0 ; d<depth ; ++d )
  {
    for( auto const nodeIndex : nodeAdjacencySet )
    {
      for( localIndex b=0 ; b<nodeToElementRegionList[nodeIndex].size() ; ++b )
      {
        localIndex const regionIndex = nodeToElementRegionList[nodeIndex][b];
        localIndex const subRegionIndex = nodeToElementSubRegionList[nodeIndex][b];
        localIndex const elementIndex = nodeToElementList[nodeIndex][b];
        elementAdjacencySet[regionIndex][subRegionIndex].insert(elementIndex);
      }
    }

    for( typename dataRepository::indexType kReg=0 ; kReg<elemManager->numRegions() ; ++kReg  )
    {
      ElementRegion const * const elemRegion = elemManager->GetRegion(kReg);

      for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
      {
        CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);

        lArray2d const & elemsToNodes = subRegion->getReference<FixedOneToManyRelation>(subRegion->viewKeys().nodeList);
        for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
        {
          arrayView1d<localIndex const> const nodeList = elemsToNodes[elementIndex];
          for( localIndex a=0 ; a<elemsToNodes.size(1) ; ++a )
          {
            nodeAdjacencySet.insert(nodeList[a]);
          }
        }

        lArray2d const & elemsToFaces = subRegion->getReference<FixedOneToManyRelation>(subRegion->viewKeys().faceList);
        for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
        {
          arrayView1d<localIndex const> const faceList = elemsToFaces[elementIndex];
          for( localIndex a=0 ; a<elemsToFaces.size(1) ; ++a )
          {
            faceAdjacencySet.insert(faceList[a]);
          }
        }
      }
    }
    nodeAdjacencyList.clear();
    nodeAdjacencyList.resize( integer_conversion<localIndex>(nodeAdjacencySet.size()));
    std::copy(nodeAdjacencySet.begin(), nodeAdjacencySet.end(), nodeAdjacencyList.begin() );

  }

  nodeAdjacencyList.clear();
  nodeAdjacencyList.resize(integer_conversion<localIndex>(nodeAdjacencySet.size()));
  std::copy(nodeAdjacencySet.begin(), nodeAdjacencySet.end(), nodeAdjacencyList.begin() );

  faceAdjacencyList.clear();
  faceAdjacencyList.resize(integer_conversion<localIndex>(faceAdjacencySet.size()));
  std::copy(faceAdjacencySet.begin(), faceAdjacencySet.end(), faceAdjacencyList.begin() );


  for( typename dataRepository::indexType kReg=0 ; kReg<elemManager->numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(kReg);

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      elementAdjacencyList[kReg][kSubReg]->clear();
      elementAdjacencyList[kReg][kSubReg]->resize( integer_conversion<localIndex>(elementAdjacencySet[kReg][kSubReg].size()) );
      std::copy( elementAdjacencySet[kReg][kSubReg].begin(),
                 elementAdjacencySet[kReg][kSubReg].end(),
                 elementAdjacencyList[kReg][kSubReg]->begin() );

    }
  }

}

} /* namespace geosx */
