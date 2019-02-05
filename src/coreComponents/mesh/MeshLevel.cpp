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
{}


void MeshLevel::GenerateAdjacencyLists( localIndex_array & seedNodeList,
                                        localIndex_array & nodeAdjacencyList,
                                        localIndex_array & edgeAdjacencyList,
                                        localIndex_array & faceAdjacencyList,
                                        ElementRegionManager::ElementViewAccessor<ReferenceWrapper<localIndex_array>>& elementAdjacencyList,
                                        integer const depth )
{
  NodeManager * const nodeManager = getNodeManager();

  array1d<array1d<localIndex>> const & nodeToElementRegionList = nodeManager->elementRegionList();

  array1d<array1d<localIndex>> const & nodeToElementSubRegionList = nodeManager->elementSubRegionList();

  array1d<array1d<localIndex>> const & nodeToElementList = nodeManager->elementList();


  FaceManager * const faceManager = this->getFaceManager();
  array1d< array1d< localIndex > > const & faceToEdges = faceManager->edgeList();

  ElementRegionManager * const elemManager = this->getElemManager();

  localIndex_set nodeAdjacencySet;
  localIndex_set edgeAdjacencySet;
  localIndex_set faceAdjacencySet;
  array1d< array1d< localIndex_set > > elementAdjacencySet;
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
        CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion<CellBlockSubRegion>(kSubReg);

        array2d<localIndex> const & elemsToNodes = subRegion->nodeList();
        array2d<localIndex> const & elemsToFaces = subRegion->faceList();
        array2d<localIndex> const & elemsToEdges = subRegion->edgeList();
        for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
        {
          for( localIndex a=0 ; a<elemsToNodes.size(1) ; ++a )
          {
            nodeAdjacencySet.insert(elemsToNodes[elementIndex][a]);
          }

          for( localIndex a=0 ; a<elemsToFaces.size(1) ; ++a )
          {
            faceAdjacencySet.insert(elemsToFaces[elementIndex][a]);

            array1d<localIndex> const & edgeList = faceToEdges[elemsToFaces[elementIndex][a]];
            for( localIndex b=0 ; b<edgeList.size() ; ++b )
            {
              edgeAdjacencySet.insert(edgeList[b]);
            }

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

  edgeAdjacencyList.clear();
  edgeAdjacencyList.resize(integer_conversion<localIndex>(edgeAdjacencySet.size()));
  std::copy(edgeAdjacencySet.begin(), edgeAdjacencySet.end(), edgeAdjacencyList.begin() );

  faceAdjacencyList.clear();
  faceAdjacencyList.resize(integer_conversion<localIndex>(faceAdjacencySet.size()));
  std::copy(faceAdjacencySet.begin(), faceAdjacencySet.end(), faceAdjacencyList.begin() );


  for( typename dataRepository::indexType kReg=0 ; kReg<elemManager->numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(kReg);

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      elementAdjacencyList[kReg][kSubReg].get().clear();
      elementAdjacencyList[kReg][kSubReg].get().resize( integer_conversion<localIndex>(elementAdjacencySet[kReg][kSubReg].size()) );
      std::copy( elementAdjacencySet[kReg][kSubReg].begin(),
                 elementAdjacencySet[kReg][kSubReg].end(),
                 elementAdjacencyList[kReg][kSubReg].get().begin() );

    }
  }

}

} /* namespace geosx */
