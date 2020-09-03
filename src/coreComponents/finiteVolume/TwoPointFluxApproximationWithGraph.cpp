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
 * @file TwoPointFluxApproximationWithGraph.cpp
 *
 */

#include <unistd.h>
#include <limits.h>

#include "TwoPointFluxApproximationWithGraph.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "mesh/GraphFromText.hpp"
#include "mesh/GraphEdge.hpp"


namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximationWithGraph::TwoPointFluxApproximationWithGraph( std::string const & name,
                                                      Group * const parent )
  : TwoPointFluxApproximation( name, parent )
{
    registerWrapper( viewKeyStruct::graphString, &m_graphString )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the graph to link to the graph." );


}

void TwoPointFluxApproximationWithGraph::computeCellStencil( MeshLevel & mesh) const
{
  const Group * tmp = this->GetGroupByPath( m_graphString );
  GEOSX_ERROR_IF( tmp == nullptr, "Can't find group with path " << m_graphString );
  const GraphFromText * graph = Group::group_cast< const GraphFromText * >( tmp );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  CellElementStencilTPFA & stencil = getStencil< CellElementStencilTPFA >( mesh, viewKeyStruct::cellStencilString );

  
  std::vector<GraphEdge*> edges = graph->getEdges();
  std::vector<GraphVertex*> vertices = graph->getVertices();
  stencil.reserve( edges.size() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.ConstructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const localToGlobalMap =
    elemManager.ConstructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString );
  /*
  std::vector <GraphVertex*> new_vertex;
  for( long unsigned int i=0; i<vertices.size(); i++)
  {
    new_vertex.push_back(vertices[i]);
  }
  for (long int i = 0; i<elemGhostRank[0][0].size(); i++)
  {
    std::cout<<i<<"\n";
    std::cout<<localToGlobalMap[0][0][i]<<"\n";
    std::cout<<"\n";
    for (long unsigned int j = 0; j<new_vertex.size();j++)
    {
      if (localToGlobalMap[0][0][i] == new_vertex[i]->getEdgeIndex())
      {
        new_vertex.erase(new_vertex.begin()+i);
      }
    }
  }
  for ( long unsigned int i=0; i<new_vertex.size(); i++)
  {
    std::cout<<new_vertex[i]->getEdgeIndex()<<" ";
    for (long unsigned int j = 0; j<vertices.size();j++)
    {
      if (vertices[j]->getEdgeIndex == new_vertex[i]->getEdgeIndex())
      {
        graph->RemoveVertex(vertices[j]);
      }
    }
  std::cout<<"\n";
  }
  */
 
  for( long unsigned int i=0; i<edges.size(); i++)
  {
    //std::cout<<localToGlobalMap[0][0][edges[i]->getVertex1()->getIndice()]<<"\n";

        //std::cout<< elemGhostRank[edges[i]->getVertex1()->getRegionIndex()][edges[i]->getVertex1()->getSubRegionIndex()][edges[i]->getVertex1()->getIndice()] << " " << elemGhostRank[edges[i]->getVertex2()->getRegionIndex()][edges[i]->getVertex2()->getSubRegionIndex()][edges[i]->getVertex2()->getIndice()] << "\n";
    //std::cout<< edges[i]->getVertex1()->getIndice() << " " << edges[i]->getVertex2()->getIndice() << "\n\n";

    if( edges[i]->getVertex1()->getGhostIndex() >= 0 &&
        edges[i]->getVertex2()->getGhostIndex() >= 0 )
    {
      //std::cout<<"Ghosted\n";
    }
    else
    {
      stackArray1d< localIndex, 2 > regionIndex( 2 );
      stackArray1d< localIndex, 2 > subRegionIndex( 2 );
      stackArray1d< localIndex, 2 > elementIndex( 2 );
      stackArray1d< real64, 2 > stencilWeights( 2 );
      regionIndex[0]=LvArray::integerConversion< localIndex >(edges[i]->getVertex1()->getRegionIndex());
      subRegionIndex[0]=LvArray::integerConversion< localIndex >(edges[i]->getVertex1()->getSubRegionIndex());
      elementIndex[0]=LvArray::integerConversion< localIndex >(edges[i]->getVertex1()->getLocalVertexIndex());
      stencilWeights[0] = edges[i]->getTransmissibility();
 
      regionIndex[1]=LvArray::integerConversion< localIndex >(edges[i]->getVertex2()->getRegionIndex());
      subRegionIndex[1]=LvArray::integerConversion< localIndex >(edges[i]->getVertex2()->getSubRegionIndex());
      elementIndex[1]=LvArray::integerConversion< localIndex >(edges[i]->getVertex2()->getLocalVertexIndex());
      stencilWeights[1]=-edges[i]->getTransmissibility();
    localIndex faceIndex = LvArray::integerConversion< localIndex >(edges[i]->getEdgeIndex());
      stencil.add( 2,
                 regionIndex.data(),
                 subRegionIndex.data(),
                 elementIndex.data(),
                 stencilWeights.data(),
                 faceIndex);
    }
  }

}



REGISTER_CATALOG_ENTRY( FluxApproximationBase, TwoPointFluxApproximationWithGraph, std::string const &, Group * const )

}
