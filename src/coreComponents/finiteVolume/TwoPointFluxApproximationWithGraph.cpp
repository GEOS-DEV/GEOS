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

#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "mesh/GraphFromText.hpp"
#include "mesh/GraphEdge.hpp"

#include "mesh/GraphVertexFace.hpp"


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

  
  array1d<GraphEdge*> edges = graph->getEdges();
  array1d<std::shared_ptr<GraphVertex>> vertices = graph->getVertices();
  stencil.reserve( edges.size() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.ConstructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const localToGlobalMap =
    elemManager.ConstructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString );
   
  for( localIndex i=0; i<edges.size(); i++)
  {
    //std::cout<<localToGlobalMap[0][0][edges[i]->getVertex1()->getIndice()]<<"\n";

        //std::cout<< elemGhostRank[edges[i]->getVertex1()->getRegionIndex()][edges[i]->getVertex1()->getSubRegionIndex()][edges[i]->getVertex1()->getIndice()] << " " << elemGhostRank[edges[i]->getVertex2()->getRegionIndex()][edges[i]->getVertex2()->getSubRegionIndex()][edges[i]->getVertex2()->getIndice()] << "\n";
    //std::cout<< edges[i]->getVertex1()->getIndice() << " " << edges[i]->getVertex2()->getIndice() << "\n\n";

    if( edges[i]->getVertex1()->getGhostIndex() >= 0 && 
        edges[i]->getVertex2()->getGhostIndex() >= 0)
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

void TwoPointFluxApproximationWithGraph::computeBoundaryStencil( MeshLevel & mesh,
                                                        string const & setName,
                                                        SortedArrayView< localIndex const > const & faceSet ) const
{
const Group * tmp = this->GetGroupByPath( m_graphString );
  GEOSX_ERROR_IF( tmp == nullptr, "Can't find group with path " << m_graphString );
  const GraphFromText * graph = Group::group_cast< const GraphFromText * >( tmp );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  BoundaryStencil & stencil = getStencil< BoundaryStencil >( mesh, setName );
  
  array1d<GraphEdge*> edges = graph->getBoundaryEdges();
  array1d<std::shared_ptr<GraphVertex>> vertices = graph->getVertices();
  stencil.reserve( faceSet.size() );

  constexpr localIndex numPts = BoundaryStencil::NUM_POINT_IN_FLUX;

  stackArray1d< localIndex, numPts > stencilRegionIndices( numPts );
  stackArray1d< localIndex, numPts > stencilSubRegionIndices( numPts );
  stackArray1d< localIndex, numPts > stencilElemOrFaceIndices( numPts );
  stackArray1d< real64, numPts > stencilWeights( numPts );


  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.ConstructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const localToGlobalMap =
    elemManager.ConstructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString );
  
  for( localIndex kf : faceSet )
  {
    for( localIndex i=0; i<edges.size(); i++)
    {
      std::shared_ptr<GraphVertexFace> face = std::dynamic_pointer_cast<GraphVertexFace>(edges[i]->getVertex2());
      //std::cout<<edges[i]->getVertex1()->getGlobalVertexIndex() << " " << face->getCorrespondingId() << "\n";

      if( (edges[i]->getVertex1()->getGhostIndex() >= 0) || kf != face->getCorrespondingId())
      {
        //std::cout<<"Ghosted\n";
      }
      else
      {
        //std::cout << kf <<" "<<edges[i]->getVertex1()->getLocalVertexIndex()<<"\n";
        stencilRegionIndices[BoundaryStencil::Order::ELEM] = LvArray::integerConversion< localIndex >(edges[i]->getVertex1()->getRegionIndex());
        stencilSubRegionIndices[BoundaryStencil::Order::ELEM] = LvArray::integerConversion< localIndex >(edges[i]->getVertex1()->getSubRegionIndex());
        stencilElemOrFaceIndices[BoundaryStencil::Order::ELEM] = LvArray::integerConversion< localIndex >(edges[i]->getVertex1()->getLocalVertexIndex());
        stencilWeights[BoundaryStencil::Order::ELEM] = edges[i]->getTransmissibility();

        stencilRegionIndices[BoundaryStencil::Order::FACE] = LvArray::integerConversion< localIndex >(edges[i]->getVertex2()->getRegionIndex());
        stencilSubRegionIndices[BoundaryStencil::Order::FACE] = LvArray::integerConversion< localIndex >(edges[i]->getVertex2()->getSubRegionIndex());
        stencilElemOrFaceIndices[BoundaryStencil::Order::FACE] = kf;
        stencilWeights[BoundaryStencil::Order::FACE] = -edges[i]->getTransmissibility();

        //std::cout << stencilElemOrFaceIndices[0] << " " << stencilElemOrFaceIndices[1] << " " << "\n";

        stencil.add( stencilRegionIndices.size(),
                     stencilRegionIndices.data(),
                     stencilSubRegionIndices.data(),
                     stencilElemOrFaceIndices.data(),
                     stencilWeights.data(),
                     kf);
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( FluxApproximationBase, TwoPointFluxApproximationWithGraph, std::string const &, Group * const )

}
