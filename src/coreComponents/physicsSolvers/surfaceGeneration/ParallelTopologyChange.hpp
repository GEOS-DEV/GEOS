/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParallelTopologyChange.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_
#define GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"

#define PARALLEL_TOPOLOGY_CHANGE_METHOD 1

namespace geos
{
class MeshLevel;
class NeighborCommunicator;

struct ModifiedObjectLists
{
  std::set< localIndex > newNodes;
  std::set< localIndex > newEdges;
  std::set< localIndex > newFaces;
  std::set< localIndex > modifiedNodes;
  std::set< localIndex > modifiedEdges;
  std::set< localIndex > modifiedFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newElements;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > modifiedElements;

  void clearNewFromModified()
  {
    for( localIndex const a : newNodes )
    {
      modifiedNodes.erase( a );
    }

    for( localIndex const a : newEdges )
    {
      modifiedEdges.erase( a );
    }

    for( localIndex const a : newFaces )
    {
      modifiedFaces.erase( a );
    }
  }

  void insert( ModifiedObjectLists const & modifiedObjects )
  {
    newNodes.insert( modifiedObjects.newNodes.begin(),
                     modifiedObjects.newNodes.end() );
    modifiedNodes.insert( modifiedObjects.modifiedNodes.begin(),
                          modifiedObjects.modifiedNodes.end() );

    newEdges.insert( modifiedObjects.newEdges.begin(),
                     modifiedObjects.newEdges.end() );
    modifiedEdges.insert( modifiedObjects.modifiedEdges.begin(),
                          modifiedObjects.modifiedEdges.end() );

    newFaces.insert( modifiedObjects.newFaces.begin(),
                     modifiedObjects.newFaces.end() );
    modifiedFaces.insert( modifiedObjects.modifiedFaces.begin(),
                          modifiedObjects.modifiedFaces.end() );

    for( auto & iter : modifiedObjects.newElements )
    {
      std::pair< localIndex, localIndex > const & key = iter.first;
      std::set< localIndex > const & values = iter.second;
      newElements[key].insert( values.begin(), values.end() );
    }

    for( auto & iter : modifiedObjects.modifiedElements )
    {
      std::pair< localIndex, localIndex > const & key = iter.first;
      std::set< localIndex > const & values = iter.second;
      modifiedElements[key].insert( values.begin(), values.end() );
    }

  }
};

namespace parallelTopologyChange
{

void synchronizeTopologyChange( MeshLevel * const mesh,
                                stdVector< NeighborCommunicator > & neighbors,
                                ModifiedObjectLists & modifiedObjects,
                                ModifiedObjectLists & receivedObjects,
                                int mpiCommOrder );



struct TopologyChangeStepData
{
  void init( ElementRegionManager const & elemManager )
  {
    m_nodes.resize( 0 );
    m_edges.resize( 0 );
    m_faces.resize( 0 );
    m_elements.resize( elemManager.numRegions() );
    m_elementsView.resize( elemManager.numRegions() );
    m_elementsData.resize( elemManager.numRegions() );
    m_size = 0;

    for( localIndex er=0; er<elemManager.numRegions(); ++er )
    {
      ElementRegionBase const & elemRegion = elemManager.getRegion( er );
      m_elements[er].resize( elemRegion.numSubRegions());
      m_elementsView[er].resize( elemRegion.numSubRegions());
      m_elementsData[er].resize( elemRegion.numSubRegions());
      for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
      {
        m_elementsData[er][esr].resize( 0 );
        m_elements[er][esr].set( m_elementsData[er][esr] );
        m_elementsView[er][esr] = m_elementsData[er][esr];
      }
    }
  }


  localIndex_array m_nodes;
  localIndex_array m_edges;
  localIndex_array m_faces;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > m_elements;
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > m_elementsView;

  array1d< array1d< localIndex_array > > m_elementsData;
  buffer_type::size_type m_size;

};

struct TopologyChangeUnpackStepData : public TopologyChangeStepData
{
  void init( buffer_type const & receiveBuffer,
             ElementRegionManager const & elemManager )
  {
    m_bufferPtr = receiveBuffer.data();
    TopologyChangeStepData::init( elemManager );
  }

  buffer_unit_type const * m_bufferPtr;
};

}

}

#endif /* GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_ */
