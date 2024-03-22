/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "MeshBlockMappings.hpp"

namespace geos
{

using namespace dataRepository;

MeshBlockMappings::MeshBlockMappings( string const & name,
                                      Group * const parent )
  :
  CellBlockManagerABC( name, parent )
{
  this->registerGroup< Group >( viewKeyStruct::cellBlocks() );
  this->registerGroup< Group >( viewKeyStruct::faceBlocks() );
  this->registerGroup< Group >( viewKeyStruct::lineBlocks() );
}

Group & MeshBlockMappings::getCellBlocks()
{
  return this->getGroup( viewKeyStruct::cellBlocks() );
}

Group & MeshBlockMappings::getFaceBlocks()
{
  return this->getGroup( viewKeyStruct::faceBlocks() );
}

LineBlockABC const & MeshBlockMappings::getLineBlock( string name ) const
{
  return this->getGroup( viewKeyStruct::lineBlocks() ).getGroup< LineBlockABC >( name );
}

Group const & MeshBlockMappings::getCellBlocks() const
{
  return this->getGroup( viewKeyStruct::cellBlocks() );
}

Group const & MeshBlockMappings::getFaceBlocks() const
{
  return this->getGroup( viewKeyStruct::faceBlocks() );
}

localIndex MeshBlockMappings::numNodes() const
{
  return m_numNodes;
}

localIndex MeshBlockMappings::numEdges() const
{
  return m_numEdges;
}

localIndex MeshBlockMappings::numFaces() const
{
  return m_numFaces;
}

array2d< real64, nodes::REFERENCE_POSITION_PERM > MeshBlockMappings::getNodePositions() const
{
  return {};
}

ArrayOfArrays< localIndex > MeshBlockMappings::getNodeToEdges() const
{
  return {};
}

ArrayOfArrays< localIndex > MeshBlockMappings::getNodeToFaces() const
{
  return {};
}

ToCellRelation< ArrayOfArrays< localIndex>> MeshBlockMappings::getNodeToElements() const
{
  return {};
}

array2d< localIndex > MeshBlockMappings::getEdgeToNodes() const
{
  return {};
}

ArrayOfArrays< localIndex > MeshBlockMappings::getEdgeToFaces() const
{
  return {};
}

ArrayOfArrays< localIndex > MeshBlockMappings::getFaceToNodes() const
{
  return {};
}

ArrayOfArrays< localIndex > MeshBlockMappings::getFaceToEdges() const
{
  return {};
}

ToCellRelation< array2d< localIndex>> MeshBlockMappings::getFaceToElements() const
{
  return {};
}

array1d< globalIndex > MeshBlockMappings::getNodeLocalToGlobal() const
{
  return {};
}

std::map< string, SortedArray< localIndex>> const & MeshBlockMappings::getNodeSets() const
{
  return m_nodeSets;
}

real64 MeshBlockMappings::getGlobalLength() const
{
  return m_globalLength;
}

void MeshBlockMappings::generateHighOrderMaps( localIndex const order,
                                               globalIndex const maxVertexGlobalID,
                                               globalIndex const maxEdgeGlobalID,
                                               globalIndex const maxFaceGlobalID,
                                               arrayView1d< globalIndex const > const edgeLocalToGlobal,
                                               arrayView1d< globalIndex const > const faceLocalToGlobal )
{
  GEOS_ERROR( "`MeshBlockMappings::generateHighOrderMaps` is not implemented and should not be called." );
}

} // geos