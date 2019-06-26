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
 * @file FaceElementRegion.cpp
 */

#include "FaceElementRegion.hpp"

namespace geosx
{
using namespace dataRepository;

FaceElementRegion::FaceElementRegion( string const & name, ManagedGroup * const parent ):
  ElementRegion( name, parent ),
  m_edgesToFractureConnectors(),
  m_fractureConnectorsToEdges(),
  m_fractureConnectorsToFaceElements(),
  m_faceElementsToCells()
{
  this->GetGroup(viewKeyStruct::elementSubRegions)->RegisterGroup<FaceElementSubRegion>("default");

  RegisterViewWrapper( viewKeyStruct::edgesTofractureConnectorsString, &m_edgesToFractureConnectors, 0 )->
    setRestartFlags( RestartFlags::NO_WRITE)->
    setDescription( "A map of edge local indices to the fracture connector local indices.")->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::fractureConnectorsToEdgesString, &m_fractureConnectorsToEdges, 0 )->
    setRestartFlags( RestartFlags::NO_WRITE)->
    setDescription( "A map of fracture connector local indices to edge local indices.")->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::fractureConnectorsToFaceElementsString, &m_fractureConnectorsToFaceElements, 0 )->
    setRestartFlags( RestartFlags::NO_WRITE)->
    setDescription( "A map of fracture connector local indices face element local indices")->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::faceElementsToCellsString, &m_faceElementsToCells, 0 )->
    setRestartFlags( RestartFlags::NO_WRITE)->
    setDescription( "A map of face element local indices to the cell local indices")->
    setSizedFromParent(0);


}

FaceElementRegion::~FaceElementRegion()
{}



localIndex FaceElementRegion::AddToFractureMesh( FaceManager const * const faceManager,
                                                 array1d< array1d<localIndex> > const & originalFaceToEdges,
                                                 string const & subRegionName,
                                                 localIndex const faceIndices[2]  )
{
  localIndex rval = -1;

  set<localIndex> connectedEdges;

  array2d<localIndex > const & faceToElementRegion = faceManager->elementRegionList();
  array2d<localIndex > const & faceToElementSubRegion = faceManager->elementSubRegionList();
  array2d<localIndex > const & faceToElementIndex = faceManager->elementList();

  ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);

  FaceElementSubRegion * subRegion = elementSubRegions->GetGroup<FaceElementSubRegion>(subRegionName);
  subRegion->resize( subRegion->size() + 1 );
  rval = subRegion->size() - 1;

  m_faceElementsToCells.resize( subRegion->size(), 2 );

  FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
  FaceElementSubRegion::EdgeMapType & edgeMap = subRegion->edgeList();
  FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

  OrderedVariableOneToManyRelation const & facesToNodesMap = faceManager->nodeList();
  array1d< array1d<localIndex> > const & facesToEdgesMap = originalFaceToEdges;

  localIndex const kfe = subRegion->size() - 1;

  m_newFractureElements.insert(kfe);

  faceMap[kfe][0] = faceIndices[0];
  faceMap[kfe][1] = faceIndices[1];

  // Add the nodes that compose the new FaceElement to the nodeList
  arrayView1d<localIndex const> const & faceToNodesMap0 = facesToNodesMap[faceIndices[0]];
  arrayView1d<localIndex const> const & faceToNodesMap1 = facesToNodesMap[faceIndices[1]];
  nodeMap[kfe].resize( faceToNodesMap0.size() * 2 );
  for( localIndex a=0 ; a<faceToNodesMap0.size() ; ++a )
  {
    localIndex const aa = a < 2 ? a : faceToNodesMap0.size() - a + 1;
    localIndex const bb = aa == 0 ? aa : faceToNodesMap0.size() - aa;

    // TODO HACK need to generalize to something other than quads
    nodeMap[kfe][a]   = faceToNodesMap0[aa];
    nodeMap[kfe][a+4] = faceToNodesMap1[bb];
  }

  // Add the edges that compose the faceElement to the edge map. This is essentially a copy of
  // the facesToEdges entry.
  arrayView1d<localIndex const> const & faceToEdgesMap = facesToEdgesMap[faceIndices[0]];
  edgeMap[kfe].resize( faceToEdgesMap.size() );
  for( localIndex a=0 ; a<faceToEdgesMap.size() ; ++a )
  {
    edgeMap[kfe][a] = faceToEdgesMap[a];
    connectedEdges.insert( faceToEdgesMap[a] );
  }

  // Add the cell region/subregion/index to the faceElementToCells map

  for( localIndex ke=0 ; ke<2 ; ++ke )
  {
    m_faceElementsToCells.m_toElementRegion[kfe][ke] = faceToElementRegion[faceIndices[ke]][ke];
    m_faceElementsToCells.m_toElementSubRegion[kfe][ke] = faceToElementSubRegion[faceIndices[ke]][ke];
    m_faceElementsToCells.m_toElementIndex[kfe][ke] = faceToElementIndex[faceIndices[ke]][ke];
  }

  // Fill the connectivity between FaceElement entries. This is essentially a copy of the
  // edgesToFaces map, but with differing offsets.
  for( auto const & edge : connectedEdges )
  {
    // check to see if the edgesToFractureConnectors already have an entry
    if( m_edgesToFractureConnectors.count(edge)==0 )
    {
      // if not, then fill increase the size of the fractureConnectors to face elements map and
      // fill the fractureConnectorsToEdges map with the current edge....and the inverse map too.
      m_fractureConnectorsToFaceElements.push_back({});
      m_fractureConnectorsToEdges.push_back(edge);
      m_edgesToFractureConnectors[edge] = m_fractureConnectorsToEdges.size()-1;
    }
    // now fill the fractureConnectorsToFaceElements map. This is analogous to the edge to face map
    localIndex const connectorIndex = m_edgesToFractureConnectors[edge];
    localIndex const numCells = m_fractureConnectorsToFaceElements[connectorIndex].size() + 1;
    m_fractureConnectorsToFaceElements[connectorIndex].resize( numCells );
    m_fractureConnectorsToFaceElements[connectorIndex][ numCells-1 ] = kfe;

    // And fill the list of connectors that will need stencil modifications
    m_recalculateConnectors.insert( connectorIndex );
  }

  subRegion->CalculateElementGeometricQuantities( kfe, faceManager->faceArea() );

  // update the sets
  for( auto const & setIter : faceManager->sets()->wrappers() )
  {
    set<localIndex> const & faceSet = faceManager->sets()->getReference<set<localIndex> >( setIter.first );
    set<localIndex> & faceElementSet = subRegion->sets()->RegisterViewWrapper< set<localIndex> >( setIter.first )->reference();
    for( localIndex a=0 ; a<faceMap.size(0) ; ++a )
    {
      localIndex const faceIndex = faceMap[a][0];
      if( faceSet.count( faceIndex ) )
      {
        faceElementSet.insert( a );
      }
    }
  }

  return rval;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceElementRegion, std::string const &, ManagedGroup * const )

} /* namespace geosx */
