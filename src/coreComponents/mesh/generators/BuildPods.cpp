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

#include "BuildPods.hpp"


namespace geos::ghosting
{

void buildPods( MeshGraph const & owned,
                MeshGraph const & present,
                MeshGraph const & ghosts,
                GhostRecv const & recv,
                GhostSend const & send )
{
//  std::size_t const numNodes = std::size( ownerships.nodes );
//  std::size_t const numEdges = std::size( ownerships.edges );
//  std::size_t const numFaces = std::size( ownerships.faces );
//
//  auto [ghostRank, l2g] = buildGhostRankAndL2G( ownerships.edges );
//
//  NodeMgrImpl const nodeMgr( NodeLocIdx{ intConv< localIndex >( numNodes ) } );
//  EdgeMgrImpl const edgeMgr( EdgeLocIdx{ intConv< localIndex >( numEdges ) }, std::move( ghostRank ), std::move( l2g ) );
//  FaceMgrImpl const faceMgr( FaceLocIdx{ intConv< localIndex >( numFaces ) } );
}

}
