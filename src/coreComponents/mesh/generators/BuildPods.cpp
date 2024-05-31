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

#include "Pods.hpp"

//#include "common/MpiWrapper.hpp"
#include "codingUtilities/Utilities.hpp"

#include "LvArray/src/ArrayOfArraysView.hpp"

namespace geos::ghosting
{

struct GlobalNumberings
{
  std::map< NodeGlbIdx, NodeLocIdx > ng2l;
  std::map< NodeLocIdx, NodeGlbIdx > nl2g;

  std::map< EdgeGlbIdx, EdgeLocIdx > eg2l;
  std::map< EdgeLocIdx, EdgeGlbIdx > el2g;  // TODO Do we want to make this a vector already? It's surely the most efficient way to build it.

  std::map< FaceGlbIdx, FaceLocIdx > fg2l;
  std::map< FaceLocIdx, FaceGlbIdx > fl2g;

  std::map< CellGlbIdx, CellLocIdx > cg2l;
  std::map< CellLocIdx, CellGlbIdx > cl2g;
};


template< class GI, class LI >
void buildL2GMappings( std::set< GI > const & owned,
                       std::set< GI > const & present,
                       std::set< GI > const & ghosts,
                       std::map< GI, LI > & g2l,
                       std::map< LI, GI > & l2g ) // TODO we can make it a vector -> simple push_backs will be enough to build it
{
  g2l.clear();
  l2g.clear();

  LI li{ 0 };
  for( std::set< GI > const & gis: { owned, present, ghosts } )
  {
    for( GI const & gi: gis )
    {
      g2l.insert( { gi, li } );
      l2g.insert( { li, gi } );
      ++li;
    }
  }

  // TODO add a check to see if we have the full range of local indices.
}

GlobalNumberings buildL2GMappings( MeshGraph const & owned,
                                   MeshGraph const & present,
                                   MeshGraph const & ghosts )
{
  GlobalNumberings result;

  // Use `std::ranges::views::keys` when switching to C++20
  auto const keys = []( auto const & m )
  {
    return mapKeys< std::set >( m );
  };

  buildL2GMappings( owned.n, present.n, ghosts.n, result.ng2l, result.nl2g );
  buildL2GMappings( keys( owned.e2n ), keys( present.e2n ), keys( ghosts.e2n ), result.eg2l, result.el2g );
  buildL2GMappings( keys( owned.f2e ), keys( present.f2e ), keys( ghosts.f2e ), result.fg2l, result.fl2g );
  buildL2GMappings( keys( owned.c2f ), keys( present.c2f ), keys( ghosts.c2f ), result.cg2l, result.cl2g );

  return result;
}

template< class GI, class LI >
std::vector< integer > buildGhostRank( std::map< GI, LI > const & g2l,
                                       std::map< GI, std::set< MpiRank > > const & send,
                                       std::map< GI, MpiRank > const & recv )
{
  std::vector< integer > ghostRank( std::size( g2l ), -2 );

  for( auto const & [gi, _]: send )
  {
    ghostRank[g2l.at( gi ).get()] = -1;
  }
  for( auto const & [gi, rank]: recv )
  {
    ghostRank[g2l.at( gi ).get()] = rank.get();
  }

  return ghostRank;
}

EdgeMgrImpl makeFlavorLessEdgeMgrImpl( std::size_t const & numEdges,
                                       std::vector< integer > const & ghostRank,
                                       std::map< EdgeLocIdx, std::tuple< NodeLocIdx, NodeLocIdx > > const & e2n,
                                       std::map< EdgeLocIdx, std::vector< FaceLocIdx > > const & e2f,
                                       std::map< EdgeGlbIdx, EdgeLocIdx > const & eg2l,
                                       std::map< EdgeLocIdx, EdgeGlbIdx > const & el2g )
{
  array2d< localIndex > e2n_( numEdges, 2 );
  for( int i = 0; i < intConv< int >( numEdges ); ++i )
  {
    std::tuple< NodeLocIdx, NodeLocIdx > nlis = e2n.at( EdgeLocIdx{ i } );
    e2n_[i][0] = std::get< 0 >( nlis ).get();
    e2n_[i][1] = std::get< 1 >( nlis ).get();
  }

  ArrayOfArrays< localIndex > e2f_;
  std::vector< int > sizes;
  sizes.reserve( numEdges );
  for( auto const & [_, flis]: e2f )
  {
    sizes.emplace_back( std::size( flis ) );
  }
  GEOS_ASSERT_EQ( std::size( sizes ), numEdges );
  e2f_.resizeFromCapacities< serialPolicy >( numEdges, sizes.data() );
  for( auto const & [eli, flis]: e2f )
  {
    for( FaceLocIdx const & fli: flis )
    {
      e2f_.emplaceBack( eli.get(), fli.get() );
    }
  }

  array1d< globalIndex > l2g( numEdges );
  unordered_map< globalIndex, localIndex > g2l;
  for( auto const & [egi, eli]: eg2l )
  {
    l2g[eli.get()] = egi.get();
    g2l[egi.get()] = eli.get();
  }

  array1d< integer > ghostRank_( numEdges );
  for( integer const & gr: ghostRank )
  {
    ghostRank_.emplace_back( gr );
  }

  return EdgeMgrImpl( numEdges, std::move( ghostRank_ ), std::move( e2n_ ), std::move( e2f_ ), std::move( g2l ), std::move( l2g ) );
}

EdgeMgrImpl buildEdgeMgr( GlobalNumberings const & numberings,
                          MeshGraph const & owned,
                          MeshGraph const & present,
                          MeshGraph const & ghosts,
                          GhostRecv const & recv,
                          GhostSend const & send )
{
  // Total number of edges available in the rank (including the ghosted edges).
  std::size_t const numEdges = std::size( numberings.el2g );

  // Building the ghost rank.
  std::vector< integer > const ghostRank = buildGhostRank( numberings.eg2l, send.edges, recv.edges );

  // Building the edges to nodes mapping
  std::map< EdgeLocIdx, std::tuple< NodeLocIdx, NodeLocIdx > > e2n;
  for( MeshGraph const & graph: { owned, present, ghosts } )
  {
    for( auto const & [egi, ngis]: graph.e2n )
    {
      NodeLocIdx const nli0 = numberings.ng2l.at( std::get< 0 >( ngis ) );
      NodeLocIdx const nli1 = numberings.ng2l.at( std::get< 1 >( ngis ) );
      e2n[numberings.eg2l.at( egi )] = { nli0, nli1 };
    }
  }

  // Building the edges to nodes mapping
  std::map< EdgeLocIdx, std::vector< FaceLocIdx > > e2f;
  for( MeshGraph const & graph: { owned, present, ghosts } )
  {
    for( auto const & [fgi, edgeInfos]: graph.f2e )
    {
      for( EdgeInfo const & edgeInfo: edgeInfos )
      {
        e2f[numberings.eg2l.at( edgeInfo.index )].emplace_back( numberings.fg2l.at( fgi ) );
      }
    }
  }

  return makeFlavorLessEdgeMgrImpl( numEdges, ghostRank, e2n, e2f, numberings.eg2l, numberings.el2g );
}

void buildPods( MeshGraph const & owned,
                MeshGraph const & present,
                MeshGraph const & ghosts,
                GhostRecv const & recv,
                GhostSend const & send )
{
  GlobalNumberings const numberings = buildL2GMappings( owned, present, ghosts );
//  GEOS_LOG_RANK( "numberings ng2l = " << json( numberings.ng2l ) );
//  GEOS_LOG_RANK( "owned mesh graph nodes = " << json( owned.n ) );
//  GEOS_LOG_RANK( "present mesh graph nodes = " << json( present.n ) );
//  GEOS_LOG_RANK( "ghosts mesh graph nodes = " << json( ghosts.n ) );
//  GEOS_LOG_RANK( "owned mesh graph e2n = " << json( owned.e2n ) );
//  GEOS_LOG_RANK( "present mesh graph e2n = " << json( present.e2n ) );
//  GEOS_LOG_RANK( "ghosts mesh graph e2n = " << json( ghosts.e2n ) );
//  auto const EdgeMgr = buildEdgeMgr( owned, present, ghosts, recv, send );

//  MpiWrapper::barrier();

  buildEdgeMgr( numberings, owned, present, ghosts, recv, send );
}

}
