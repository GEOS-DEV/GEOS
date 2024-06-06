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

#include "NewGlobalNumbering.hpp"

//#include "common/MpiWrapper.hpp"
#include "codingUtilities/Utilities.hpp"

#include "LvArray/src/ArrayOfArraysView.hpp"

namespace geos::ghosting
{

struct GlobalToLocal
{
  std::map< NodeGlbIdx, NodeLocIdx > nodes;
  std::map< EdgeGlbIdx, EdgeLocIdx > edges;
  std::map< FaceGlbIdx, FaceLocIdx > faces;
  std::map< CellGlbIdx, CellLocIdx > cells;
};

struct LocalToGlobal
{
  std::map< NodeLocIdx, NodeGlbIdx > nodes;
  std::map< EdgeLocIdx, EdgeGlbIdx > edges;  // TODO Do we want to make this a vector already? It's surely the most efficient way to build it.
  std::map< FaceLocIdx, FaceGlbIdx > faces;
  std::map< CellLocIdx, CellGlbIdx > cells;
};

template< class GI, class LI >
void buildL2GMappings( std::set< GI > const & gis,
                       std::map< GI, LI > & g2l,
                       std::map< LI, GI > & l2g ) // TODO we can make it a vector -> simple push_backs will be enough to build it
{
  g2l.clear();
  l2g.clear();

  LI li{ 0 };
  for( GI const & gi: gis )
  {
    g2l.insert( { gi, li } );
    l2g.insert( { li, gi } );
    ++li;
  }

  // TODO add a check to see if we have the full range of local indices.
}

std::tuple< GlobalToLocal, LocalToGlobal > buildL2GMappings( MeshGraph const & graph )
{
  GlobalToLocal g2l;
  LocalToGlobal l2g;

  // Use `std::ranges::views::keys` when switching to C++20
  auto const keys = []( auto const & m )
  {
    return mapKeys< std::set >( m );
  };

  buildL2GMappings( keys( graph.n2pos ), g2l.nodes, l2g.nodes );
  buildL2GMappings( keys( graph.e2n ), g2l.edges, l2g.edges );
  buildL2GMappings( keys( graph.f2e ), g2l.faces, l2g.faces );
  buildL2GMappings( keys( graph.c2f ), g2l.cells, l2g.cells );

  return { g2l, l2g };
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


MeshGraph mergeMeshGraph( MeshGraph const & owned,
                          MeshGraph const & present,
                          MeshGraph const & ghosts )
{
  MeshGraph result{ owned };
  for( MeshGraph const & graph: { present, ghosts } )
  {
    result.c2f.insert( std::cbegin( graph.c2f ), std::cend( graph.c2f ) );
    result.f2e.insert( std::cbegin( graph.f2e ), std::cend( graph.f2e ) );
    result.e2n.insert( std::cbegin( graph.e2n ), std::cend( graph.e2n ) );
    result.n2pos.insert( std::cbegin( graph.n2pos ), std::cend( graph.n2pos ) );
  }
  return result;
}

struct DownwardMappings
{
  std::map< EdgeLocIdx, std::tuple< NodeLocIdx, NodeLocIdx > > e2n;
  std::map< FaceLocIdx, std::vector< NodeLocIdx > > f2n;
  std::map< FaceLocIdx, std::vector< EdgeLocIdx > > f2e;
  std::map< CellLocIdx, std::vector< NodeLocIdx > > c2n;
  std::map< CellLocIdx, std::vector< FaceLocIdx > > c2f;
  std::map< CellLocIdx, std::vector< EdgeLocIdx > > c2e;
};


struct UpwardMappings
{
  std::map< EdgeLocIdx, std::vector< FaceLocIdx > > e2f;
  std::map< FaceLocIdx, std::vector< CellLocIdx > > f2c;
  std::map< NodeLocIdx, std::vector< EdgeLocIdx > > n2e;
  std::map< NodeLocIdx, std::vector< FaceLocIdx > > n2f;
  std::map< NodeLocIdx, std::vector< CellLocIdx > > n2c;
};

template< class T, class U >
ArrayOfArrays< localIndex > convertToAoA( std::map< T, std::vector< U > > const & t2u )
{
  static_assert( !std::is_same_v< T, U > );
  ArrayOfArrays< localIndex > t2u_;

  std::size_t const numTs = std::size( t2u );

  std::vector< int > sizes;
  sizes.reserve( numTs );
  for( auto const & [_, us]: t2u )
  {
    sizes.emplace_back( std::size( us ) );
  }
  t2u_.resizeFromCapacities< serialPolicy >( numTs, sizes.data() );
  for( auto const & [t, us]: t2u )
  {
    for( U const & u: us )
    {
      t2u_.emplaceBack( t.get(), u.get() );
    }
  }

  return t2u_;
}

template< class T, class U >
array2d< localIndex > convertToA2d( std::map< T, std::vector< U > > const & t2u,
                                    int dimU )
{
  static_assert( !std::is_same_v< T, U > );
  array2d< localIndex > t2u_( std::size( t2u ), dimU );
  t2u_.setValues< serialPolicy >( -1 );

  for( auto const & [t, us]: t2u )
  {
    for( int i = 0; i < intConv< int >( std::size( us ) ); ++i )
    {
      t2u_( t.get(), i ) = us[i].get();
    }
  }
  return t2u_;
}

template< class GI, class LI >
std::tuple< array1d< globalIndex >, unordered_map< globalIndex, localIndex > >
convertGlbLoc( std::map< GI, LI > const & g2l )
{
  static_assert( !std::is_same_v< GI, LI > );
  array1d< globalIndex > l2g_( std::size( g2l ) );
  unordered_map< globalIndex, localIndex > g2l_;
  for( auto const & [gi, li]: g2l )
  {
    l2g_[li.get()] = gi.get();
    g2l_[gi.get()] = li.get();
  }

  return { l2g_, g2l_ };
}

EdgeMgrImpl makeFlavorlessEdgeMgrImpl( std::size_t const & numEdges,
                                       std::vector< integer > const & ghostRank,
                                       std::map< EdgeLocIdx, std::tuple< NodeLocIdx, NodeLocIdx > > const & e2n,
                                       std::map< EdgeLocIdx, std::vector< FaceLocIdx > > const & e2f,
                                       std::map< EdgeGlbIdx, EdgeLocIdx > const & eg2l,
                                       std::map< EdgeLocIdx, EdgeGlbIdx > const & el2g )
{
  GEOS_ASSERT_EQ( numEdges, std::size( ghostRank ) );
  GEOS_ASSERT_EQ( numEdges, std::size( e2n ) );
  GEOS_ASSERT_EQ( numEdges, std::size( e2f ) );
  GEOS_ASSERT_EQ( numEdges, std::size( eg2l ) );
  GEOS_ASSERT_EQ( numEdges, std::size( el2g ) );

  array2d< localIndex > e2n_( numEdges, 2 );
  for( int i = 0; i < intConv< int >( numEdges ); ++i )
  {
    std::tuple< NodeLocIdx, NodeLocIdx > nlis = e2n.at( EdgeLocIdx{ i } );
    e2n_[i][0] = std::get< 0 >( nlis ).get();
    e2n_[i][1] = std::get< 1 >( nlis ).get();
  }

  auto [l2g, g2l] = convertGlbLoc( eg2l );

  array1d< integer > ghostRank_( numEdges );
  for( integer const & gr: ghostRank )
  {
    ghostRank_.emplace_back( gr );
  }

  return EdgeMgrImpl( numEdges, std::move( ghostRank_ ), std::move( e2n_ ), convertToAoA( e2f ), std::move( g2l ), std::move( l2g ) );
}

FaceMgrImpl makeFlavorlessFaceMgrImpl( std::size_t const & numFaces,
                                       std::vector< integer > const & ghostRank,
                                       std::map< FaceLocIdx, std::vector< NodeLocIdx > > const & f2n,
                                       std::map< FaceLocIdx, std::vector< EdgeLocIdx > > const & f2e,
                                       std::map< FaceLocIdx, std::vector< CellLocIdx > > const & f2c,
                                       std::map< FaceGlbIdx, FaceLocIdx > const & fg2l,
                                       std::map< FaceLocIdx, FaceGlbIdx > const & fl2g )
{
  GEOS_ASSERT_EQ( numFaces, std::size( ghostRank ) );
  GEOS_ASSERT_EQ( numFaces, std::size( f2n ) );
  GEOS_ASSERT_EQ( numFaces, std::size( f2e ) );
  GEOS_ASSERT_EQ( numFaces, std::size( f2c ) );
  GEOS_ASSERT_EQ( numFaces, std::size( fg2l ) );
  GEOS_ASSERT_EQ( numFaces, std::size( fl2g ) );

  auto [l2g, g2l] = convertGlbLoc( fg2l );

  array1d< integer > ghostRank_( numFaces );
  for( integer const & gr: ghostRank )
  {
    ghostRank_.emplace_back( gr );
  }

  return FaceMgrImpl( numFaces, std::move( ghostRank_ ), convertToAoA( f2n ), convertToAoA( f2e ), convertToA2d( f2c, 2 ), std::move( g2l ), std::move( l2g ) );
}

NodeMgrImpl makeFlavorlessNodeMgrImpl( std::size_t const & numNodes,
                                       std::vector< integer > const & ghostRank,
                                       std::map< NodeGlbIdx, std::array< double, 3 > > const & n2pos,
                                       std::map< NodeLocIdx, std::vector< EdgeLocIdx > > const & n2e,
                                       std::map< NodeLocIdx, std::vector< FaceLocIdx > > const & n2f,
                                       std::map< NodeLocIdx, std::vector< CellLocIdx > > const & n2c,
                                       std::map< NodeGlbIdx, NodeLocIdx > const & ng2l,
                                       std::map< NodeLocIdx, NodeGlbIdx > const & nl2g )
{
  GEOS_ASSERT_EQ( numNodes, std::size( ghostRank ) );
  GEOS_ASSERT_EQ( numNodes, std::size( n2pos ) );
  GEOS_ASSERT_EQ( numNodes, std::size( n2e ) );
  GEOS_ASSERT_EQ( numNodes, std::size( n2f ) );
  GEOS_ASSERT_EQ( numNodes, std::size( n2c ) );
  GEOS_ASSERT_EQ( numNodes, std::size( ng2l ) );
  GEOS_ASSERT_EQ( numNodes, std::size( nl2g ) );

  // TODO MISSING cell type. Should be OK, the information is conveyed.
  // TODO MISSING get the cell -> numNodesPerElement... from the original CellBlock
  auto [l2g, g2l] = convertGlbLoc( ng2l );

  array1d< integer > ghostRank_( numNodes );
  for( integer const & gr: ghostRank )
  {
    ghostRank_.emplace_back( gr );
  }

  array2d< real64, nodes::REFERENCE_POSITION_PERM > positions( numNodes, 3 );
  for( auto const & [ngi, pos]: n2pos )
  {
    NodeLocIdx const nli = ng2l.at( ngi );
    positions( nli.get(), 0 ) = pos[0];
    positions( nli.get(), 1 ) = pos[2];
    positions( nli.get(), 2 ) = pos[2];
  }

  return NodeMgrImpl( intConv< localIndex >( numNodes ), std::move( positions ), std::move( ghostRank_ ), convertToAoA( n2e ), convertToAoA( n2f ), convertToAoA( n2c ), std::move( l2g ) );
}

CellBlkImpl makeFlavorlessCellBlkImpl( std::size_t const & numCells,
                                       std::vector< integer > const & ghostRank,
                                       std::map< CellLocIdx, std::vector< NodeLocIdx > > const & c2n,
                                       std::map< CellLocIdx, std::vector< EdgeLocIdx > > const & c2e,
                                       std::map< CellLocIdx, std::vector< FaceLocIdx > > const & c2f,
                                       std::map< CellGlbIdx, CellLocIdx > const & cg2l,
                                       std::map< CellLocIdx, CellGlbIdx > const & cl2g )
{
  GEOS_ASSERT_EQ( numCells, std::size( ghostRank ) );
  GEOS_ASSERT_EQ( numCells, std::size( c2n ) );
  GEOS_ASSERT_EQ( numCells, std::size( c2e ) );
  GEOS_ASSERT_EQ( numCells, std::size( c2f ) );
  GEOS_ASSERT_EQ( numCells, std::size( cg2l ) );
  GEOS_ASSERT_EQ( numCells, std::size( cl2g ) );

  // TODO MISSING cell type. Should be OK, the information is conveyed.
  // TODO MISSING get the cell -> numNodesPerElement... from the original CellBlock
  auto [l2g, g2l] = convertGlbLoc( cg2l );

  array1d< integer > ghostRank_( numCells );
  for( integer const & gr: ghostRank )
  {
    ghostRank_.emplace_back( gr );
  }

  return CellBlkImpl( intConv< localIndex >( numCells ), ghostRank_, convertToA2d( c2n, 8 ), convertToA2d( c2e, 12 ), convertToA2d( c2f, 6 ), std::move( l2g ) );
}

DownwardMappings buildDownwardMappings( GlobalToLocal const & g2l,
                                        MeshGraph const & graph )
{
  DownwardMappings res;

  // Building the `e2n` (edges to nodes) mapping
  for( auto const & [egi, ngis]: graph.e2n )
  {
    NodeLocIdx const nli0 = g2l.nodes.at( std::get< 0 >( ngis ) );
    NodeLocIdx const nli1 = g2l.nodes.at( std::get< 1 >( ngis ) );
    res.e2n[g2l.edges.at( egi )] = { nli0, nli1 };
  }

  // Building the `f2n` (face to nodes) and `f2e` (faces to edges) mappings
  for( auto const & [fgi, edgeInfos]: graph.f2e )
  {
    FaceLocIdx const & fli = g2l.faces.at( fgi );

    std::vector< NodeLocIdx > & nodes = res.f2n[fli];
    std::vector< EdgeLocIdx > & edges = res.f2e[fli];

    nodes.reserve( std::size( edgeInfos ) );
    edges.reserve( std::size( edgeInfos ) );

    for( EdgeInfo const & edgeInfo: edgeInfos )
    {
      std::tuple< NodeGlbIdx, NodeGlbIdx > const & ngis = graph.e2n.at( edgeInfo.index );
      nodes.emplace_back( edgeInfo.start == 0 ? g2l.nodes.at( std::get< 0 >( ngis ) ) : g2l.nodes.at( std::get< 1 >( ngis ) ) );

      edges.emplace_back( g2l.edges.at( edgeInfo.index ) );
    }
  }

  // Building the `c2n` (cell to nodes), `c2e` (cell to edges) and `c2f` (cell to faces) mappings
  for( auto const & [cgi, faceInfos]: graph.c2f )
  {
    CellLocIdx const & cli = g2l.cells.at( cgi );

    std::vector< FaceLocIdx > & faces = res.c2f[cli];
    std::vector< EdgeLocIdx > & edges = res.c2e[cli];

    faces.reserve( std::size( faceInfos ) );
    edges.reserve( std::size( faceInfos ) );

    // c2f
    for( FaceInfo const & faceInfo: faceInfos )
    {
      faces.emplace_back( g2l.faces.at( faceInfo.index ) );
    }

    // c2e
    std::set< EdgeLocIdx > tmpEdges;
    for( FaceLocIdx const & fli: faces )
    {
      std::vector< EdgeLocIdx > const & es = res.f2e.at( fli );
      tmpEdges.insert( std::cbegin( es ), std::cend( es ) );
    }
    edges.assign( std::cbegin( tmpEdges ), std::cend( tmpEdges ) );

    // c2n
    FaceInfo const & bottomFace = faceInfos.at( 4 ); // (0, 3, 2, 1) // TODO depends on element type.
    FaceInfo const & topFace = faceInfos.at( 5 ); // (4, 5, 6, 7)

    std::vector< NodeLocIdx > const & bottomNodes = res.f2n.at( g2l.faces.at( bottomFace.index ) );
    std::vector< NodeLocIdx > const & topNodes = res.f2n.at( g2l.faces.at( topFace.index ) );

    std::vector< NodeLocIdx > const bn = resetFaceNodes( bottomNodes, bottomFace.isFlipped, bottomFace.start );
    std::vector< NodeLocIdx > const tn = resetFaceNodes( topNodes, topFace.isFlipped, topFace.start );

    // TODO carefully check the ordering...
    res.c2n[cli] = { bn[0], bn[3], bn[2], bn[1], tn[0], tn[1], tn[2], tn[3] };
//    nodes.insert( std::end( nodes ), std::cbegin( bn ), std::cend( bn ) );
//    nodes.insert( std::end( nodes ), std::cbegin( tn ), std::cend( tn ) );
  }

  return res;
}

/**
 * @brief Simple inversions of the downward mappings.
 * @param downwardMappings The downward mappings.
 * @return The mappings.
 */
UpwardMappings buildUpwardMappings( DownwardMappings const & downwardMappings )
{
  UpwardMappings res;

  // Filling `e2f`
  for( auto const & [fli, elis]: downwardMappings.f2e )
  {
    for( EdgeLocIdx const & eli: elis )
    {
      res.e2f[eli].emplace_back( fli );
    }
  }

  // Filling `f2c`
  for( auto const & [cli, flis]: downwardMappings.c2f )
  {
    for( FaceLocIdx const & fli: flis )
    {
      res.f2c[fli].emplace_back( cli );
    }
  }

  // Filling `n2e`
  for( auto const & [eli, nlis]: downwardMappings.e2n )
  {
    res.n2e[std::get< 0 >( nlis )].emplace_back( eli );
    res.n2e[std::get< 1 >( nlis )].emplace_back( eli );
  }

  // Filling `n2f`
  for( auto const & [fli, nlis]: downwardMappings.f2n )
  {
    for( NodeLocIdx const & nli: nlis )
    {
      res.n2f[nli].emplace_back( fli );
    }
  }

  // Filling `n2c`
  for( auto const & [cli, nlis]: downwardMappings.c2n )
  {
    for( NodeLocIdx const & nli: nlis )
    {
      res.n2c[nli].emplace_back( cli );
    }
  }

  return res;
}

void buildPods( MeshGraph const & owned,
                MeshGraph const & present,
                MeshGraph const & ghosts,
                GhostRecv const & recv,
                GhostSend const & send,
                MeshMappingImpl & meshMappings )
{
  MeshGraph const graph = mergeMeshGraph( owned, present, ghosts );

  auto const [g2l, l2g] = buildL2GMappings( graph );

  DownwardMappings const downwardMappings = buildDownwardMappings( g2l, graph );
  UpwardMappings const upwardMappings = buildUpwardMappings( downwardMappings );

  NodeMgrImpl const nodeMgr = makeFlavorlessNodeMgrImpl( std::size( g2l.nodes ),
                                                         buildGhostRank( g2l.nodes, send.nodes, recv.nodes ),
                                                         graph.n2pos,
                                                         upwardMappings.n2e,
                                                         upwardMappings.n2f,
                                                         upwardMappings.n2c,
                                                         g2l.nodes,
                                                         l2g.nodes );

  EdgeMgrImpl const edgeMgr = makeFlavorlessEdgeMgrImpl( std::size( g2l.edges ),
                                                         buildGhostRank( g2l.edges, send.edges, recv.edges ),
                                                         downwardMappings.e2n,
                                                         upwardMappings.e2f,
                                                         g2l.edges,
                                                         l2g.edges );

  FaceMgrImpl const faceMgr = makeFlavorlessFaceMgrImpl( std::size( g2l.faces ),
                                                         buildGhostRank( g2l.faces, send.faces, recv.faces ),
                                                         downwardMappings.f2n,
                                                         downwardMappings.f2e,
                                                         upwardMappings.f2c,
                                                         g2l.faces,
                                                         l2g.faces );

  CellBlkImpl const cellBlock = makeFlavorlessCellBlkImpl( std::size( g2l.cells ),
                                                           buildGhostRank( g2l.cells, send.cells, recv.cells ),
                                                           downwardMappings.c2n,
                                                           downwardMappings.c2e,
                                                           downwardMappings.c2f,
                                                           g2l.cells,
                                                           l2g.cells );

  CellMgrImpl const cellMgr( cellBlock );

  meshMappings.setCellMgr( cellMgr );
  meshMappings.setEdgeMgr( edgeMgr );
  meshMappings.setFaceMgr( faceMgr );
  meshMappings.setNodeMgr( nodeMgr );
}

}
