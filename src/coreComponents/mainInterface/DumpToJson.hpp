/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_DUMPTOJSON_HPP
#define GEOSX_DUMPTOJSON_HPP

#include "mesh/generators/CellBlockManagerABC.hpp"

#include <nlohmann/json.hpp>

namespace geosx
{

using json = nlohmann::json;

std::vector< std::vector< localIndex > > convert( ArrayOfArrays< localIndex > const & v )
{
  std::vector< std::vector< localIndex > > result;

  for( localIndex i = 0; i < v.size(); ++i )
  {
    auto const & vv = v[i];
    std::vector< localIndex > tmp;
    for( localIndex ii = 0; ii < vv.size(); ++ii )
    {
      tmp.push_back( vv[ii] );
    }
    result.push_back( tmp );
  }

  return result;
}

std::vector< std::vector< real64 > > convert( array2d< real64, nodes::REFERENCE_POSITION_PERM > const & v )
{
  std::vector< std::vector< real64 > > result;

  for( localIndex i = 0; i < v.size( 0 ); ++i )
  {
    auto const & vv = v[i];
    std::vector< real64 > tmp;
    for( localIndex ii = 0; ii < vv.size(); ++ii )
    {
      tmp.push_back( vv[ii] );
    }
    result.push_back( tmp );
  }

  return result;
}

std::vector< std::vector< localIndex > > convert( array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & v )
{
  std::vector< std::vector< localIndex > > result;

  for( localIndex i = 0; i < v.size( 0 ); ++i )
  {
    auto const & vv = v[i];
    std::vector< localIndex > tmp;
    for( localIndex ii = 0; ii < vv.size(); ++ii )
    {
      tmp.push_back( vv[ii] );
    }
    result.push_back( tmp );
  }

  return result;
}

//std::vector< std::vector< localIndex > > convert( array2d< localIndex > const & v )
//{
//  std::vector< std::vector< localIndex > > result;
//
//  for( localIndex i = 0; i < v.size(0); ++i )
//  {
//    auto const & vv = v[i];
//    std::vector< localIndex > tmp;
//    for( localIndex ii = 0; ii < vv.size(); ++ii )
//    {
//      tmp.push_back( vv[ii] );
//    }
//    result.push_back( tmp );
//  }
//
//  return result;
//}

std::vector< localIndex > convert( SortedArray< localIndex > const & v )
{
  std::vector< localIndex > result;
  for( localIndex i = 0; i < v.size(); ++i )
  {
    result.push_back( v[i] );
  }
  return result;
}

void to_json( json & j,
              const CellBlockABC & cb )
{
  j["num_nodes_per_elem"] = cb.numNodesPerElement();
  j["num_edges_per_elem"] = cb.numEdgesPerElement();
  j["num_faces_per_elem"] = cb.numFacesPerElement();
  j["num_elems"] = cb.numElements();
  j["elem_to_nodes"] = convert( cb.getElemToNodes() );
  j["elem_to_edges"] = convert( cb.getElemToEdges() );
  j["elem_to_faces"] = convert( cb.getElemToFaces() );
}

void to_json( json & j,
              const CellBlockManagerABC & cbm )
{
  std::cout << "json" << std::endl;
  j["num_nodes"] = cbm.numNodes();
  j["num_edges"] = cbm.numEdges();
  j["num_faces"] = cbm.numFaces();
  j["node_positions"] = convert( cbm.getNodePositions() );
  j["node_to_edges"] = convert( cbm.getNodeToEdges() );
  j["node_to_faces"] = convert( cbm.getNodeToFaces() );
  auto const n2e = cbm.getNodeToElements();
  j["node_to_elems"]["block_index"] = convert( n2e.toBlockIndex );
  j["node_to_elems"]["cell_index"] = convert( n2e.toCellIndex );
  j["edge_to_nodes"] = convert( cbm.getEdgeToNodes() );
  j["edge_to_faces"] = convert( cbm.getEdgeToFaces() );
  j["face_to_nodes"] = convert( cbm.getFaceToNodes() );
  j["face_to_edges"] = convert( cbm.getFaceToEdges() );
  auto const f2e = cbm.getFaceToElements();
  j["face_to_elems"]["block_index"] = convert( f2e.toBlockIndex );
  j["face_to_elems"]["cell_index"] = convert( f2e.toCellIndex );

  for( auto const & kv: cbm.getNodeSets() )
  {
    j["node_sets"][kv.first] = convert( kv.second );
  }

  auto const & cbs = cbm.getCellBlocks();
  for( localIndex i = 0; i < cbs.numSubGroups(); ++i )
  {
    j["cell_blocks"].push_back( cbs.getGroup< CellBlockABC >( i ) );
  }

}

}

#endif //GEOSX_DUMPTOJSON_HPP
