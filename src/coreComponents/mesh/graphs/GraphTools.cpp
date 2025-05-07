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
 * @file GraphTools.cpp
 *
 */

#include "GraphTools.hpp"
#include <cstdlib>
#include <ctime>
#include <algorithm>


namespace geos
{
namespace graph
{

size_t getGraphNodeDegree( idx_t node, const std::vector< idx_t > & xadj )
{
  return xadj[node + 1] - xadj[node];
}


std::unordered_set< idx_t > getGraphNodeNeighbors( idx_t node, const std::vector< idx_t > & xadj, const std::vector< idx_t > & adjncy )
{
  std::unordered_set< idx_t > neighbors;
  for( idx_t i = xadj[node]; i < xadj[node + 1]; ++i )
  {
    neighbors.insert( adjncy[i] );
  }
  return neighbors;
}


bool isGraphValid( const std::vector< idx_t > & xadj,
                   const std::vector< idx_t > & adjncy )
{
  idx_t num_nodes = xadj.size() - 1;

  for( idx_t i = 0; i < num_nodes; ++i )
  {
    std::unordered_set< idx_t > neighbors;
    for( idx_t j = xadj[i]; j < xadj[i + 1]; ++j )
    {
      idx_t neighbor = adjncy[j];

      // Check for out-of-bounds indices
      if( neighbor >= num_nodes )
      {
        //std::cerr << "Neighbor index " << neighbor << " is out of bounds for node " << i << "." << std::endl;
        return false;
      }

      // Check for duplicate nodes
      if( neighbors.find( neighbor ) != neighbors.end())
      {
        //std::cerr << "Duplicate neighbor " << neighbor << " found for node " << i << "." << std::endl;
        return false;
      }
      neighbors.insert( neighbor );

      // Check for bidirectional connection
      bool bidirectional = false;
      for( idx_t k = xadj[neighbor]; k < xadj[neighbor + 1]; ++k )
      {
        if( adjncy[k] == i )
        {
          bidirectional = true;
          break;
        }
      }
      if( !bidirectional )
      {
        //std::cerr << "Node " << neighbor << " does not have a bidirectional connection with node " << i << "." << std::endl;
        return false;
      }
    }
  }

  return true;
}


std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphRandom( size_t numVertices, size_t numEdges )
{
  std::vector< std::pair< size_t, size_t > > edges;
  srand( static_cast< unsigned int >(time( 0 )));

  // Generate random edges
  for( size_t i = 0; i < numEdges; ++i )
  {
    size_t u = rand() % numVertices;
    size_t v = rand() % numVertices;
    if( u != v && std::find( edges.begin(), edges.end(), std::make_pair( v, u )) == edges.end())
    {
      edges.push_back( {u, v} );
      edges.push_back( {v, u} );
    }
  }

  // Sort edges by the first vertex
  std::sort( edges.begin(), edges.end());

  // Initialize xadj and adjncy
  std::vector< idx_t > xadj( numVertices + 1, 0 );
  std::vector< idx_t > adjncy;
  adjncy.reserve( edges.size());

  // Fill xadj and adjncy
  size_t currentVertex = 0;
  for( const auto & edge : edges )
  {
    while( currentVertex <= edge.first )
    {
      xadj[currentVertex] = adjncy.size();
      ++currentVertex;
    }
    adjncy.push_back( edge.second );
  }
  while( currentVertex <= numVertices )
  {
    xadj[currentVertex] = adjncy.size();
    ++currentVertex;
  }

  //isGraphValid( xadj, adjncy );
  return {xadj, adjncy};
}


std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphCartPartitionning3D( idx_t nx, idx_t ny, idx_t nz, const std::vector< std::array< int, 3 > > & neighbor_offsets )
{
  idx_t num_nodes = nx * ny * nz;
  std::vector< idx_t > xadj( num_nodes + 1, 0 );
  std::vector< idx_t > adjncy;

  auto getNodeIndex = [nx, ny] ( const idx_t x, const idx_t y, const idx_t z )
  {
    return x + nx * (y + ny * z);
  };

  idx_t node_counter = 0;
  for( idx_t z = 0; z < nz; ++z )
  {
    for( idx_t y = 0; y < ny; ++y )
    {
      for( idx_t x = 0; x < nx; ++x )
      {
        std::vector< idx_t > neighbors;
        for( const auto & offset : neighbor_offsets )
        {

          int nx_pos = static_cast< int >(x) + offset[0];
          int ny_pos = static_cast< int >(y) + offset[1];
          int nz_pos = static_cast< int >(z) + offset[2];

          if( nx_pos >= 0 && nx_pos < static_cast< int >(nx) && ny_pos >= 0
              && ny_pos < static_cast< int >(ny) && nz_pos >= 0
              && nz_pos < static_cast< int >(nz))
          {
            neighbors.push_back( getNodeIndex( nx_pos, ny_pos, nz_pos ));
          }
        }

        xadj[node_counter + 1] = xadj[node_counter] + neighbors.size();
        adjncy.insert( adjncy.end(), neighbors.begin(), neighbors.end());
        ++node_counter;
      }
    }
  }

  return {xadj, adjncy};
}

std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphCartPartitionning3D6( idx_t nx, idx_t ny, idx_t nz )
{
  std::vector< std::array< int, 3 > > neighbor_offsets = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}
  };
  return generateGraphCartPartitionning3D( nx, ny, nz, neighbor_offsets );
}

std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphCartPartitionning3D26( idx_t nx, idx_t ny, idx_t nz )
{
  std::vector< std::array< int, 3 > > neighbor_offsets;
  for( int dz = -1; dz <= 1; ++dz )
  {
    for( int dy = -1; dy <= 1; ++dy )
    {
      for( int dx = -1; dx <= 1; ++dx )
      {
        if( dx != 0 || dy != 0 || dz != 0 )
        {
          neighbor_offsets.push_back( {dx, dy, dz} );
        }
      }
    }
  }
  return generateGraphCartPartitionning3D( nx, ny, nz, neighbor_offsets );
}


} // namespace graph
} // namespace geos
