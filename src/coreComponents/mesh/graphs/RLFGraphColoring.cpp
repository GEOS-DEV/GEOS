/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License: LGPL-2.1-only
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
 * @file RLFGraphColoring.cpp
 */

#include "RLFGraphColoring.hpp"
#include "GraphTools.hpp"
#include <algorithm>


namespace geos
{
namespace graph
{

RLFGraphColoring::RLFGraphColoring(): GraphColoringBase()
{}

RLFGraphColoring::~RLFGraphColoring()
{}


int RLFGraphColoring::colorGraph( const std::vector< camp::idx_t > & GEOS_UNUSED_PARAM( adjncy ))
{
  return -1;
}


std::vector< int > RLFGraphColoring::colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy )
{
  return RecursiveLargestFirstColoring( xadj, adjncy );
}

std::vector< int > RLFGraphColoring::RecursiveLargestFirstColoring( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy )
{
  std::vector< int > color( xadj.size() - 1, -1 );
  std::unordered_set< camp::idx_t > uncoloredNodes;
  for( size_t i = 0; i < (xadj.size() - 1); ++i )
  {
    uncoloredNodes.insert( i );
  }
  int currentColor = 0;

  while( !uncoloredNodes.empty())
  {
    // Find the uncolored node with the largest number of connections
    idx_t largestDegreeNode = *std::max_element( uncoloredNodes.begin(), uncoloredNodes.end(),
                                                 [xadj]( const idx_t a, const idx_t b ){ return getGraphNodeDegree( a, xadj ) < getGraphNodeDegree( b, xadj ); } );

    // Save its neighbors
    std::unordered_set< camp::idx_t > largestDegreeNodeNeighbors = getGraphNodeNeighbors( largestDegreeNode, xadj, adjncy );

    // Initialize the node to be colored with the current color
    int selectedNode = static_cast< int >(largestDegreeNode);

    // Declare the set of forbidden nodes (nodes adjacent to nodes_to_color)
    std::unordered_set< camp::idx_t > forbiddenNodes;

    while( selectedNode != -1 )
    {
      // Color the selected node
      color[selectedNode] = currentColor;
      std::unordered_set< camp::idx_t > selectedNodeNeighbors = getGraphNodeNeighbors( selectedNode, xadj, adjncy );
      forbiddenNodes.insert( selectedNodeNeighbors.begin(), selectedNodeNeighbors.end());
      uncoloredNodes.erase( selectedNode );

      // Find new nodes to color with the current color
      std::unordered_set< camp::idx_t > candidateNodesToColor;
      for( camp::idx_t node : uncoloredNodes )
      {
        if( forbiddenNodes.find( node ) == forbiddenNodes.end())
        {
          candidateNodesToColor.insert( node );
        }
      }

      // Among the candidates, select one node
      size_t maxCommonNeighbors = 0;
      selectedNode = -1;
      for( camp::idx_t node : candidateNodesToColor )
      {
        std::unordered_set< camp::idx_t > nodeNeighbors = getGraphNodeNeighbors( node, xadj, adjncy );

        size_t commonNeighbors = 0;
        for( camp::idx_t neighbor : nodeNeighbors )
        {
          if( largestDegreeNodeNeighbors.find( neighbor ) != largestDegreeNodeNeighbors.end())
          {
            ++commonNeighbors;
          }
        }
        size_t nodeDegree = getGraphNodeDegree( node, xadj );

        if( commonNeighbors > maxCommonNeighbors || (commonNeighbors == maxCommonNeighbors && nodeDegree < getGraphNodeDegree( selectedNode, xadj )))
        {
          selectedNode = static_cast< int >(node);
          maxCommonNeighbors = commonNeighbors;
        }
      }
    }

    // Increment the current color
    ++currentColor;
  }

  return color;
}

size_t RLFGraphColoring::getNumberOfColors( const std::vector< int > & colors ) const
{
  return GraphColoringBase::getNumberOfColors( colors );
}

bool RLFGraphColoring::isColoringValid( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy, const std::vector< int > & colors ) const
{
  return GraphColoringBase::isColoringValid( xadj, adjncy, colors );
}


} // namespace graph
} // namespace geos
