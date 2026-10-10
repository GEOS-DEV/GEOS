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


int RLFGraphColoring::colorGraph( const stdVector< size_t > & GEOS_UNUSED_PARAM( adjncy ))
{
  return -1;
}


stdVector< int > RLFGraphColoring::colorGraph( const stdVector< size_t > & xadj, const stdVector< size_t > & adjncy )
{
  return RecursiveLargestFirstColoring( xadj, adjncy );
}


stdVector< int > RLFGraphColoring::RecursiveLargestFirstColoring( const stdVector< size_t > & xadj, const stdVector< size_t > & adjncy )
{
  if( xadj.empty())
  {
    return {};
  }

  stdVector< int > color( xadj.size() - 1, -1 );
  std::unordered_set< size_t > uncoloredNodes;
  for( size_t i = 0; i < (xadj.size() - 1); ++i )
  {
    uncoloredNodes.insert( i );
  }
  int currentColor = 0;

  while( !uncoloredNodes.empty())
  {
    // Find the uncolored node with the largest number of connections
    auto cmp = [&]( size_t a, size_t b ) {
      const auto degreeA = getGraphNodeDegree( a, xadj );
      const auto degreeB = getGraphNodeDegree( b, xadj );
      if( degreeA != degreeB )
      {
        return degreeA < degreeB;
      }
      return a > b;
    };
    size_t startNode = *std::max_element( uncoloredNodes.begin(), uncoloredNodes.end(), cmp );

    // Declare the set of forbidden nodes (nodes adjacent to nodes_to_color)
    std::unordered_set< size_t > forbiddenNodes;

    // Create the candidate set once per color.
    // At this stage, candidates are all uncolored nodes... that set will shrink as we go
    std::unordered_set< size_t > candidatesForThisColor = uncoloredNodes;

    size_t selectedNode = startNode;

    while( selectedNode != std::numeric_limits< size_t >::max())
    {
      // Color the selected node and remove it from future consideration
      color[selectedNode] = currentColor;
      uncoloredNodes.erase( selectedNode );
      candidatesForThisColor.erase( selectedNode );

      // Update the set of forbidden nodes for the current color
      // Any node adjacent to the selectedNode cannot be colored with currentColor
      const auto selectedNodeNeighbors = getGraphNodeNeighbors( selectedNode, xadj, adjncy );
      forbiddenNodes.insert( selectedNodeNeighbors.begin(), selectedNodeNeighbors.end());

      for( const auto & neighbor : selectedNodeNeighbors )
      {
        candidatesForThisColor.erase( neighbor );
      }

      // Find the next best node to color from the remaining candidates
      size_t bestCandidate = std::numeric_limits< size_t >::max();
      size_t maxCommonNeighbors = 0;
      size_t bestCandidateDegree = std::numeric_limits< size_t >::max();

      // Iterate over the shrinking set of candidates
      for( const auto & node : candidatesForThisColor )
      {
        const auto nodeNeighbors = getGraphNodeNeighbors( node, xadj, adjncy );

        // Count common neighbors with the forbidden set
        size_t commonNeighborsCount = 0;
        for( const auto & neighbor : nodeNeighbors )
        {
          if( forbiddenNodes.count( neighbor ))
          {
            ++commonNeighborsCount;
          }
        }

        // Apply RLF heuristic
        bool shouldSelect = false;
        if( bestCandidate == std::numeric_limits< size_t >::max() || commonNeighborsCount > maxCommonNeighbors )
        {
          shouldSelect = true;
        }
        else if( commonNeighborsCount == maxCommonNeighbors )
        {
          const size_t nodeDegree = getGraphNodeDegree( node, xadj );
          if( nodeDegree < bestCandidateDegree )
          {
            shouldSelect = true;
          }
        }

        if( shouldSelect )
        {
          bestCandidate = node;
          maxCommonNeighbors = commonNeighborsCount;
          bestCandidateDegree = getGraphNodeDegree( node, xadj );
        }
      }
      selectedNode = bestCandidate;
    }

    // Move to the next color
    ++currentColor;
  }

  return color;
}


size_t RLFGraphColoring::getNumberOfColors( const stdVector< int > & colors ) const
{
  return GraphColoringBase::getNumberOfColors( colors );
}

bool RLFGraphColoring::isColoringValid( const stdVector< size_t > & xadj, const stdVector< size_t > & adjncy, const stdVector< int > & colors ) const
{
  return GraphColoringBase::isColoringValid( xadj, adjncy, colors );
}

} // namespace graph
} // namespace geos
