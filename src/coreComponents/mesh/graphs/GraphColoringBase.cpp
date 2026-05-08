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
 * @file GraphColoringBase.cpp
 */

#include "GraphColoringBase.hpp"
#include "GraphTools.hpp"
#include <unordered_set>


namespace geos
{
namespace graph
{

bool GraphColoringBase::isColoringValid( const std::vector< size_t > & xadj,
                                         const std::vector< size_t > & adjncy,
                                         const std::vector< int > & coloring )
{
  for( size_t node = 0; node < coloring.size(); ++node )
  {
    int node_color = coloring[node];
    std::unordered_set< size_t > neighbors = getGraphNodeNeighbors( node, xadj, adjncy );
    for( size_t neighbor : neighbors )
    {
      if( coloring[neighbor] == node_color )
      {
        return false;
      }
    }
  }
  return true;
}


size_t GraphColoringBase::getNumberOfColors( const std::vector< int > & colors )
{
  std::unordered_set< int > uniqueColors;
  for( int color : colors )
  {
    if( color >= 0 )
    {
      uniqueColors.insert( color );
    }
  }
  return uniqueColors.size();
}


// Assume only one node per rank.
bool GraphColoringBase::isColoringValid( const std::vector< size_t > & adjncy,
                                         const int color,
                                         MPI_Comm comm )
{
  // Use a collective allgather to retrieve every rank's color.
  // The original point-to-point iSend/iRecv approach with a fixed tag is
  // fragile in Release builds: messages from successive calls to getColor()
  // (one per SurfaceGenerator time-step) can linger in the MPI queue and be
  // matched by the wrong iRecv, producing false color conflicts.
  // Allgather has no tag and is always safe.
  int const size = MpiWrapper::commSize( comm );
  std::vector< int > allColors( size );
  MpiWrapper::allgather( &color, 1, allColors.data(), 1, comm );

  // Check for color conflicts with each listed neighbor.
  bool isLocalColoringValid = true;
  for( camp::idx_t const neighborRank : adjncy )
  {
    if( allColors[neighborRank] == color )
    {
      isLocalColoringValid = false;
      break;
    }
  }

  // Reduce to a single global answer.
  bool isColoringValid = MpiWrapper::allReduce( &isLocalColoringValid, MpiWrapper::Reduction::LogicalAnd, comm );

  return isColoringValid;
}


size_t GraphColoringBase::getNumberOfColors( const int color, MPI_Comm comm )
{
  std::vector< int > colors = {color};
  return getNumberOfColors( colors, comm );
}


size_t GraphColoringBase::getNumberOfColors( const std::vector< int > & colors, MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const size = MpiWrapper::commSize( comm );

  std::set< int > localDistinctColors = std::set< int >( colors.begin(), colors.end());
  std::vector< int > localDistinctColorsVector( localDistinctColors.begin(), localDistinctColors.end());
  int const localSize = localDistinctColorsVector.size();

  // Gather the sizes of the local color vectors from all ranks
  std::vector< int > allSizes( size );
  MpiWrapper::gather( &localSize, 1, allSizes.data(), 1, 0, comm );

  // Calculate the total number of colors and the displacements for gathering
  int totalSize = 0;
  std::vector< int > displacements( size, 0 );
  if( rank == 0 )
  {
    for( int i = 0; i < size; ++i )
    {
      displacements[i] = totalSize;
      totalSize += allSizes[i];
    }
  }

  // Gather all colors from all ranks to rank 0
  std::vector< int > allColors( totalSize );
  MpiWrapper::gatherv( localDistinctColorsVector.data(), localSize,
                       allColors.data(), allSizes.data(), displacements.data(),
                       0, comm );

  // Determine the number of distinct colors on rank 0
  int numDistinctColors = 0;
  if( rank == 0 )
  {
    std::set< int > uniqueColors( allColors.begin(), allColors.end());
    numDistinctColors = uniqueColors.size();
  }

  // Broadcast the number of distinct colors to all ranks
  MpiWrapper::bcast( &numDistinctColors, 1, 0, comm );

  return numDistinctColors;
}


} // namespace graph
} // namespace geos
