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
 * @file GraphPartitioner.cpp
 */

#include "GraphPartitioner.hpp"
#include "common/MpiWrapper.hpp"

#ifdef GEOS_USE_TRILINOS
#include "mesh/graphs/ZoltanGraphColoring.hpp"
#else
#include "mesh/graphs/RLFGraphColoringMPI.hpp"
#endif

namespace geos
{

using namespace dataRepository;

using camp::idx_t;

static_assert( std::is_same< idx_t, pmet_idx_t >::value, "Non-matching index types. ParMETIS must be built with 64-bit indices." );


void GraphPartitioner::processCommandLineOverrides( unsigned int xparCL, unsigned int yparCL, unsigned int zparCL )
{
  // If user provided command-line overrides for partition counts...
  if( xparCL != 0 || yparCL != 0 || zparCL != 0 )
  {
    // ...warn that they are ignored for non-geometric partitioners
    GEOS_WARNING( "Partition counts from the command line (-x, -y, -z) are ignored when using a graph-based partitioner." );
  }

  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );
  setPartitionCounts( 1, 1, mpiSize );
}

void GraphPartitioner::setPartitionCounts( unsigned int xPartitions, unsigned int yPartitions, unsigned int zPartitions )
{
  // Only the total number of partitions matters
  m_numPartitions = xPartitions * yPartitions * zPartitions;
}

void GraphPartitioner::setNeighborsRank( const std::vector< int > & neighborsRank )
{
  m_neighborsRank = neighborsRank;
  // Trigger coloring whenever the neighbor list is updated
  color();
}

void GraphPartitioner::color()
{
  // Build adjacency list from neighbor ranks
  std::vector< camp::idx_t > adjncy;
  adjncy.reserve( m_neighborsRank.size() );
  std::copy( m_neighborsRank.begin(), m_neighborsRank.end(), std::back_inserter( adjncy ) );

  // Use the appropriate graph coloring algorithm based on third-party library availability
#ifdef GEOS_USE_TRILINOS
  geos::graph::ZoltanGraphColoring coloring;
#else
  geos::graph::RLFGraphColoringMPI coloring;
#endif

  m_color = coloring.colorGraph( adjncy );

  // Verify that the coloring is valid (no two neighbors have the same color)
  if( !coloring.isColoringValid( adjncy, m_color ) )
  {
    GEOS_ERROR( "Graph coloring failed: an adjacent partition has the same color." );
  }

  m_numColors = coloring.getNumberOfColors( m_color );
}

} // namespace geos
