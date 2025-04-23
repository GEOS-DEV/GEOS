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
 * @file RLFGraphColoringMPI.cpp
 */

#include "RLFGraphColoringMPI.hpp"
#include "RLFGraphColoring.hpp"
#include "GraphToolsMPI.hpp"
#include <algorithm>


namespace geos
{
namespace graph
{

RLFGraphColoringMPI::RLFGraphColoringMPI( MPI_Comm comm ): GraphColoringBase( comm )
{}

RLFGraphColoringMPI::~RLFGraphColoringMPI()
{}


int RLFGraphColoringMPI::colorGraph( const std::vector< camp::idx_t > & localAdjncy )
{
  std::vector< camp::idx_t > localXadj = createXadjFromAdjncy( localAdjncy, m_comm );
  std::vector< int > localColors = RLFGraphColoringMPI::colorGraph( localXadj, localAdjncy );
  return localColors[0];
}


std::vector< int > RLFGraphColoringMPI::colorGraph( const std::vector< camp::idx_t > & localXadj,
                                                    const std::vector< camp::idx_t > & localAdjncy )
{
  int const rank = MpiWrapper::commRank( m_comm );
  int const size = MpiWrapper::commSize( m_comm );

  // Perform coloring on rank 0
  auto [xadj, adjncy] = gatherGraphData( localXadj, localAdjncy, m_comm );
  std::vector< int > colors;
  if( rank == 0 )
  {
    geos::graph::RLFGraphColoring graphColoring;
    colors = graphColoring.colorGraph( xadj, adjncy );
  }

  // Scatter colors back to original ranks
  std::vector< int > sendCounts;
  std::vector< int > displacements;

  if( rank==0 )
  {
    sendCounts.resize( size );
  }
  int localNodeCounts = static_cast< int >(localXadj.size()-1);
  MpiWrapper::gather( &localNodeCounts, 1, sendCounts.data(), 1, 0, m_comm );

  // Calculate displacements for the scattered data
  if( rank == 0 )
  {
    displacements.resize( size );
    displacements[0] = 0;
    for( int i = 1; i < size; ++i )
    {
      displacements[i] = displacements[i - 1] + sendCounts[i - 1];
    }
  }

  std::vector< int > localColors( localNodeCounts );
  MPI_Scatterv( colors.data(), sendCounts.data(), displacements.data(), MPI_INT,
                localColors.data(), localNodeCounts, MPI_INT, 0, m_comm );

  return localColors;
}


size_t RLFGraphColoringMPI::getNumberOfColors( const int color ) const
{
  return GraphColoringBase::getNumberOfColors( color, m_comm );
}

size_t RLFGraphColoringMPI::getNumberOfColors( const std::vector< int > & colors ) const
{
  return GraphColoringBase::getNumberOfColors( colors, m_comm );
}


bool RLFGraphColoringMPI::isColoringValid( const std::vector< camp::idx_t > & adjncy, const int color ) const
{
  return GraphColoringBase::isColoringValid( adjncy, color, m_comm );
}


} // namespace graph
} // namespace geos
