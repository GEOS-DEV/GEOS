
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
 * @file PartitionerBase.hpp
 */

#ifndef GEOS_PARTITIONER_PARTITIONERBASE_HPP_
#define GEOS_PARTITIONER_PARTITIONERBASE_HPP_

//#include "common/StdContainerWrappers.hpp"   // for stdVector
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"



namespace geos
{

#if defined(GEOS_USE_HIP) // still need int32 hypre for the current hip-capable build
/// Typedef to allow us to specify required parmetis integer type.
using pmet_idx_t = int32_t;
#else
/// Typedef to allow us to specify required parmetis integer type.
using pmet_idx_t = int64_t;
#endif

class PartitionerBase
{
public:
  /**
   * @brief Partition a mesh according to its dual graph.
   * @param graph the input graph (edges of locally owned nodes)
   * @param vertDist the parallel distribution of vertices: vertex index offset on each rank
   * @param numParts target number of partitions
   * @param comm the MPI communicator of processes to partition over
   * @param numRefinements number of partition refinement iterations
   * @return an array of target partitions for each element in local mesh
   */
  virtual array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                           arrayView1d< pmet_idx_t const > const & vertDist,
                                           pmet_idx_t const numParts,
                                           MPI_Comm comm,
                                           int const numRefinements ) = 0;

  virtual void setPartitionCounts( unsigned int xPartitions,
                                   unsigned int yPartitions,
                                   unsigned int zPartitions ) = 0;


  virtual void setNeighborsRank( const std::vector< int > & neighborsRank ) = 0;

  std::vector< int > const & getNeighborsRank() const
  {
    return m_neighborsRank;
  }

  std::vector< int > & getNeighborsRank()
  {
    return m_neighborsRank;
  }


  void buildNeighbors()
  {
    m_neighbors.clear();

    for( int rank : m_neighborsRank )
    {
      m_neighbors.push_back( NeighborCommunicator( rank ));
    }
  }


  stdVector< NeighborCommunicator > & getNeighbors()
  { return m_neighbors; }

  stdVector< NeighborCommunicator > const & getNeighbors() const
  { return m_neighbors; };


  virtual void color() = 0;

  int getColor() const
  {
    return m_color;
  }
  int getColor()
  {
    return m_color;
  }

  int getNumColors() const
  {
    return m_numColors;
  }
  int getNumColors()
  {
    return m_numColors;
  }


  int getNumPartitions() const
  {
    return m_numPartitions;
  }
  int getNumPartitions()
  {
    return m_numPartitions;
  }



  virtual ~PartitionerBase() = default;

  PartitionerBase(): m_numColors( 1 ), m_color( 0 ), m_numPartitions( 1 ) {}



protected:
  stdVector< NeighborCommunicator > m_neighbors;
  std::vector< int > m_neighborsRank;


  int m_numColors;
  int m_color;
  int m_numPartitions;
};

}
#endif // GEOS_PARTITIONER_PARTITIONERBASE_HPP_
