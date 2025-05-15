
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
 * @file PTScotchPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_PTSCOTCHPARTITIONER_HPP_
#define GEOS_PARTITIONER_PTSCOTCHPARTITIONER_HPP_

#include "PartitionerBase.hpp"


namespace geos
{


class PTScotchPartitioner : public PartitionerBase
{
public:
  PTScotchPartitioner();
  ~PTScotchPartitioner();


  void setNeighborsRank( const std::vector< int > & neighborsRank ) override;


  void setPartitionCounts( unsigned int xPartitions, unsigned int yPartitions, unsigned int zPartitions ) override;



/**
 * @brief Partition a mesh according to its dual graph.
 * @param graph the input graph (edges of locally owned nodes)
 * @param numParts target number of partitions
 * @param comm the MPI communicator of processes to partition over
 * @return an array of target partitions for each element in local mesh
 */

//array1d< int64_t >
//partition( ArrayOfArraysView< int64_t const, int64_t > const & graph,
//           int64_t const numParts,
//           MPI_Comm comm )

  array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                   arrayView1d< pmet_idx_t const > const & GEOS_UNUSED_PARAM( vertDist ),
                                   pmet_idx_t const numParts,
                                   MPI_Comm comm,
                                   int const numRefinements ) override;



private:

  void color() override;

};


} // namespace geos


#endif // GEOS_PARTITIONER_PTSCOTCHPARTITIONER_HPP_
