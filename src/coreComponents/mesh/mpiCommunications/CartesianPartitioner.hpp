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
 * @file CartesianPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_CARTESIANPARTITIONER_HPP_
#define GEOS_PARTITIONER_CARTESIANPARTITIONER_HPP_


#include "PartitionerBase.hpp"
#include "common/MpiWrapper.hpp"


namespace geos
{

class CartesianPartitioner : public PartitionerBase
{
public:

  CartesianPartitioner();
  ~CartesianPartitioner();

  static constexpr int m_ndim = 3;

  void partition( const real64 (& globalGridMin)[m_ndim], const real64 (& globalGridMax)[m_ndim] );

  array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                   arrayView1d< pmet_idx_t const > const & vertDist,
                                   pmet_idx_t const numParts,
                                   MPI_Comm comm,
                                   int const numRefinements ) override;

  void setPartitionCounts( unsigned int xPartitions, unsigned int yPartitions, unsigned int zPartitions ) override;

  void setNeighborsRank( const std::vector< int > & neighborsRank ) override;


  bool isCoordInPartition( const real64 & coord, const int dir ) const;


  array1d< int > const & getPartitionCounts() const
  {
    return m_partitionCounts;
  }

  array1d< int > & getPartitionCounts()
  {
    return m_partitionCounts;
  }


  void setPeriodicity( int index, int value )
  {
    m_periodic[index] = value;
  }

  array1d< int > & getPeriodicity()
  {
    return m_periodic;
  }

  array1d< int > const & getPeriodicity() const
  {
    return m_periodic;
  }


private:
  void color() override;

  void initializeCartesianCommunicator();
  void determineNeighborsRank();
  void validatePartitionSize( int size ) const;
  void setGlobalGridValues( const real64 (& globalGridMin)[m_ndim], const real64 (& globalGridMax)[m_ndim] );
  void setLocalPartitionValues( const real64 (& globalGridMin)[m_ndim] );

  void addNeighbors( const unsigned int idim, int * ncoords );


protected:
  MPI_Comm m_cartComm;
  array1d< int > m_periodic;
  /// ijk partition indexes
  array1d< int > m_coords;

  array1d< int > m_partitionCounts;


  /// Minimum extent of partition dimensions (excluding ghost objects)
  real64 m_localMin[m_ndim];
  /// Maximum extent of partition dimensions (excluding ghost objects)
  real64 m_localMax[m_ndim];
  /// Length of partition dimensions (excluding ghost objects).
  real64 m_localSize[m_ndim];


  /// Minimum extent of problem dimensions (excluding ghost objects).
  real64 m_globalGridMin[m_ndim];
  /// Maximum extent of problem dimensions (excluding ghost objects).
  real64 m_globalGridMax[m_ndim];
  /// Total length of problem dimensions (excluding ghost objects).
  real64 m_globalGridSize[m_ndim];

  int m_cartRank;
};

} // namespace geos

#endif // GEOS_PARTITIONER_CARTESIANPARTITIONER_HPP_
