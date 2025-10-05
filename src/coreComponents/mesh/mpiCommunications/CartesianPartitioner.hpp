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


#include "GeometricPartitioner.hpp"
#include "common/MpiWrapper.hpp"


namespace geos
{

class CartesianPartitioner : public GeometricPartitioner
{
public:
  static constexpr int m_ndim = 3;

  CartesianPartitioner( string const & name,
                        Group * const parent );

  /**
   * @brief Structure to hold scoped key names
   */
  struct viewKeyStruct
  {
    constexpr static char const * partitionCountsString() { return "partitionCounts"; }
  };

  /**
   * @brief Return the name of the CartesianPartitioner in object Catalog.
   * @return string that contains the key name to CartesianPartitioner in the Catalog
   */
  static string catalogName() { return "Cartesian"; }

  void initializeDomain( const R1Tensor & globalMin, const R1Tensor & globalMax ) override;

  bool isCoordInLocalPartition( const R1Tensor & coords ) const override;

  bool isCoordInPartition( const real64 & coord, const int dir ) const;

  void processCommandLineOverrides( unsigned int xPartitionsCL,
                                    unsigned int yPartitionsCL,
                                    unsigned int zPartitionsCL ) override;

  void setPartitionCounts( unsigned int xPartitions,
                           unsigned int yPartitions,
                           unsigned int zPartitions ) override;

  void setNeighborsRank( const std::vector< int > & neighborsRank ) override;



  array1d< int > const & getPartitionCounts() const
  {
    return m_partitionCounts;
  }


  void setPeriodicity( int index, int value )
  {
    m_periodic[index] = value;
  }

  array1d< int > const & getPeriodicity() const
  {
    return m_periodic;
  }


  real64 const * getLocalMin() const
  {
    return m_localMin;
  }

  real64 const * getLocalMax() const
  {
    return m_localMax;
  }

  real64 const * getGlobalMin() const
  {
    return m_globalGridMin;
  }

  real64 const * getGlobalMax() const
  {
    return m_globalGridMax;
  }


private:
  void color() override;

  void initializeCartesianCommunicator( MPI_Comm & cartComm );
  void determineNeighborsRank( MPI_Comm cartComm );
  void addNeighbors( const unsigned int idim, int * ncoords, MPI_Comm cartComm );

  void validatePartitionSize( int size ) const;
  void setGlobalGridValues( const R1Tensor & globalGridMin, const R1Tensor & globalGridMax );
  void setLocalPartitionValues( const R1Tensor & globalGridMin );

protected:
  /// Periodicity in each dimension (0 = non-periodic, 1 = periodic)
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
