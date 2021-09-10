/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_


#include "PartitionBase.hpp"

#include <map>

constexpr int nsdof = 3;
namespace geosx
{

/**
 * @brief Concrete (cartesian?) partitioning.
 */
class SpatialPartition : public PartitionBase
{
public:
  SpatialPartition();
  virtual ~SpatialPartition();

  virtual bool isCoordInPartition( const real64 & coord, const int dir );
  virtual bool isCoordInPartition( real64 const ( &coordinates )[ 3 ] );
  virtual bool isCoordInPartition( real64 const ( &coordinates )[ 3 ],
                                   const int numDistPartition );

  /**
   * @brief Variant of IsCoordInPartition with intervals closed at both ends
   * @param coordinates The point coordinates.
   * @return The predicate result.
   *
   * @note Unused.
   */
  virtual bool isCoordInPartitionClosed( real64 const ( &coordinates )[ 3 ] );

  /**
   * @brief Checks if the point located inside the current partition bounding box.
   * @param coordinates The element center coordinates.
   * @return The predicate result.
   *
   * @note Unused.
   */
  virtual bool isCoordInPartitionBoundingBox( real64 const ( &coordinates )[ 3 ] );

  virtual bool isCoordInContactGhostRange( real64 const ( &coordinates )[ 3 ] );

  void setSizes( real64 const ( &min )[ 3 ],
                 real64 const ( &max )[ 3 ] );

  /**
   * @brief Defines the boundaries of the partition
   * @param min The minimum.
   * @param max The maximum.
   */
  void setPartitionGeometricalBoundary( real64 const ( &min )[ 3 ],
                                        real64 const ( &max )[ 3 ] );

  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions )
  {
    m_Partitions.resize( 3 );
    m_Partitions( 0 ) = xPartitions;
    m_Partitions( 1 ) = yPartitions;
    m_Partitions( 2 ) = zPartitions;
    m_size = 1;
    for( int i = 0; i < nsdof; i++ )
      m_size *= m_Partitions( i );
    setContactGhostRange( 0.0 );
  }

  /**
   * @brief Defines periodicity along the three (x, y, z) axis. An argument of 1 indicates periodicity.
   * @param xPeriodic Periodicity in x.
   * @param yPeriodic Periodicity in y.
   * @param zPeriodic Periodicity in z.
   */
  void setPeriodic( unsigned int xPeriodic,
                    unsigned int yPeriodic,
                    unsigned int zPeriodic )
  {
    m_Periodic( 0 ) = xPeriodic;
    m_Periodic( 1 ) = yPeriodic;
    m_Periodic( 2 ) = zPeriodic;
  }

  virtual void setContactGhostRange( const real64 bufferSize );

  int getColor();

protected:
  void initializePostSubGroups();

public:
  /// number of partitions
  array1d< int > m_Partitions;
  /**
   * @brief Boolean like array of length 3 (space dimensions).
   *
   * 1 means periodic.
   */
  array1d< int > m_Periodic;
  /// ijk partition indexes
  array1d< int > m_coords;

  /// Minimum extent of partition dimensions (excluding ghost objects)
  real64 m_min[3];
  /// Maximum extent of partition dimensions (excluding ghost objects)
  real64 m_max[3];

  /**
   * @brief Bounding box (minimum value) along the x direction.
   * @note Probably unused.
   */
  real64 m_xBoundingBoxMin[3];
  /**
   * @brief Bounding box (maximum value) along the x direction.
   * @note Probably unused.
   */
  real64 m_xBoundingBoxMax[3];

  /// Locations of partition boundaries
  array1d< real64 > m_PartitionLocations[3];

  /// Length of partition dimensions (excluding ghost objects).
  real64 m_blockSize[3];

  /// Total length of problem dimensions (excluding ghost objects).
  real64 m_gridSize[3];
  /// Minimum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMin[3];
  /// Maximum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMax[3];

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
  void addNeighbors( const unsigned int idim, MPI_Comm & cartcomm,
                     int * ncoords );

};
}
#endif /* GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
