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

#ifndef GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_


#include "PartitionBase.hpp"
#include "mesh/DomainPartition.hpp"

#include <array>
#include <map>

constexpr int nsdof = 3;
namespace geos
{

/**
 * @brief Concrete (cartesian?) partitioning.
 */
class SpatialPartition : public PartitionBase
{
public:
  SpatialPartition();

  ~SpatialPartition() override;

  bool isCoordInPartition( const real64 & coord, const int dir ) const override;

  bool isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                      const real64 & boundaryRadius ) const;

  void updateSizes( arrayView1d< real64 > const domainL,
                    real64 const dt );

  void setSizes( real64 const ( &min )[ 3 ],
                 real64 const ( &max )[ 3 ] ) override;

  real64 * getLocalMin()
  {
    return m_min;
  }

  real64 * getLocalMax()
  {
    return m_max;
  }

  real64 * getGlobalMin()
  {
    return m_gridMin;
  }

  real64 * getGlobalMax()
  {
    return m_gridMax;
  }

  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions ) override;

  int getColor() override;

  void repartitionMasterParticles( ParticleSubRegion & subRegion,
                                   MPI_iCommData & commData );

  void getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                   MPI_iCommData & commData,
                                                   const real64 & boundaryRadius );

  /**
   * @brief Send coordinates to neighbors as part of repartition.
   * @param[in] particleCoordinatesSendingToNeighbors Single list of coordinates sent to all neighbors
   * @param[in] commData Solver's MPI communicator
   * @param[in] particleCoordinatesReceivedFromNeighbors List of lists of coordinates received from each neighbor
   */
  void sendCoordinateListToNeighbors( arrayView1d< R1Tensor > const & particleCoordinatesSendingToNeighbors,
                                      MPI_iCommData & commData,
                                      std::vector< array1d< R1Tensor > > & particleCoordinatesReceivedFromNeighbors
                                      );

  template< typename indexType >
  void sendListOfIndicesToNeighbors( std::vector< array1d< indexType > > & listSendingToEachNeighbor,
                                     MPI_iCommData & commData,
                                     std::vector< array1d< indexType > > & listReceivedFromEachNeighbor );

  void sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                std::vector< int > const & newParticleStartingIndices,
                                std::vector< int > const & numberOfIncomingParticles,
                                MPI_iCommData & commData,
                                std::vector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor );

  /**
   * @brief Get the metis neighbors indices, const version. @see DomainPartition#m_metisNeighborList
   * @return Container of global indices.
   */
  std::set< int > const & getMetisNeighborList() const
  {
    return m_metisNeighborList;
  }

  /**
   * @brief Sets the list of metis neighbor list.
   * @param metisNeighborList A reference to the Metis neighbor list.
   */
  void setMetisNeighborList( std::set< int > const & metisNeighborList )
  {
    m_metisNeighborList = metisNeighborList;
  }

  /**
   * @brief Get the number of domains in each dimension for a regular partition with InternalMesh.
   * @return An array containing number of partition in X, Y and Z directions.
   */
  array1d< int > const & getPartitions() const
  {
    return m_Partitions;
  }

  void setGrid( std::array< real64, 9 > const & grid )
  {
    m_gridSize[0] = grid[0];
    m_gridSize[1] = grid[1];
    m_gridSize[2] = grid[2];

    m_gridMin[0] = grid[3];
    m_gridMin[1] = grid[4];
    m_gridMin[2] = grid[5];

    m_gridMax[0] = grid[6];
    m_gridMax[1] = grid[7];
    m_gridMax[2] = grid[8];
  }

  void setBlockSize( std::array< real64, 3 > const & blockSize )
  {
    m_blockSize[0] = blockSize[0];
    m_blockSize[1] = blockSize[1];
    m_blockSize[2] = blockSize[2];
  }

  void setBoundingBox( std::array< real64, 6 > const & bb )
  {
//    real64 const min[3]{ bb[0], bb[1], bb[2] };
//    real64 const max[3]{ bb[3], bb[4], bb[5] };

//    setSizes( min, max );

    m_min[0] = bb[0];
    m_min[1] = bb[1];
    m_min[2] = bb[2];

    m_max[0] = bb[3];
    m_max[1] = bb[4];
    m_max[2] = bb[5];
  }

  /**
   * @brief Boolean like array of length 3 (space dimensions).
   *
   * 1 means periodic.
   */
  array1d< int > m_Periodic;
  /// ijk partition indexes
  array1d< int > m_coords;

private:

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );

  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  void setContactGhostRange( const real64 bufferSize );

  /// Minimum extent of partition dimensions (excluding ghost objects)
  real64 m_min[3];
  /// Maximum extent of partition dimensions (excluding ghost objects)
  real64 m_max[3];

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
   * @brief Ghost position (min).
   */
  real64 m_contactGhostMin[3];

  /**
   * @brief Ghost position (max).
   */
  real64 m_contactGhostMax[3];

  /// number of partitions
  array1d< int > m_Partitions;

  /**
   * @brief Contains the global indices of the metis neighbors in case `metis` is used. Empty otherwise.
   */
  std::set< int > m_metisNeighborList;

};

}
#endif /* GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
