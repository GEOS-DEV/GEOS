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
#include "mesh/DomainPartition.hpp"


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

  ~SpatialPartition() override;

  bool isCoordInPartition( const real64 & coord, const int dir ) override;

  bool isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                      const real64 & boundaryRadius );

  void updateSizes( arrayView1d< real64 > const domainL,
                    real64 const dt );

  void setSizes( real64 const ( &min )[ 3 ],
                 real64 const ( &max )[ 3 ] ) override;

  real64* getLocalMin()
  {
    return m_min;
  }

  real64* getLocalMax()
  {
    return m_max;
  }

  real64* getGlobalMin()
  {
    return m_gridMin;
  }

  real64* getGlobalMax()
  {
    return m_gridMax;
  }

  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions ) override;

  int getColor() override;

  void repartitionMasterParticles( ParticleSubRegion & subRegion,
                                   MPI_iCommData & icomm,
                                   string solidMaterialName );

  void getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                   MPI_iCommData & icomm,
                                                   const real64 & boundaryRadius );

  void sendCoordinateListToNeighbors( arrayView1d< R1Tensor > const & particleCoordinatesSendingToNeighbors,           // Single list of coordinates sent to all neighbors
                                      MPI_iCommData & icomm,                                                                      // Solver's MPI communicator
                                     std::vector<array1d<R1Tensor>>& particleCoordinatesReceivedFromNeighbors       // List of lists of coordinates received from each neighbor.
  );

  template <typename indexType>
  void sendListOfIndicesToNeighbors( std::vector< array1d< indexType > > & listSendingToEachNeighbor,
                                     MPI_iCommData & icomm,
                                     std::vector< array1d< indexType > > & listReceivedFromEachNeighbor );

  void sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                std::vector<int> const & newParticleStartingIndices,
                                std::vector<int> const & numberOfIncomingParticles,
                                MPI_iCommData & icomm,
                                std::vector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor );


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
};

}
#endif /* GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
