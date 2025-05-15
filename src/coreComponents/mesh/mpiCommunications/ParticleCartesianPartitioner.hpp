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
 * @file ParticleCartesianPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_PARTICLECARTESIANPARTITIONER_HPP_
#define GEOS_PARTITIONER_PARTICLECARTESIANPARTITIONER_HPP_

#include "CartesianPartitioner.hpp"

#include "mesh/DomainPartition.hpp"


namespace geos
{

class ParticleCartesianPartitioner : public CartesianPartitioner
{
public:
  ParticleCartesianPartitioner();
  ~ParticleCartesianPartitioner();

  bool isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                      const real64 & boundaryRadius ) const;
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
                                      stdVector< array1d< R1Tensor > > & particleCoordinatesReceivedFromNeighbors
                                      );

  template< typename indexType >
  void sendListOfIndicesToNeighbors( stdVector< array1d< indexType > > & listSendingToEachNeighbor,
                                     MPI_iCommData & commData,
                                     stdVector< array1d< indexType > > & listReceivedFromEachNeighbor );

  void sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                stdVector< int > const & newParticleStartingIndices,
                                stdVector< int > const & numberOfIncomingParticles,
                                MPI_iCommData & commData,
                                stdVector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor );



};

} // namespace geos

#endif // GEOS_PARTITIONER_PARTICLECARTESIANPARTITIONER_HPP_
