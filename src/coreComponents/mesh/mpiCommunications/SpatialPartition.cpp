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

#include "SpatialPartition.hpp"
#include "codingUtilities/Utilities.hpp"
#include "LvArray/src/genericTensorOps.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"

#include <cmath>

namespace geosx
{

  namespace
  {

    // Modulo
    // returns a positive value regardless of the sign of numerator
    real64 Mod( real64 num, real64 denom )
    {
      if( fabs( denom )<fabs( num )*1.0e-14 )
      {
        return num;
      }

      return num - denom * std::floor( num/denom );
    }

    // MapValueToRange
    // returns a periodic value in the range [min, max)
    real64 MapValueToRange( real64 value, real64 min, real64 max )
    {
      return Mod( value-min, max-min )+min;
    }

  }

  SpatialPartition::SpatialPartition():
    PartitionBase(),
    m_Partitions(),
    m_Periodic( nsdof ),
    m_coords( nsdof ),
    m_min{ 0.0 },
    m_max{ 0.0 },
    m_blockSize{ 1.0 },
    m_gridSize{ 0.0 },
    m_gridMin{ 0.0 },
    m_gridMax{ 0.0 }
  {
    m_size = 0;
    m_rank = 0;
    m_numColors = 8,
    setPartitions( 1, 1, 1 );
  }

  SpatialPartition::~SpatialPartition()
  {}

  void SpatialPartition::setPartitions( unsigned int xPartitions,
                                        unsigned int yPartitions,
                                        unsigned int zPartitions )
  {
    m_Partitions.resize( 3 );
    m_Partitions( 0 ) = xPartitions;
    m_Partitions( 1 ) = yPartitions;
    m_Partitions( 2 ) = zPartitions;
    m_size = 1;
    for( int i = 0; i < nsdof; i++ )
    {
      m_size *= m_Partitions( i );
    }
    setContactGhostRange( 0.0 );
  }

  int SpatialPartition::getColor()
  {
    int color = 0;

    if( isOdd( m_coords[0] ) )
    {
      color += 1;
    }

    if( isOdd( m_coords[1] ) )
    {
      color += 2;
    }

    if( isOdd( m_coords[2] ) )
    {
      color += 4;
    }

    m_numColors = 8;

    return color;
  }

  void SpatialPartition::addNeighbors( const unsigned int idim,
                                       MPI_Comm & cartcomm,
                                       int * ncoords )
  {

    if( idim == nsdof )
    {
      bool me = true;
      for( int i = 0; i < nsdof; i++ )
      {
        if( ncoords[i] != this->m_coords( i ))
        {
          me = false;
          break;
        }
      }
      if( !me )
      {
        int const rank = MpiWrapper::cartRank( cartcomm, ncoords );
        m_neighbors.push_back( NeighborCommunicator( rank ) );
      }
    }
    else
    {
      const int dim = this->m_Partitions( LvArray::integerConversion< localIndex >( idim ) );
      const bool periodic = this->m_Periodic( LvArray::integerConversion< localIndex >( idim ) );
      for( int i = -1; i < 2; i++ )
      {
        ncoords[idim] = this->m_coords( LvArray::integerConversion< localIndex >( idim ) ) + i;
        bool ok = true;
        if( periodic )
        {
          if( ncoords[idim] < 0 )
            ncoords[idim] = dim - 1;
          else if( ncoords[idim] >= dim )
            ncoords[idim] = 0;
        }
        else
        {
          ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
        }
        if( ok )
        {
          addNeighbors( idim + 1, cartcomm, ncoords );
        }
      }
    }
  }

  void SpatialPartition::setSizes( real64 const ( &min )[ 3 ],
                                   real64 const ( &max )[ 3 ] )
  {

    {
      //get size of problem and decomposition
      m_size = MpiWrapper::commSize( MPI_COMM_GEOSX );

      //check to make sure our dimensions agree
      {
        int check = 1;
        for( int i = 0; i < nsdof; i++ )
        {
          check *= this->m_Partitions( i );
        }
        GEOSX_ERROR_IF_NE( check, m_size );
      }

      //get communicator, rank, and coordinates
      MPI_Comm cartcomm;
      {
        int reorder = 0;
        MpiWrapper::cartCreate( MPI_COMM_GEOSX, nsdof, m_Partitions.data(), m_Periodic.data(), reorder, &cartcomm );
      }
      m_rank = MpiWrapper::commRank( cartcomm );
      MpiWrapper::cartCoords( cartcomm, m_rank, nsdof, m_coords.data());

      //add neighbors
      {
        int ncoords[nsdof];
        m_neighbors.clear();
        addNeighbors( 0, cartcomm, ncoords );
      }

      MpiWrapper::commFree( cartcomm );
    }

    // global values
    LvArray::tensorOps::copy< 3 >( m_gridMin, min );
    LvArray::tensorOps::copy< 3 >( m_gridMax, max );
    LvArray::tensorOps::copy< 3 >( m_gridSize, max );
    LvArray::tensorOps::subtract< 3 >( m_gridSize, min );

    // block values
    LvArray::tensorOps::copy< 3 >( m_blockSize, m_gridSize );

    LvArray::tensorOps::copy< 3 >( m_min, min );
    for( int i = 0; i < nsdof; ++i )
    {
      const int nloc = m_Partitions( i ) - 1;
      const localIndex nlocl = static_cast< localIndex >(nloc);
      if( m_PartitionLocations[i].empty() )
      {
        // the default "even" spacing
        m_blockSize[ i ] /= m_Partitions( i );
        m_min[ i ] += m_coords( i ) * m_blockSize[ i ];
        m_max[ i ] = min[ i ] + (m_coords( i ) + 1) * m_blockSize[ i ];

        m_PartitionLocations[i].resize( nlocl );
        for( localIndex j = 0; j < m_PartitionLocations[ i ].size(); ++j )
        {
          m_PartitionLocations[ i ][ j ] = (j+1) * m_blockSize[ i ];
        }
      }
      else if( nlocl == m_PartitionLocations[i].size() )
      {
        const int parIndex = m_coords[i];
        if( parIndex == 0 )
        {
          m_min[i] = min[i];
          m_max[i] = m_PartitionLocations[i][parIndex];
        }
        else if( parIndex == nloc )
        {
          m_min[i] = m_PartitionLocations[i][parIndex-1];
          m_max[i] = max[i];
        }
        else
        {
          m_min[i] = m_PartitionLocations[i][parIndex-1];
          m_max[i] = m_PartitionLocations[i][parIndex];
        }
      }
      else
      {
        GEOSX_ERROR( "SpatialPartition::setSizes(): number of partition locations does not equal number of partitions - 1\n" );
      }
    }
  }

  bool SpatialPartition::isCoordInPartition( const real64 & coord, const int dir )
  {
    bool rval = true;
    const int i = dir;
    if( m_Periodic( i ))
    {
      if( m_Partitions( i ) != 1 )
      {
        real64 localCenter = MapValueToRange( coord, m_gridMin[ i ], m_gridMax[ i ] );
        rval = rval && localCenter >= m_min[ i ] && localCenter < m_max[ i ];
      }

    }
    else
    {
      rval = rval && (m_Partitions[ i ] == 1 || (coord >= m_min[ i ] && coord < m_max[ i ]));
    }

    return rval;
  }

  void SpatialPartition::setContactGhostRange( const real64 bufferSize )
  {
    LvArray::tensorOps::copy< 3 >( m_contactGhostMin, m_min );
    LvArray::tensorOps::addScalar< 3 >( m_contactGhostMin, -bufferSize );

    LvArray::tensorOps::copy< 3 >( m_contactGhostMax, m_max );
    LvArray::tensorOps::addScalar< 3 >( m_contactGhostMax, bufferSize );
  }

  void SpatialPartition::repartitionMasterParticlesToNeighbors(DomainPartition & domain,
                                                               MPI_iCommData & icomm)
  {

    /*
     * Search for any particles owned by this partition, which are no longer in the
     * partition domain.  Send a copy of these particles to their new partition, but
     * keep them as ghosts on the current partition.
     *
     * This assumes that particles have already been partitioned consistent with the background
     * grid, but does not assume any specific partition topology, only that the neighbor list
     * is complete, and each partition can evaluate its own isCoordinateInPartition() function.
     *
     * After this function, each particle should be in its correct partition, and the
     * ghostRank of particles that were moved from the current partition will be correct,
     * but the master-ghost map in the neighbor list still needs to updated.  Particles that
     * were ghosts of an object that has been repartitioned will need to have their ghost
     * rank updated, and may need to be deleted/added elsewhere.
     *
     */

    // MPM-specific code where we assume there are 2 mesh bodies and only one of them has particles
    dataRepository::Group & meshBodies = domain.getMeshBodies();
    MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >(0);
    MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >(1);
    MeshBody & particles = meshBody1.m_hasParticles ? meshBody1 : meshBody2;
    ParticleManager & particleManager = particles.getMeshLevel(0).getParticleManager();

    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {

      // (1) Identify any particles that are master on the current partition, but whose center lies
      // outside of the partition domain.  ghostRank() for particles is defined such that it always
      // equals the rank of the owning process. Thus a particle is master iff ghostRank==partition.rank
      //
      // Temporarily set the ghost rank of any particle to be moved to "-1".  If the particle is
      // requested by another partition, its ghost rank will be updated.  Any particle that still
      // has a ghostRank=-1 at the end of this function is lost and needs to be deleted.  This
      // should only happen if it has left the global domain (hopefully at an outflow b.c.).

      arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
      arrayView1d< int > const ghostRank = subRegion.getGhostRank();
      array1d< R1Tensor > outOfDomainParticleCoordinates;
      std::vector< localIndex > outOfDomainParticleLocalIndices;
      unsigned int nn = m_neighbors.size();

      // Number of partition neighbors.
      for (localIndex pp = 0; pp < subRegion.size(); pp++)
      {
        bool inPartition = true;
        R1Tensor p_x;
        for(int i=0; i<3; i++)
        {
          p_x[i] = particleCenter[pp][i];
          inPartition = inPartition && isCoordInPartition( p_x[i], i);
        }
        if( ghostRank[pp]==this->m_rank && !inPartition )
        {
          outOfDomainParticleCoordinates.emplace_back(p_x);  // Store the coordinate of the out-of-domain particle
          outOfDomainParticleLocalIndices.push_back(pp);     // Store the local index "pp" for the current coordinate.
          ghostRank[pp] = -1;                                // Temporarily set ghostRank of out-of-domain particle to -1 until it is requested by someone.
        }
      }


      // (2) Pack the list of particle center coordinates to each neighbor, and send/receive the list to neighbors.
      std::vector<array1d<R1Tensor>> particleCoordinatesReceivedFromNeighbors(nn);

      sendCoordinateListToNeighbors(outOfDomainParticleCoordinates.toView(),      // input: Single list of coordinates sent to all neighbors
                                    icomm,                                        // input: Solver MPI communicator
                                    particleCoordinatesReceivedFromNeighbors      // output: List of lists of coordinates received from each neighbor.
                                    );

    } );

  }


  void SpatialPartition::sendCoordinateListToNeighbors(arrayView1d<R1Tensor> const & particleCoordinatesSendingToNeighbors,           // Single list of coordinates sent to all neighbors
                                                       MPI_iCommData & icomm,                                                         // Solver's MPI communicator
                                                       std::vector<array1d<R1Tensor>>& particleCoordinatesReceivedFromNeighbors   // List of lists of coordinates received from each neighbor.
  )
  {
    // Number of neighboring partitions
    unsigned int nn = m_neighbors.size();

    // The send buffer is identical for each neighbor, and contains the packed coordinate list.
    unsigned int sizeToBePacked = 0;                                                      // size of the outgoing data with packing=false (we need to run through it first without packing so we can size the buffer)
    unsigned int sizeOfPacked = 0;                                                        // size of the outgoing data with packing=true
    buffer_unit_type* junk;                                                               // junk buffer, stores nothing since we're just getting the buffer size on the first pass
    sizeToBePacked += bufferOps::Pack< false >( junk,
                                                particleCoordinatesSendingToNeighbors );  // get the size of the coordinate list
    buffer_type sendBuffer(sizeToBePacked);                                               // the actual sized buffer that we pack into
    buffer_unit_type* sendBufferPtr = sendBuffer.data();                                  // get a pointer to the buffer
    sizeOfPacked += bufferOps::Pack< true >( sendBufferPtr,
                                             particleCoordinatesSendingToNeighbors );     // pack the coordinate list into the buffer
    GEOSX_ERROR_IF_NE( sizeToBePacked, sizeOfPacked );                                    // make sure the packer is self-consistent

    // Declare the receive buffers
    array1d<unsigned int> sizeOfReceived(nn);
    array1d< buffer_type > receiveBuffer(nn);

    // send the coordinate list to each neighbor.  Using an asynchronous send,
    // the mpi request will be different for each send, but the buffer is the same
    {
      array1d<MPI_Request>   sendRequest(nn);
      array1d<MPI_Status>    sendStatus(nn);
      array1d<MPI_Request>   receiveRequest(nn);
      array1d<MPI_Status>    receiveStatus(nn);

      // Send/receive the size of the packed buffer
      for(size_t n=0; n<nn; n++ )
      {
        m_neighbors[n].mpiISendReceive( &sizeOfPacked,
                                        1,
                                        sendRequest[n],
                                        &(sizeOfReceived[n]),
                                        1,
                                        receiveRequest[n],
                                        icomm.commID(),
                                        MPI_COMM_GEOSX );
      }
      MPI_Waitall(nn, sendRequest.data(), sendStatus.data());
      MPI_Waitall(nn, receiveRequest.data(), receiveStatus.data());
    }

    // Send/receive the buffer containing the list of coordinates
    {
      array1d<MPI_Request>   sendRequest(nn);
      array1d<MPI_Status>    sendStatus(nn);
      array1d<MPI_Request>   receiveRequest(nn);
      array1d<MPI_Status>    receiveStatus(nn);

      for(size_t n=0; n<nn; n++ )
      {
        receiveBuffer[n].resize(sizeOfReceived[n]);
        m_neighbors[n].mpiISendReceive( sendBuffer.data(), // TODO: This can't be sendBufferPtr, why not? I guess cuz sendBufferPtr gets incremented (moved) during packing.
                                        sizeOfPacked,
                                        sendRequest[n],
                                        receiveBuffer[n].data(),
                                        sizeOfReceived[n],
                                        receiveRequest[n],
                                        icomm.commID(),
                                        MPI_COMM_GEOSX );
      }
      MPI_Waitall(nn, sendRequest.data(), sendStatus.data());
      MPI_Waitall(nn, receiveRequest.data(), receiveStatus.data());
    }

    // Unpack the received coordinate list from each neighbor
    for(size_t n=0; n<nn; n++ )
    {
      // Unpack the buffer to an array of coordinates.
      const buffer_unit_type* receiveBufferPtr = receiveBuffer[n].data(); // needed for const cast
      bufferOps::Unpack( receiveBufferPtr, particleCoordinatesReceivedFromNeighbors[n] );
    }
  }

}
