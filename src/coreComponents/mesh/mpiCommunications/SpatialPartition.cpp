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

#include "SpatialPartition.hpp"
// #include "mesh/DomainPartition.hpp"
#include "codingUtilities/Utilities.hpp"
#include "LvArray/src/genericTensorOps.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "mesh/generators/CellBlockManager.hpp"

#ifdef GEOS_USE_TRILINOS
#include "mesh/graphs/ZoltanGraphColoring.hpp"
#else
#include "mesh/graphs/RLFGraphColoringMPI.hpp"
#endif

#include <cmath>
#include <string>
#include <array>
#include <cstring>
#include <functional>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

namespace geos
{

using namespace dataRepository;

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

/**
 * @brief Compact metadata record used during MPM ghost particle discovery.
 *
 * The receiver needs:
 *  - x: position, to decide whether this particle lies in its ghost box.
 *  - globalID: to decide whether this particle already exists locally.
 *  - ownerLocalIndex: to request the full particle payload from the owner without
 *    forcing the owner to do a globalID -> localID lookup later.
 *
 * This struct is intentionally file-local so SpatialPartition.hpp does not need
 * to change.
 */
struct GhostParticleCandidate
{
  real64 x[3];
  globalIndex globalID;
  localIndex ownerLocalIndex;
};

static_assert( std::is_trivially_copyable< GhostParticleCandidate >::value,
               "GhostParticleCandidate must remain trivially copyable because it is exchanged as raw bytes." );


/**
 * @brief Convert a byte count to a buffer_unit_type count.
 *
 * GEOS packed buffers are counted in buffer_unit_type units, not necessarily
 * bytes. This helper keeps the raw-POD exchange below independent of the exact
 * buffer_unit_type size.
 */
inline unsigned int bufferUnitsForBytes( size_t const numberOfBytes )
{
  size_t const unitSize = sizeof( buffer_unit_type );
  size_t const numberOfUnits = numberOfBytes == 0 ? 0 : 1 + ( numberOfBytes - 1 ) / unitSize;

  GEOS_ERROR_IF( numberOfUnits > std::numeric_limits< unsigned int >::max(),
                 "SpatialPartition: raw POD message is too large for unsigned int buffer counts." );

  return static_cast< unsigned int >( numberOfUnits );
}


/**
 * @brief Pack an array1d of trivially-copyable records into a GEOS buffer.
 *
 * This is used only for file-local POD metadata. Full particle data is still
 * packed through ParticleSubRegionBase::particlePack().
 */
template< typename POD_TYPE >
unsigned int packPODArrayToBuffer( array1d< POD_TYPE > const & values,
                                   buffer_type & buffer )
{
  static_assert( std::is_trivially_copyable< POD_TYPE >::value,
                 "packPODArrayToBuffer only supports trivially-copyable types." );

  size_t const numberOfValues = static_cast< size_t >( values.size() );

  GEOS_ERROR_IF( numberOfValues > std::numeric_limits< size_t >::max() / sizeof( POD_TYPE ),
                 "SpatialPartition: POD buffer byte-count overflow." );

  size_t const numberOfBytes = numberOfValues * sizeof( POD_TYPE );
  unsigned int const numberOfUnits = bufferUnitsForBytes( numberOfBytes );

  buffer.resize( numberOfUnits );

  // Zero the whole buffer so any padding bytes in the last buffer_unit_type are deterministic.
  if( numberOfUnits > 0 )
  {
    std::memset( buffer.data(), 0, static_cast< size_t >( numberOfUnits ) * sizeof( buffer_unit_type ) );
  }

  if( numberOfBytes > 0 )
  {
    std::memcpy( buffer.data(), values.data(), numberOfBytes );
  }

  return numberOfUnits;
}


/**
 * @brief Unpack a GEOS buffer into an array1d of trivially-copyable records.
 */
template< typename POD_TYPE >
void unpackPODArrayFromBuffer( buffer_type const & buffer,
                               unsigned int const numberOfValues,
                               array1d< POD_TYPE > & values )
{
  static_assert( std::is_trivially_copyable< POD_TYPE >::value,
                 "unpackPODArrayFromBuffer only supports trivially-copyable types." );

  size_t const n = static_cast< size_t >( numberOfValues );

  GEOS_ERROR_IF( n > std::numeric_limits< size_t >::max() / sizeof( POD_TYPE ),
                 "SpatialPartition: POD buffer byte-count overflow." );

  size_t const numberOfBytes = n * sizeof( POD_TYPE );
  unsigned int const expectedUnits = bufferUnitsForBytes( numberOfBytes );

  GEOS_ERROR_IF( static_cast< size_t >( buffer.size() ) < static_cast< size_t >( expectedUnits ),
                 "SpatialPartition: received POD buffer is smaller than expected." );

  values.resize( numberOfValues );

  if( numberOfBytes > 0 )
  {
    std::memcpy( values.data(), buffer.data(), numberOfBytes );
  }
}


/**
 * @brief Exchange one POD list per neighbor.
 *
 * This mirrors the communication pattern used by sendListOfIndicesToNeighbors(),
 * but it does not rely on bufferOps serialization for the file-local
 * GhostParticleCandidate type.
 *
 * The exchange is two-stage:
 *  1. exchange number of POD records,
 *  2. exchange raw POD payload buffers.
 *
 * This is intended for homogeneous MPI jobs. Full particle payloads are not sent
 * through this helper.
 */
template< typename POD_TYPE >
void exchangePODListsWithNeighbors( stdVector< array1d< POD_TYPE > > const & listSendingToEachNeighbor,
                                    MPI_iCommData & commData,
                                    stdVector< array1d< POD_TYPE > > & listReceivedFromEachNeighbor,
                                    stdVector< NeighborCommunicator > & neighbors )
{
  static_assert( std::is_trivially_copyable< POD_TYPE >::value,
                 "exchangePODListsWithNeighbors only supports trivially-copyable types." );

  GEOS_ERROR_IF( neighbors.size() > std::numeric_limits< unsigned int >::max(),
                 "SpatialPartition: too many neighbors for unsigned int message counts." );

  unsigned int const nn = static_cast< unsigned int >( neighbors.size() );

  GEOS_ERROR_IF( listSendingToEachNeighbor.size() != nn,
                 "SpatialPartition: send-list count does not match neighbor count." );

  listReceivedFromEachNeighbor.resize( nn );

  if( nn == 0 )
  {
    return;
  }

  GEOS_ERROR_IF( nn > static_cast< unsigned int >( std::numeric_limits< int >::max() ),
                 "SpatialPartition: too many neighbors for MPI_Waitall count." );

  int const mpiNeighborCount = static_cast< int >( nn );

  stdVector< unsigned int > numberSending( nn );
  stdVector< unsigned int > numberReceived( nn );

  stdVector< unsigned int > sizeOfPackedSending( nn );
  stdVector< unsigned int > sizeOfPackedReceived( nn );

  stdVector< buffer_type > sendBuffer( nn );
  stdVector< buffer_type > receiveBuffer( nn );

  for( size_t n = 0; n < nn; ++n )
  {
    size_t const currentSize = static_cast< size_t >( listSendingToEachNeighbor[n].size() );

    GEOS_ERROR_IF( currentSize > std::numeric_limits< unsigned int >::max(),
                   "SpatialPartition: too many POD records for unsigned int message counts." );

    numberSending[n] = static_cast< unsigned int >( currentSize );
    sizeOfPackedSending[n] = packPODArrayToBuffer( listSendingToEachNeighbor[n], sendBuffer[n] );
  }

  // First exchange the number of POD records each neighbor will send.
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status > sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status > receiveStatus( nn );

    for( size_t n = 0; n < nn; ++n )
    {
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      neighbors[n].mpiISendReceive( &( numberSending[n] ),
                                    1,
                                    sendRequest[n],
                                    &( numberReceived[n] ),
                                    1,
                                    receiveRequest[n],
                                    commData.commID(),
                                    MPI_COMM_GEOS );
    }

    MPI_Waitall( mpiNeighborCount, sendRequest.data(), sendStatus.data() );
    MPI_Waitall( mpiNeighborCount, receiveRequest.data(), receiveStatus.data() );
  }

  // Allocate receive buffers using the received POD counts.
  for( size_t n = 0; n < nn; ++n )
  {
    size_t const numberOfBytes = static_cast< size_t >( numberReceived[n] ) * sizeof( POD_TYPE );
    sizeOfPackedReceived[n] = bufferUnitsForBytes( numberOfBytes );
    receiveBuffer[n].resize( sizeOfPackedReceived[n] );
  }

  // Exchange the raw POD payloads.
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status > sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status > receiveStatus( nn );

    for( size_t n = 0; n < nn; ++n )
    {
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      neighbors[n].mpiISendReceive( sendBuffer[n].data(),
                                    sizeOfPackedSending[n],
                                    sendRequest[n],
                                    receiveBuffer[n].data(),
                                    sizeOfPackedReceived[n],
                                    receiveRequest[n],
                                    commData.commID(),
                                    MPI_COMM_GEOS );
    }

    MPI_Waitall( mpiNeighborCount, sendRequest.data(), sendStatus.data() );
    MPI_Waitall( mpiNeighborCount, receiveRequest.data(), receiveStatus.data() );
  }

  for( size_t n = 0; n < nn; ++n )
  {
    unpackPODArrayFromBuffer( receiveBuffer[n],
                              numberReceived[n],
                              listReceivedFromEachNeighbor[n] );
  }
}

}

SpatialPartition::SpatialPartition( string const & name,
                                    Group * const parent ):
  PartitionBase( name, parent ),
  m_min( m_nsdof ),
  m_max( m_nsdof ),
  m_blockSize( m_nsdof ),
  m_gridMin( m_nsdof ),
  m_gridMax( m_nsdof ),
  m_gridSize( m_nsdof ),
  m_coords( m_nsdof ),
  m_partitions( m_nsdof ),
  m_periodic( m_nsdof ),
  m_contactGhostMin( m_nsdof ),
  m_contactGhostMax( m_nsdof )
{
  m_size = 0;
  m_rank = 0;
  m_numColors = 8;
  setPartitions( 1, 1, 1 );

  // Do m_coords, m_partitions need to be registered?

  registerWrapper( viewKeyStruct::minString(), &m_min ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Minimum extent of partition dimensions (excluding ghost objects)" );

  registerWrapper( viewKeyStruct::maxString(), &m_max ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Maximum extent of partition dimensions (excluding ghost objects)" );

  registerWrapper( viewKeyStruct::blockSizeString(), &m_blockSize ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Length of partition dimensions (excluding ghost objects)." );

  registerWrapper( viewKeyStruct::partitionLocationsString(), &m_partitionLocations ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Locations of partition boundaries" );

  registerWrapper( viewKeyStruct::gridMinString(), &m_gridMin ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Minimum extent of problem dimensions (excluding ghost objects)." );

  registerWrapper( viewKeyStruct::gridMaxString(), &m_gridMax ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Maximum extent of problem dimensions (excluding ghost objects)." );

  registerWrapper( viewKeyStruct::gridSizeString(), &m_gridSize ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Total length of problem dimensions (excluding ghost objects)." );

  registerWrapper( viewKeyStruct::periodicString(), &m_periodic ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "periodic flag for each direction of mesh" );

  registerWrapper( viewKeyStruct::contactGhostMinString(), &m_contactGhostMin ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Ghost position min." );

  registerWrapper( viewKeyStruct::contactGhostMaxString(), &m_contactGhostMax ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Ghost position max." );
}

SpatialPartition::~SpatialPartition()
{}

void SpatialPartition::postInputInitialization()
{
  PartitionBase::postInputInitialization();

  // Do LvArrays explicitly need to be resized to 0?
  if( m_partitionLocations.size() == 0 )
  {
    m_partitionLocations.resize( 3 );
    for( int i= 0; i < 3; i++ )
    {
      m_partitionLocations[i].resize( 0 );
    }
  }

  if( m_periodic.size() == 0 )
  {
    m_periodic.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_periodic, 0 );
  }
  else
  {
    GEOS_ERROR_IF( m_periodic.size() !=3, "Periodic flags must have size 3" );
  }

  if( m_min.size() == 0 )
  {
    m_min.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_min, 0.0 );
  }

  if( m_max.size() == 0 )
  {
    m_max.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_max, 0.0 );
  }

  if( m_blockSize.size() == 0 )
  {
    m_blockSize.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_blockSize, 1.0 );
  }
  if( m_gridSize.size() == 0 )
  {
    m_gridSize.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_gridSize, 0.0 );
  }

  if( m_gridMin.size() == 0 )
  {
    m_gridMin.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_gridMin, 0.0 );
  }

  if( m_gridMax.size() == 0 )
  {
    m_gridMax.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_gridMax, 0.0 );
  }

  if( m_contactGhostMin.size() == 0 )
  {
    m_contactGhostMin.resize( 3 );
  }

  if( m_contactGhostMax.size() == 0 )
  {
    m_contactGhostMax.resize( 3 );
  }
}

void SpatialPartition::setPartitions( unsigned int xPartitions,
                                      unsigned int yPartitions,
                                      unsigned int zPartitions )
{
  m_partitions.resize( 3 );
  m_partitions( 0 ) = xPartitions;
  m_partitions( 1 ) = yPartitions;
  m_partitions( 2 ) = zPartitions;
  m_size = 1;
  for( int i = 0; i < m_nsdof; i++ )
  {
    m_size *= m_partitions( i );
  }
  setContactGhostRange( 0.0 );
}

int SpatialPartition::getColor()
{
  if( m_metisNeighborList.empty() )
  {
    // Internal cartesian partitioner (for internal mesh)
    int color=0;
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

    // With this algorithm, numbering may have gaps.
    // In that case m_numColors is an upper bound, not the exact number of distinct colors used.
    m_numColors = MpiWrapper::max( color )+1;
    return color;
  }
  else
  {
    // External partitioner such as ParMetis or PTScotch (for VTK external mesh)
    std::vector< camp::idx_t > adjncy;
    adjncy.reserve( m_metisNeighborList.size());
    std::copy( m_metisNeighborList.begin(), m_metisNeighborList.end(), std::back_inserter( adjncy ));
#ifdef GEOS_USE_TRILINOS
    geos::graph::ZoltanGraphColoring coloring;
#else
    geos::graph::RLFGraphColoringMPI coloring;
#endif
    int color = coloring.colorGraph( adjncy );

    if( !coloring.isColoringValid( adjncy, color ))
    {
      GEOS_ERROR( "Invalid partition coloring: two neighboring partitions share the same color" );
    }
    m_numColors = coloring.getNumberOfColors( color );

    return color;
  }
}

void SpatialPartition::addNeighbors( const unsigned int idim,
                                     MPI_Comm & cartcomm,
                                     int * ncoords )
{
  if( idim == m_nsdof )
  {
    bool me = true;
    for( int i = 0; i < m_nsdof; i++ )
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
    const int dim = this->m_partitions( LvArray::integerConversion< localIndex >( idim ) );
    const bool periodic = this->m_periodic( LvArray::integerConversion< localIndex >( idim ) );
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
    m_size = MpiWrapper::commSize( MPI_COMM_GEOS );

    //check to make sure our dimensions agree
    {
      string_view partitionsLogMessage =
        "The total number of processes = {} does not correspond to the total number of partitions = {}.\n"
        "The number of cells in an axis cannot be lower that the partition count of this axis\n";

      int nbPartitions = 1;
      for( int i = 0; i < m_nsdof; i++ )
      {
        nbPartitions *= this->m_partitions( i );
      }
      GEOS_ERROR_IF_NE_MSG( nbPartitions, m_size, GEOS_FMT( partitionsLogMessage, m_size, nbPartitions )  );
    }

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOS, m_nsdof, m_partitions.data(), m_periodic.data(), reorder, &cartcomm );
    }
    m_rank = MpiWrapper::commRank( cartcomm );
    MpiWrapper::cartCoords( cartcomm, m_rank, m_nsdof, m_coords.data());

    //add neighbors
    {
      int ncoords[m_nsdof];
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
  for( int i = 0; i < m_nsdof; ++i )
  {
    const int nloc = m_partitions( i ) - 1;
    const localIndex nlocl = static_cast< localIndex >(nloc);
    if( m_partitionLocations[i].empty() )
    {
      // the default "even" spacing
      m_blockSize[ i ] /= m_partitions( i );
      m_min[ i ] += m_coords( i ) * m_blockSize[ i ];
      m_max[ i ] = min[ i ] + (m_coords( i ) + 1) * m_blockSize[ i ];

      m_partitionLocations[i].resize( nlocl );
      for( localIndex j = 0; j < m_partitionLocations[ i ].size(); ++j )
      {
        m_partitionLocations[ i ][ j ] = (j+1) * m_blockSize[ i ];
      }
    }
    else if( nlocl == m_partitionLocations[i].size() )
    {
      const int parIndex = m_coords[i];
      if( parIndex == 0 )
      {
        m_min[i] = min[i];
        m_max[i] = m_partitionLocations[i][parIndex];
      }
      else if( parIndex == nloc )
      {
        m_min[i] = m_partitionLocations[i][parIndex-1];
        m_max[i] = max[i];
      }
      else
      {
        m_min[i] = m_partitionLocations[i][parIndex-1];
        m_max[i] = m_partitionLocations[i][parIndex];
      }
    }
    else
    {
      GEOS_ERROR( "SpatialPartition::setSizes(): number of partition locations does not equal number of partitions - 1\n" );
    }
  }
}

void SpatialPartition::updateSizes( arrayView1d< real64 > const domainL,
                                    real64 const dt )
{
  for( int i=0; i<3; i++ )
  {
    real64 ratio = 1.0 + domainL[i] * dt;
    m_min[i] *= ratio;
    m_max[i] *= ratio;
    //m_partitionLocations[i] *= ratio; ?
    m_blockSize[i] *= ratio;
    m_gridSize[i] *= ratio;
    m_gridMin[i] *= ratio;
    m_gridMax[i] *= ratio;
    m_contactGhostMin[i] *= ratio;
    m_contactGhostMax[i] *= ratio;
  }
}

bool SpatialPartition::isCoordInPartition( const real64 & coord, const int dir ) const
{
  bool rval = true;
  const int i = dir;
  if( m_periodic( i ) )
  {
    if( m_partitions( i ) != 1 )
    {
      real64 localCenter = MapValueToRange( coord, m_gridMin[ i ], m_gridMax[ i ] );
      rval = rval && localCenter >= m_min[ i ] && localCenter < m_max[ i ];
    }

  }
  else
  {
    rval = rval && (m_partitions[ i ] == 1 || (coord >= m_min[ i ] && coord < m_max[ i ]));
  }

  return rval;
}

bool SpatialPartition::isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                                      const real64 & boundaryRadius ) const
// test a point relative to a boundary box. If non-zero buffer specified, expand the box.
{
  for( int i = 0; i < m_nsdof; i++ )
  {
    // Is particle already in bounds of partition?
    if( !(m_partitions( i )==1 || ( elemCenter[i] >= (m_min[i] - boundaryRadius) && elemCenter[i] <= (m_max[i] + boundaryRadius) ) ) )
    {
      // Particle not in bounds, check if direction has a periodic boundary
      if( m_periodic( i ) && (m_coords[i] == 0 || m_coords[i] == m_partitions[i] - 1) )
      {
        // Partition minimum boundary is periodic
        if( m_coords[i] == 0 && ( (elemCenter[i] - m_gridSize[i]) < (m_min[i] - boundaryRadius) ) )
        {
          return false;
        }
        // Partition maximum boundary is periodic
        if( m_coords[i] == m_partitions[i] - 1 && ( (elemCenter[i] + m_gridSize[i]) > (m_max[i] + boundaryRadius) ) )
        {
          return false;
        }
      }
      else
      {
        // No periodic boundary then particle is not in partition
        return false;
      }
    }
  }
  return true;
}

void SpatialPartition::setContactGhostRange( const real64 bufferSize )
{
  LvArray::tensorOps::copy< 3 >( m_contactGhostMin, m_min );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMin, -bufferSize );

  LvArray::tensorOps::copy< 3 >( m_contactGhostMax, m_max );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMax, bufferSize );
}

//CC: overrides global indices on periodic faces so they are matched when finding neighboring nodes
void SpatialPartition::setPeriodicDomainBoundaryObjects( MeshBody & grid,
                                                         NodeManager & nodeManager,
                                                         EdgeManager & edgeManager,
                                                         FaceManager & faceManager )
{
  GEOS_LOG_RANK( "Set periodic domain boundary objects" );
  arrayView1d< globalIndex > localToGlobalMap = nodeManager.localToGlobalMap();
  // unordered_map< globalIndex, localIndex > const & globalToLocalMap = nodeManager.globalToLocalMap(); // CC: need this for single
  // partition case
  const arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > gridPosition = nodeManager.referencePosition();

  // CC: Should we be using periodicSets? Old geos used periodic sets in the input file, we don't here
  CellBlockManager & cellBlockManager = grid.getGroup< CellBlockManager >( dataRepository::keys::cellManager );
  auto & nodeSets = cellBlockManager.getNodeSets();

  //Get cartesian communicator to get rank of neighbor periodic partition
  MPI_Comm cartcomm;
  {
    int reorder = 0;
    MpiWrapper::cartCreate( MPI_COMM_GEOS, 3, m_partitions.data(), m_periodic.data(), reorder, &cartcomm );
    GEOS_ERROR_IF( cartcomm == MPI_COMM_NULL, "Fail to run MPI_Cart_create and establish communications" );
  }

  // Check for periodic boundaries in each direction
  for( unsigned int dimension =0; dimension < 3; dimension++ )
  {
    if( m_periodic[dimension] )
    {
      // Is this partition on a boundary of domain?
      if( (m_coords[dimension] == 0)  ||
          (m_coords[dimension] == m_partitions[dimension]-1) )
      {
        // Reset global id numbers
        ///////////////////////////

        // Pick sets based on direction
        string setnames[2];
        switch( dimension )
        {
          case 0:
            setnames[0] = "xneg";
            setnames[1] = "xpos";
            break;
          case 1:
            setnames[0] = "yneg";
            setnames[1] = "ypos";
            break;
          case 2:
            setnames[0] = "zneg";
            setnames[1] = "zpos";
            break;
          default:
            GEOS_ERROR( "SpatialPartition::setPeriodicDomainBoundaryObjects() unrecognized direction!\n" );
        }

        SortedArray< localIndex > * theSets[2];
        theSets[0] = &(nodeSets[setnames[0]]);
        theSets[1] = &(nodeSets[setnames[1]]);

        PlanarSorter planarSorter( gridPosition, dimension );

        if( m_partitions[dimension] > 1 )
        {
          // Multiple partitions

          // Find periodic neighbor partition coordinates
          array1d< int > nbr_coords = m_coords;
          if( m_coords[dimension] == 0 )
          {
            nbr_coords[dimension] = m_partitions[dimension]-1;
          }
          else
          {
            nbr_coords[dimension] = 0;
          }

          int mySetId = (theSets[0]->size() > 0)? 0 : 1;
          int nbrSetId = 1-mySetId;
          if( theSets[nbrSetId]->size() > 0 )
          {
            GEOS_ERROR( "SpatialPartition::SetPeriodicDomainBoundaryObjects: " + setnames[0] + " and " + setnames[1] + " present on same partition\n" );
          }
          SortedArray< localIndex > & mySet =  *(theSets[mySetId]);

          // gather local and global ids
          std::vector< std::pair< localIndex, localIndex > > myLocalAndGlobalIds;

          for( int i = 0; i < mySet.size(); ++i )
          {
            localIndex globalId = localToGlobalMap[ mySet[i] ];   //nodeGlobalIds[*itr];
            myLocalAndGlobalIds.push_back( std::pair< localIndex, localIndex >( mySet[i], globalId ) );
          }

          // Sort local/global ids by position in plane
          array1d< localIndex > mySortedGlobalIds( myLocalAndGlobalIds.size());
          array1d< localIndex > nbrSortedGlobalIds;

          std::sort( myLocalAndGlobalIds.begin(), myLocalAndGlobalIds.end(), planarSorter );
          for( unsigned int ii = 0; ii <myLocalAndGlobalIds.size(); ++ii )
          {
            mySortedGlobalIds[ii]  = myLocalAndGlobalIds[ii].second;
          }

          int neighbor_rank = MpiWrapper::cartRank( cartcomm, nbr_coords.data());  // Get rank of periodic neighbor
          int neighborsTag = 54;

          // Perform manual MPI communication with neighbor
          MPI_Request mpiRequest = MPI_REQUEST_NULL;
          MPI_Status mpiStatus;

          MpiWrapper::iSend( mySortedGlobalIds,
                             neighbor_rank,
                             neighborsTag,
                             MPI_COMM_GEOS,
                             &mpiRequest );

          MpiWrapper::recv( nbrSortedGlobalIds,
                            neighbor_rank,
                            neighborsTag,
                            MPI_COMM_GEOS,
                            &mpiStatus );

          MpiWrapper::waitAll( 1, &mpiRequest, &mpiStatus );  //does the count refer to siz

          // should have same number of nodes in both sets
          if( nbrSortedGlobalIds.size() !=  mySortedGlobalIds.size() )
          {
            GEOS_ERROR( "SpatialPartition::SetPeriodicDomainBoundaryObjects: Size of " + setnames[mySetId] + " does not match size of " + setnames[nbrSetId] + " on neighboring partition\n" );
          }

          // assign new global ids
          for( unsigned int ii = 0; ii < myLocalAndGlobalIds.size(); ++ii )
          {
            localIndex & nd =  myLocalAndGlobalIds[ii].first;
            localToGlobalMap[nd] = std::min( mySortedGlobalIds[ii], nbrSortedGlobalIds[ii] );
            nodeManager.updateGlobalToLocalMap( nd ); //Update global to local map so it doesn't crash when matching domain boundary objects
          }

        }
        else
        {
          //CC: Logic for single partition periodic boundaries is unimplemented

          //          // Single partition
          //          //-----------------

          //          // Nodes
          //          {
          //            std::vector< std::vector<std::pair<localIndex, localIndex>  >  > setLocalAndGlobalIds(2);
          //            for(int a =0; a<2; ++a){
          //              // Gather local/global ids
          //              for( lSet::iterator itr=theSets[a]->begin() ; itr!=theSets[a]->end() ; ++itr )
          //              {
          //                localIndex globalId = nodeGlobalIds[*itr];
          //                setLocalAndGlobalIds[a].push_back(std::pair<localIndex , localIndex>( *itr,globalId) );
          //              }
          //              // Sort local/global ids by position in plane
          //              std::sort(setLocalAndGlobalIds[a].begin(),setLocalAndGlobalIds[a].end(),planarSorter);
          //            }

          //            // should have same number of nodes in both sets
          //            if(setLocalAndGlobalIds[0].size() !=  setLocalAndGlobalIds[1].size() )
          //            {
          //              throw GPException("SpatialPartition::SetPeriodicDomainBoundaryObjects: Size of " + setnames[0] + " does not match
          // size of " + setnames[1] + " on process " +toString(m_rank) +  "\n");
          //            }

          //            // assign new global ids and make global to local map point to nodes on min boundary
          //            for(unsigned int ii = 0 ; ii <setLocalAndGlobalIds[0].size() ; ++ii ){
          //              localIndex& nd0 =  setLocalAndGlobalIds[0][ii].first;
          //              localIndex& nd1 =  setLocalAndGlobalIds[1][ii].first;

          //              // this could be done once (all nodes in the same set should lie on the one boundary)
          //              int minBoundarySetIndx = 0;
          //              if(  (*domain.m_feNodeManager.m_refposition)[nd1][dimension] <
          // (*domain.m_feNodeManager.m_refposition)[nd0][dimension] ){
          //                minBoundarySetIndx = 1;
          //              }
          //              int maxBoundarySetIndx = 1 - minBoundarySetIndx;
          //              localIndex localTarget = (minBoundarySetIndx == 0)? nd0 : nd1;
          // //             localIndex notThelocalTarget = (minBoundarySetIndx == 0)? nd1 : nd0;

          //              // fix up local to global map
          //              localIndex minBoundGlobalId = setLocalAndGlobalIds[minBoundarySetIndx][ii].second;
          //              localIndex maxBoundGlobalId = setLocalAndGlobalIds[maxBoundarySetIndx][ii].second;

          //              nodeGlobalIds[nd0] = minBoundGlobalId;
          //              nodeGlobalIds[nd1] = minBoundGlobalId;

          //              // fix up global to local map
          //              nodeGlobalToLocalMap[minBoundGlobalId] = localTarget;

          //              // not used? in any case make old Global id point to same local target
          //              nodeGlobalToLocalMap[maxBoundGlobalId] = localTarget;
          //            }
          //          }

        }

        // CC: For periodic MPM stuff we don't need to update edge or face managers, right?
        // Only need to additively sync the values of the grid nodes
        for( int i = 0; i < theSets[0]->size(); ++i )
        {
          nodeManager.getDomainBoundaryIndicator()[(*theSets[0])[i]] = 1;
          edgeManager.getDomainBoundaryIndicator()[(*theSets[0])[i]] = 1; // CC: Do I need to do this since we only use nodes?
          faceManager.getDomainBoundaryIndicator()[(*theSets[0])[i]] = 1; // CC: Do I need to do this since we only use nodes?
        }

        for( int i = 0; i < theSets[1]->size(); ++i )
        {
          nodeManager.getDomainBoundaryIndicator()[(*theSets[1])[i]] = 1;
          edgeManager.getDomainBoundaryIndicator()[(*theSets[1])[i]] = 1; // CC: Do I need to do this since we only use nodes?
          faceManager.getDomainBoundaryIndicator()[(*theSets[1])[i]] = 1; // CC: Do I need to do this since we only use nodes?
        }
      }
    }
  }

  MpiWrapper::commFree( cartcomm );
}


void SpatialPartition::repartitionMasterParticles( DomainPartition & domain,
                                                   ParticleSubRegion & subRegion )
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
   * Rank of particles that were moved from the current partition will be correct,
   * but the master-ghost map in the neighbor list still needs to updated.  Particles that
   * were ghosts of an object that has been repartitioned will need to have their ghost
   * rank updated, and may need to be deleted/added elsewhere.
   *
   */

  // (1) Identify any particles that are master on the current partition, but whose center lies
  // outside of the partition domain.  Rank() for particles is defined such that it always
  // equals the rank of the owning process. Thus a particle is master iff Rank==partition.rank
  //
  // Temporarily set the ghost rank of any particle to be moved to "-1".  If the particle is
  // requested by another partition, its ghost rank will be updated.  Any particle that still
  // has a Rank=-1 at the end of this function is lost and needs to be deleted.  This
  // should only happen if it has left the global domain (hopefully at an outflow b.c.).

  MPI_iCommData commData;
  commData.resize( domain.getNeighbors().size() );

  arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
  arrayView1d< localIndex > const particleRank = subRegion.getParticleRank();
  array1d< R1Tensor > outOfDomainParticleCoordinates;
  stdVector< localIndex > outOfDomainParticleLocalIndices;
  unsigned int nn = m_neighbors.size();   // Number of partition neighbors.

  forAll< serialPolicy >( subRegion.size(), [&, particleCenter, particleRank] GEOS_HOST ( localIndex const pp )
    {
      bool inPartition = true;
      R1Tensor p_x;
      for( int i=0; i<3; i++ )
      {
        p_x[i] = particleCenter[pp][i];
        inPartition = inPartition && isCoordInPartition( p_x[i], i );
      }
      if( particleRank[pp]==this->m_rank && !inPartition )
      {
        outOfDomainParticleCoordinates.emplace_back( p_x ); // Store the coordinate of the out-of-domain particle
        outOfDomainParticleLocalIndices.push_back( pp );   // Store the local index "pp" for the current coordinate.
        particleRank[pp] = -1;                             // Temporarily set particleRank of out-of-domain particle to -1 until it is
                                                           // requested by someone.
      }
    } );


  // (2) Pack the list of particle center coordinates to each neighbor, and send/receive the list to neighbors.

  stdVector< array1d< R1Tensor > > particleCoordinatesReceivedFromNeighbors( nn );

  sendCoordinateListToNeighbors( outOfDomainParticleCoordinates.toView(),       // input: Single list of coordinates sent to all neighbors
                                 commData,                                      // input: Solver MPI communicator
                                 particleCoordinatesReceivedFromNeighbors );    // output: List of lists of coordinates received from each
                                                                                // neighbor.



  // (3) check the received lists for particles that are in the domain of the
  //     current partition.  make a list of the locations in the coordinate list
  //     of the particles that are to be owned by the current partition.

  stdVector< array1d< localIndex > > particleListIndicesRequestingFromNeighbors( nn );
  for( size_t n=0; n<nn; n++ )
  {
    // Loop through the unpacked list and make a list of the index of any point in partition interior domain
    for( int pp=0; pp<particleCoordinatesReceivedFromNeighbors[n].toView().size(); pp++ )
    {
      bool inPartition = true;
      for( int j=0; j<3; j++ )
      {
        inPartition = inPartition && isCoordInPartition( particleCoordinatesReceivedFromNeighbors[n].toView()[pp][j], j );
      }
      if( inPartition )
      {
        // Request particle to be transferred, and take ownership
        particleListIndicesRequestingFromNeighbors[n].emplace_back( pp );
      }
    }
  }


  // (4) Pack and send/receive list of requested indices to each neighbor.  These are the indices
  //     in the list of coordinates, not the LocalIndices on the sending processor. Unpack it
  //     and store the request list.

  stdVector< array1d< localIndex > > particleListIndicesRequestedFromNeighbors( nn );

  sendListOfIndicesToNeighbors< localIndex >( particleListIndicesRequestingFromNeighbors,
                                              commData,
                                              particleListIndicesRequestedFromNeighbors );


  // (5) Update the ghost rank of the out-of-domain particles to be equal to the rank
  //     of the partition requesting to own the particle.

  stdVector< array1d< localIndex > > particleLocalIndicesRequestedFromNeighbors( nn );
  {
    unsigned int numberOfRequestedParticles = 0;
    stdVector< int > outOfDomainParticleRequests( outOfDomainParticleLocalIndices.size(), 0 );

    for( size_t n=0; n<nn; n++ )
    {
      int ni = particleListIndicesRequestedFromNeighbors[n].size();
      numberOfRequestedParticles += ni;

      // The corresponding local index for each item in the request list is stored here:
      particleLocalIndicesRequestedFromNeighbors[n].resize( ni );

      forAll< serialPolicy >( ni, [&, particleRank] GEOS_HOST ( localIndex const k )
        {
          int const i = particleListIndicesRequestedFromNeighbors[n][k];
          outOfDomainParticleRequests[i] += 1;
          localIndex pp = outOfDomainParticleLocalIndices[i];

          particleLocalIndicesRequestedFromNeighbors[n][k] = pp;
          // Set ghost rank of the particle equal to neighbor rank.
          particleRank[pp] = m_neighbors[n].neighborRank();
        } );
    }

    // Check that there is exactly one processor requesting each out-of-domain particle.
    if( numberOfRequestedParticles != outOfDomainParticleLocalIndices.size())
    {
      std::cout << "Rank " << m_rank << " has requests for " << numberOfRequestedParticles << " out of " << outOfDomainParticleLocalIndices.size() << " out-of-domain particles" << std::endl;
    }
    for( size_t i=0; i<outOfDomainParticleRequests.size(); i++ )
    {
      if( outOfDomainParticleRequests[i] != 1 )
      {
        std::cout << "Rank " << m_rank << " particle as " << outOfDomainParticleRequests[i] << " != 1 requests!" << std::endl;
      }
    }
  }


  // (5.1) Resize particle subRegion to accommodate incoming particles.
  //       Keep track of the starting indices and number of particles coming from each neighbor.

  int oldSize = subRegion.size();
  int newSize = subRegion.size();
  stdVector< int > newParticleStartingIndices( nn );
  stdVector< int > numberOfIncomingParticles( nn );
  for( size_t n=0; n<nn; n++ )
  {
    numberOfIncomingParticles[n] = particleListIndicesRequestingFromNeighbors[n].size();
    newParticleStartingIndices[n] = newSize;
    newSize += numberOfIncomingParticles[n];
  }
  if( newSize > oldSize )
  {
    subRegion.resize( newSize ); // TODO: Does this handle constitutive fields owned by the subRegion?
  }


  // (6) Pack a buffer for the particles to be sent to each neighbor, and send/receive

  //int sizeBeforeParticleSend = subRegion.size(); // subregion size changes after this, so we need this here to use to size the deletion
  // loop
  sendParticlesToNeighbor( subRegion,
                           newParticleStartingIndices,
                           numberOfIncomingParticles,
                           commData,
                           particleLocalIndicesRequestedFromNeighbors );


  // (7) Delete any out-of-domain particles that were not requested by a neighbor.  These particles
  //     will still have Rank=-1. This should only happen if the particle has left the global domain.
  //     which will hopefully only occur at outflow boundary conditions.  If it happens for a particle in
  //     the global domain, print a warning.

  arrayView2d< real64 > const particleCenterAfter = subRegion.getParticleCenter();
  arrayView1d< int > const particleRankAfter = subRegion.getParticleRank();
  std::set< localIndex > indicesToErase;
  forAll< serialPolicy >( subRegion.size(), [&, particleRankAfter, particleCenterAfter] GEOS_HOST ( localIndex const p )
    {
      if( particleRankAfter[p] == -1 )
      {
        GEOS_LOG_RANK( "Deleting orphan out-of-domain particle during repartition at p_x = " << particleCenterAfter[p] );
        indicesToErase.insert( p );
      }
      else if( particleRankAfter[p] != m_rank )
      {
        indicesToErase.insert( p );
      }
    } );
  subRegion.erase( indicesToErase );

  // Resize particle region owning this subregion
  ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
  int newRegionSize = region.getNumberOfParticles();
  region.resize( newRegionSize );
}


// void SpatialPartition::getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
//                                                                    const real64 & boundaryRadius )
// {

//   /*
//    * Make a list of the coordinates and global IDs of all non-ghost objects on the current
//    * partition.  These should all be interior to the partition domain (excluding the ghost
//    * region).  Send this list to each neighbor.  Have each neighbor check the list for
//    * objects within its bounding box (including ghost radius) that do not already exist
//    * on the processor.  Send a request list for those objects.  Mark all ghost objects as
//    * potentially abandoned.
//    *
//    */

//   MPI_iCommData commData;
//   commData.resize( domain.getNeighbors().size() );

//   // MPM-specific code where we assume there are 2 mesh bodies and only one of them has particles
//   dataRepository::Group & meshBodies = domain.getMeshBodies();
//   MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >( 0 );
//   MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >( 1 );
//   MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
//   ParticleManager & particleManager = particles.getBaseDiscretization().getParticleManager();

//   particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
//   {
//     // (1) Identify any particles that are master on the current partition, but whose center lies
//     // outside of the partition domain.  Rank() for particles is defined such that it always
//     // equals the rank of the owning process. Thus a particle is master iff Rank==partition.rank
//     //
//     // Temporarily set the ghost rank of any particle to be moved to "-1".  If the particle is
//     // requested by another partition, its ghost rank will be updated.  Any particle that still
//     // has a Rank=-1 at the end of this function is lost and needs to be deleted.  This
//     // should only happen if it has left the global domain (hopefully at an outflow b.c.).

//     arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
//     arrayView1d< localIndex > const particleRank = subRegion.getParticleRank();
//     arrayView1d< globalIndex > const particleGlobalID = subRegion.getParticleID();
//     array1d< R1Tensor > inDomainMasterParticleCoordinates;   // Theoretically the same as particle position evaluated at
//                                                              // subRegion.nonGhostIndices()?
//     stdVector< globalIndex > inDomainMasterParticleGlobalIndices;
//     unsigned int nn = m_neighbors.size();   // Number of partition neighbors.

//     forAll< serialPolicy >( subRegion.size(), [&, particleCenter, particleRank, particleGlobalID] GEOS_HOST ( localIndex const p )
//       {
//         bool inPartition = true;
//         R1Tensor p_x;
//         for( int i=0; i<3; i++ )
//         {
//           p_x[i] = particleCenter[p][i];
//           inPartition = inPartition && isCoordInPartition( p_x[i], i );
//         }
//         if( particleRank[p]==this->m_rank && inPartition )
//         {
//           inDomainMasterParticleCoordinates.emplace_back( p_x );  // Store the coordinate of the out-of-domain particle
//           inDomainMasterParticleGlobalIndices.push_back( particleGlobalID[p] );     // Store the local index "pp" for the current
//                                                                                     // coordinate.
//         }
//       } );


//     // (2) Pack the list of particle center coordinates to each neighbor, and send/receive the list to neighbors.

//     stdVector< array1d< R1Tensor > > particleCoordinatesReceivedFromNeighbors( nn );

//     sendCoordinateListToNeighbors( inDomainMasterParticleCoordinates.toView(),          // input: Single list of coordinates sent to all
//                                                                                         // neighbors
//                                    commData,                                      // input: Solver MPI communicator
//                                    particleCoordinatesReceivedFromNeighbors    // output: List of lists of coordinates received from each
//                                                                                // neighbor.
//                                    );


//     // (3) Pack the list of particle global indices to each neighbor, and send the list to neighbors.

//     stdVector< array1d< globalIndex > > particleGlobalIndicesSendingToNeighbors( nn );
//     stdVector< array1d< globalIndex > > particleGlobalIndicesReceivedFromNeighbors( nn );

//     for( size_t n=0; n<nn; ++n )
//     {
//       particleGlobalIndicesSendingToNeighbors[n].resize( inDomainMasterParticleGlobalIndices.size() );
//       for( size_t i=0; i<inDomainMasterParticleGlobalIndices.size(); ++i )
//       {
//         particleGlobalIndicesSendingToNeighbors[n][i] = inDomainMasterParticleGlobalIndices[i];
//       }
//     }
//     sendListOfIndicesToNeighbors< globalIndex >( particleGlobalIndicesSendingToNeighbors,
//                                                  commData,
//                                                  particleGlobalIndicesReceivedFromNeighbors );


//     // (4) Look through the received coordinates and make a list of the particles that are within
//     //     the bounding box (including ghost region) and for which the globalID does not already
//     //     exist on the current partition.  Make a request list of the index (order in the list)
//     //     of the particle.  This will be sent from the master as a new ghost on the current
//     //     partition.

//     stdVector< array1d< globalIndex > > particleGlobalIndicesRequestingFromNeighbors( nn );

//     for( size_t n=0; n<nn; ++n )
//     {
//       for( localIndex i=0; i<particleCoordinatesReceivedFromNeighbors[n].size(); ++i )
//       {
//         if( isCoordInPartitionBoundingBox( particleCoordinatesReceivedFromNeighbors[n][i], boundaryRadius ) )
//         {
//           globalIndex gI = particleGlobalIndicesReceivedFromNeighbors[n][i];

//           // This particle should be a ghost on the current processor. See if it already exists here.
//           bool alreadyHere = false;
//           forAll< serialPolicy >( subRegion.size(), [&, particleGlobalID, particleRank] GEOS_HOST ( localIndex const p )
//             {
//               if( gI == particleGlobalID[p] )
//               {
//                 // The particle already exists as a ghost on this partition, so we should update its rank.
//                 particleRank[p] = m_neighbors[n].neighborRank();
//                 alreadyHere = true;
//               }
//             } );
//           if( !alreadyHere )
//           {
//             // The global index is not represented on this partition, so we should add the particle.
//             particleGlobalIndicesRequestingFromNeighbors[n].emplace_back( gI );
//           }
//         }
//       }
//     }


//     // (5) Pack and send request list of Global Indices to each neighbor, receive and unpack
//     //     this into a list requested from each neighbor.

//     stdVector< array1d< globalIndex > > particleGlobalIndicesRequestedFromNeighbors( nn );

//     sendListOfIndicesToNeighbors< globalIndex >( particleGlobalIndicesRequestingFromNeighbors,
//                                                  commData,
//                                                  particleGlobalIndicesRequestedFromNeighbors );


//     // (6) Temporarily set the ghost rank of all ghosts to "-1".  After ghosts are unpacked from the
//     //     masters, the ghost rank will be overwritten.  At the end of this function, any ghosts that
//     //     still have ghostRank=-1 are orphans and need to be deleted.

//     int partitionRank = this->m_rank;
//     forAll< parallelHostPolicy >( subRegion.size(), [=] GEOS_HOST ( localIndex const p )   // TODO: Worth moving to device?
//       {
//         if( particleRank[p] != partitionRank )
//         {
//           particleRank[p] = -1;
//         }
//       } );


//     // (6.1) Resize particle subRegion to accommodate incoming particles.
//     //       Keep track of the starting indices and number of particles coming from each neighbor.

//     int oldSize = subRegion.size();
//     int newSize = subRegion.size();
//     stdVector< int > newParticleStartingIndices( nn );
//     stdVector< int > numberOfIncomingParticles( nn );
//     for( size_t n=0; n<nn; n++ )
//     {
//       numberOfIncomingParticles[n] = particleGlobalIndicesRequestingFromNeighbors[n].size();
//       newParticleStartingIndices[n] = newSize;
//       newSize += numberOfIncomingParticles[n];
//     }
//     if( newSize > oldSize )
//     {
//       subRegion.resize( newSize );   // TODO: Does this handle constitutive fields owned by the subregion's parent region?
//     }


//     // (7) Pack/Send/Receive/Unpack particles to be sent to each neighbor.

//     {
//       stdVector< array1d< localIndex > > particleLocalIndicesRequestedFromNeighbors( nn );

//       for( size_t n=0; n<nn; n++ )
//       {
//         // make a list of the local indices corresponding to the global indices in the request list for the current neighbor.
//         particleLocalIndicesRequestedFromNeighbors[n].resize( particleGlobalIndicesRequestedFromNeighbors[n].size() );
//         for( localIndex i=0; i<particleLocalIndicesRequestedFromNeighbors[n].size(); ++i )
//         {
//           particleLocalIndicesRequestedFromNeighbors[n][i] = subRegion.globalToLocalMap( particleGlobalIndicesRequestedFromNeighbors[n][i] );
//         }
//       }
//       // Send/receive the particles, which will add missing ghosts on the partition.
//       sendParticlesToNeighbor( subRegion,
//                                newParticleStartingIndices,
//                                numberOfIncomingParticles,
//                                commData,
//                                particleLocalIndicesRequestedFromNeighbors );
//     }


//     // (8) Delete any particles that have ghostRank=-1.  These will be ghosts from
//     //     a previous step for which the master is no longer in the ghost domain,
//     // std::set< localIndex > indicesToErase;
//     // arrayView1d< localIndex > const particleRankNew = subRegion.getParticleRank();
//     // forAll< serialPolicy >( subRegion.size(), [=, &indicesToErase] GEOS_HOST ( localIndex const p ) // parallelize with atomics
//     // {
//     //   if( particleRankNew[p] == -1 )
//     //   {
//     //     indicesToErase.insert(p);
//     //   }
//     // } );
//     // subRegion.erase(indicesToErase);

//   } );
// }

void SpatialPartition::getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                                   const real64 & boundaryRadius )
{
  /*
   * Faster MPM ghost particle refresh.
   *
   * Main differences from the previous implementation:
   *
   *  1. Send only boundary-near master particles to the neighbors that can
   *     actually need them as ghosts. The old implementation sent every owned
   *     master particle to every neighbor.
   *
   *  2. Send coordinate, globalID, and owner-local-index together as one compact
   *     metadata record. This removes a separate globalID metadata exchange and
   *     avoids owner-side globalID -> localID lookup when servicing requests.
   *
   *  3. Build a globalID -> localIndex lookup map once. The old implementation
   *     scanned the entire subRegion for every received candidate.
   *
   *  4. Invalidate old ghosts before processing candidates. The old ordering
   *     refreshed existing ghost ranks and then overwrote them with -1.
   *
   *  5. Delete stale ghosts whose rank remains -1 after the refresh.
   */

  MPI_iCommData commData;
  commData.resize( domain.getNeighbors().size() );

  unsigned int const nn = m_neighbors.size();

  /*
   * Build neighbor offset records in the same traversal order as addNeighbors().
   *
   * This lets us know whether m_neighbors[n] is the -x face neighbor, +x/+y edge
   * neighbor, +x/+y/+z corner neighbor, etc., without changing NeighborCommunicator
   * or SpatialPartition.hpp.
   *
   * Important: we intentionally preserve duplicate neighbor entries in periodic
   * small decompositions because m_neighbors itself may contain those duplicates
   * and the existing communication helpers expect one list per m_neighbors entry.
   */
  stdVector< std::array< int, 3 > > neighborOffsets;
  neighborOffsets.reserve( nn );

  {
    int ncoords[3] = { 0, 0, 0 };
    int offsets[3] = { 0, 0, 0 };

    std::function< void( int ) > appendNeighborOffset;

    appendNeighborOffset = [&]( int const idim )
    {
      if( idim == m_nsdof )
      {
        bool me = true;
        for( int i = 0; i < m_nsdof; ++i )
        {
          if( ncoords[i] != this->m_coords( i ) )
          {
            me = false;
            break;
          }
        }

        if( !me )
        {
          neighborOffsets.emplace_back( std::array< int, 3 >{ { offsets[0], offsets[1], offsets[2] } } );
        }

        return;
      }

      localIndex const localDim = LvArray::integerConversion< localIndex >( idim );

      int const dim = this->m_partitions( localDim );
      bool const periodic = this->m_periodic( localDim );

      for( int delta = -1; delta <= 1; ++delta )
      {
        offsets[idim] = delta;
        ncoords[idim] = this->m_coords( localDim ) + delta;

        bool ok = true;

        if( periodic )
        {
          if( ncoords[idim] < 0 )
          {
            ncoords[idim] = dim - 1;
          }
          else if( ncoords[idim] >= dim )
          {
            ncoords[idim] = 0;
          }
        }
        else
        {
          ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
        }

        if( ok )
        {
          appendNeighborOffset( idim + 1 );
        }
      }
    };

    appendNeighborOffset( 0 );
  }

  GEOS_ERROR_IF( neighborOffsets.size() != nn,
                 GEOS_FMT( "SpatialPartition::getGhostParticlesFromNeighboringPartitions: "
                           "neighbor offset count {} does not match m_neighbors count {}.",
                           neighborOffsets.size(), nn ) );

  /*
   * Use the particle manager from whichever mesh body owns particles.
   *
   * This keeps the existing MPM-specific assumption in the original code:
   * there are two mesh bodies, and only one has particles.
   */
  dataRepository::Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >( 0 );
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >( 1 );
  MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
  ParticleManager & particleManager = particles.getBaseDiscretization().getParticleManager();

  /*
   * Coordinate used for boundary-slab tests.
   *
   * For periodic dimensions, a particle coordinate may be stored outside the
   * nominal [gridMin, gridMax) interval. isCoordInPartition() already maps for
   * ownership tests. The boundary-slab filter must do the same, otherwise a
   * periodic particle near x = gridMin could be missed if it is stored as
   * x = gridMax + epsilon, or vice versa.
   */
  auto mapCoordinateForPartitionTests =
    [&]( real64 const x, int const d ) -> real64
  {
    if( this->m_periodic[d] && this->m_partitions[d] != 1 )
    {
      return MapValueToRange( x, this->m_gridMin[d], this->m_gridMax[d] );
    }

    return x;
  };

  /*
   * Decide whether an owned particle can possibly be needed by a specific neighbor.
   *
   * For a face neighbor, this is a slab test.
   * For an edge neighbor, this is an intersection of two slabs.
   * For a corner neighbor, this is an intersection of three slabs.
   *
   * Example:
   *   offset = {+1, 0, 0}
   *     send only particles with x >= m_max[0] - boundaryRadius.
   *
   *   offset = {-1, +1, 0}
   *     send only particles with
   *       x <= m_min[0] + boundaryRadius and
   *       y >= m_max[1] - boundaryRadius.
   */
  auto particleIsCandidateForNeighbor =
    [&]( R1Tensor const & mappedCenter,
         std::array< int, 3 > const & neighborOffset ) -> bool
  {
    for( int d = 0; d < 3; ++d )
    {
      if( this->m_partitions[d] == 1 || neighborOffset[d] == 0 )
      {
        continue;
      }

      if( neighborOffset[d] < 0 )
      {
        if( mappedCenter[d] > this->m_min[d] + boundaryRadius )
        {
          return false;
        }
      }
      else
      {
        if( mappedCenter[d] < this->m_max[d] - boundaryRadius )
        {
          return false;
        }
      }
    }

    return true;
  };

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView1d< localIndex > const particleRank = subRegion.getParticleRank();
    arrayView1d< globalIndex > const particleGlobalID = subRegion.getParticleID();

    int const partitionRank = this->m_rank;

    /*
     * Build one metadata list per neighbor.
     *
     * This replaces the old all-owned-particles-to-all-neighbors exchange.
     */
    stdVector< array1d< GhostParticleCandidate > > particleCandidatesSendingToNeighbors( nn );

    forAll< serialPolicy >( subRegion.size(), [&, particleCenter, particleRank, particleGlobalID] GEOS_HOST ( localIndex const p )
    {
      if( particleRank[p] != partitionRank )
      {
        return;
      }

      bool inPartition = true;

      R1Tensor rawCenter;
      R1Tensor mappedCenter;

      for( int d = 0; d < 3; ++d )
      {
        rawCenter[d] = particleCenter[p][d];
        mappedCenter[d] = mapCoordinateForPartitionTests( rawCenter[d], d );

        // Keep the original ownership check. If repartitioning invariants are
        // made stronger later, this could become a debug-only check.
        inPartition = inPartition && this->isCoordInPartition( rawCenter[d], d );
      }

      if( !inPartition )
      {
        return;
      }

      for( size_t n = 0; n < nn; ++n )
      {
        if( particleIsCandidateForNeighbor( mappedCenter, neighborOffsets[n] ) )
        {
          GhostParticleCandidate candidate;

          candidate.x[0] = rawCenter[0];
          candidate.x[1] = rawCenter[1];
          candidate.x[2] = rawCenter[2];

          candidate.globalID = particleGlobalID[p];

          // This local index is meaningful only on the owning rank and only for
          // this exchange. The receiver sends it back if it needs the full
          // particle payload.
          candidate.ownerLocalIndex = p;

          particleCandidatesSendingToNeighbors[n].emplace_back( candidate );
        }
      }
    } );

    /*
     * Exchange compact candidate metadata.
     *
     * Each received candidate is a master particle owned by the corresponding
     * neighbor. This rank will request only the candidates that lie in its ghost
     * box and are not already present locally.
     */
    stdVector< array1d< GhostParticleCandidate > > particleCandidatesReceivedFromNeighbors( nn );

    exchangePODListsWithNeighbors( particleCandidatesSendingToNeighbors,
                                   commData,
                                   particleCandidatesReceivedFromNeighbors,
                                   m_neighbors );

    /*
     * Build a globalID -> localIndex map once.
     *
     * The old code scanned subRegion.size() for every received candidate, which
     * is typically the dominant CPU cost.
     *
     * Insert masters first so that, in the unlikely event of a duplicated
     * globalID, a local master takes precedence over a ghost.
     */
    std::unordered_map< globalIndex, localIndex > localIndexByGlobalID;
    localIndexByGlobalID.reserve( static_cast< size_t >( subRegion.size() ) * 2 );

    for( localIndex p = 0; p < subRegion.size(); ++p )
    {
      if( particleRank[p] == partitionRank )
      {
        localIndexByGlobalID.emplace( particleGlobalID[p], p );
      }
    }

    for( localIndex p = 0; p < subRegion.size(); ++p )
    {
      if( particleRank[p] != partitionRank )
      {
        localIndexByGlobalID.emplace( particleGlobalID[p], p );
      }
    }

    /*
     * Mark all existing ghosts stale before processing received candidates.
     *
     * If a ghost is still needed, the candidate-processing loop below will set
     * its rank back to the owning neighbor rank.
     */
    forAll< parallelHostPolicy >( subRegion.size(), [=] GEOS_HOST ( localIndex const p )
    {
      if( particleRank[p] != partitionRank )
      {
        particleRank[p] = -1;
      }
    } );

    /*
     * For each neighbor, request only candidates that:
     *  - lie inside this rank's ghost bounding box, and
     *  - are not already present on this rank.
     *
     * The request is the owner-local-index, not the globalID. That lets the
     * owner pack particles directly without a globalToLocalMap lookup.
     */
    stdVector< array1d< localIndex > > particleLocalIndicesRequestingFromNeighbors( nn );

    std::unordered_set< globalIndex > pendingNewGhostGlobalIDs;

    {
      size_t totalReceivedCandidates = 0;
      for( size_t n = 0; n < nn; ++n )
      {
        totalReceivedCandidates += static_cast< size_t >( particleCandidatesReceivedFromNeighbors[n].size() );
      }
      pendingNewGhostGlobalIDs.reserve( totalReceivedCandidates );
    }

    for( size_t n = 0; n < nn; ++n )
    {
      int const ownerRank = m_neighbors[n].neighborRank();

      for( localIndex i = 0; i < particleCandidatesReceivedFromNeighbors[n].size(); ++i )
      {
        GhostParticleCandidate const & candidate = particleCandidatesReceivedFromNeighbors[n][i];

        R1Tensor candidateCenter;
        candidateCenter[0] = candidate.x[0];
        candidateCenter[1] = candidate.x[1];
        candidateCenter[2] = candidate.x[2];

        if( !this->isCoordInPartitionBoundingBox( candidateCenter, boundaryRadius ) )
        {
          continue;
        }

        globalIndex const gI = candidate.globalID;

        auto const iter = localIndexByGlobalID.find( gI );

        if( iter != localIndexByGlobalID.end() )
        {
          localIndex const p = iter->second;

          // Refresh existing ghosts. Never convert a local master into a ghost.
          if( particleRank[p] != partitionRank )
          {
            particleRank[p] = ownerRank;
          }
        }
        else
        {
          /*
           * Avoid requesting the same missing ghost more than once. This can
           * happen with duplicate periodic neighbor entries or when a particle
           * appears through more than one geometric neighbor relation.
           */
          if( pendingNewGhostGlobalIDs.insert( gI ).second )
          {
            particleLocalIndicesRequestingFromNeighbors[n].emplace_back( candidate.ownerLocalIndex );
          }
        }
      }
    }

    /*
     * Exchange request lists.
     *
     * After this call:
     *  - particleLocalIndicesRequestingFromNeighbors[n] are particles this rank
     *    wants to receive from neighbor n.
     *  - particleLocalIndicesRequestedFromNeighbors[n] are local master
     *    particles this rank must send to neighbor n.
     */
    stdVector< array1d< localIndex > > particleLocalIndicesRequestedFromNeighbors( nn );

    sendListOfIndicesToNeighbors< localIndex >( particleLocalIndicesRequestingFromNeighbors,
                                                commData,
                                                particleLocalIndicesRequestedFromNeighbors );

    /*
     * Resize this subRegion to hold incoming ghosts.
     */
    int const oldSize = subRegion.size();
    int newSize = subRegion.size();

    stdVector< int > newParticleStartingIndices( nn );
    stdVector< int > numberOfIncomingParticles( nn );

    for( size_t n = 0; n < nn; ++n )
    {
      numberOfIncomingParticles[n] = particleLocalIndicesRequestingFromNeighbors[n].size();
      newParticleStartingIndices[n] = newSize;
      newSize += numberOfIncomingParticles[n];
    }

    if( newSize > oldSize )
    {
      subRegion.resize( newSize );
    }

    /*
     * Send/receive the full particle payloads.
     *
     * The outgoing local-index lists are already owner-local indices, so no
     * owner-side globalToLocalMap lookup is needed here.
     */
    sendParticlesToNeighbor( subRegion,
                             newParticleStartingIndices,
                             numberOfIncomingParticles,
                             commData,
                             particleLocalIndicesRequestedFromNeighbors );

    /*
     * Delete ghosts that were not refreshed by this exchange.
     *
     * These are ghosts from previous steps whose masters are no longer inside
     * this rank's ghost region.
     */
    arrayView1d< localIndex > const particleRankAfter = subRegion.getParticleRank();

    std::set< localIndex > indicesToErase;

    forAll< serialPolicy >( subRegion.size(), [&, particleRankAfter] GEOS_HOST ( localIndex const p )
    {
      if( particleRankAfter[p] == -1 )
      {
        indicesToErase.insert( p );
      }
    } );

    bool const subRegionSizeChanged = newSize != oldSize || !indicesToErase.empty();

    if( !indicesToErase.empty() )
    {
      subRegion.erase( indicesToErase );
    }

    /*
     * Keep region-level particle fields consistent with the subRegion size,
     * mirroring the pattern used in repartitionMasterParticles().
     */
    if( subRegionSizeChanged )
    {
      ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
      int const newRegionSize = region.getNumberOfParticles();
      region.resize( newRegionSize );
    }
  } );
}

/**
 * @brief Send coordinates to neighbors as part of repartition.
 * @param[in] particleCoordinatesSendingToNeighbors Single list of coordinates sent to all neighbors
 * @param[in] commData Solver's MPI communicator
 * @param[in] particleCoordinatesReceivedFromNeighbors List of lists of coordinates received from each neighbor
 */
void SpatialPartition::sendCoordinateListToNeighbors( arrayView1d< R1Tensor > const & particleCoordinatesSendingToNeighbors,
                                                      MPI_iCommData & commData,
                                                      stdVector< array1d< R1Tensor > > & particleCoordinatesReceivedFromNeighbors
                                                      )
{
  // Number of neighboring partitions
  unsigned int nn = m_neighbors.size();

  // The send buffer is identical for each neighbor, and contains the packed coordinate list.
  unsigned int sizeToBePacked = 0;                                                        // size of the outgoing data with packing=false
                                                                                          // (we need to run through it first without
                                                                                          // packing so we can size the buffer)
  unsigned int sizeOfPacked = 0;                                                          // size of the outgoing data with packing=true
  buffer_unit_type * junk;                                                                // junk buffer, stores nothing since we're just
                                                                                          // getting the buffer size on the first pass
  sizeToBePacked += bufferOps::Pack< false >( junk,
                                              particleCoordinatesSendingToNeighbors );    // get the size of the coordinate list
  buffer_type sendBuffer( sizeToBePacked );                                               // the actual sized buffer that we pack into
  buffer_unit_type * sendBufferPtr = sendBuffer.data();                                   // get a pointer to the buffer
  sizeOfPacked += bufferOps::Pack< true >( sendBufferPtr,
                                           particleCoordinatesSendingToNeighbors );       // pack the coordinate list into the buffer
  GEOS_ERROR_IF_NE( sizeToBePacked, sizeOfPacked );                                       // make sure the packer is self-consistent

  // Declare the receive buffers
  stdVector< unsigned int > sizeOfReceived( nn );
  stdVector< buffer_type > receiveBuffer( nn );

  // send the coordinate list to each neighbor.  Using an asynchronous send,
  // the mpi request will be different for each send, but the buffer is the same
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    // Send/receive the size of the packed buffer
    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      m_neighbors[n].mpiISendReceive( &sizeOfPacked,
                                      1,
                                      sendRequest[n],
                                      &(sizeOfReceived[n]),
                                      1,
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Send/receive the buffer containing the list of coordinates
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      receiveBuffer[n].resize( sizeOfReceived[n] );
      m_neighbors[n].mpiISendReceive( sendBuffer.data(),   // TODO: This can't be sendBufferPtr, why not? I guess cuz sendBufferPtr gets
                                                           // incremented (moved) during packing.
                                      sizeOfPacked,
                                      sendRequest[n],
                                      receiveBuffer[n].data(),
                                      sizeOfReceived[n],
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Unpack the received coordinate list from each neighbor
  for( size_t n=0; n<nn; n++ )
  {
    // Unpack the buffer to an array of coordinates.
    const buffer_unit_type * receiveBufferPtr = receiveBuffer[n].data();  // needed for const cast
    bufferOps::Unpack( receiveBufferPtr, particleCoordinatesReceivedFromNeighbors[n] );
  }
}

template< typename indexType >
void SpatialPartition::sendListOfIndicesToNeighbors( stdVector< array1d< indexType > > & listSendingToEachNeighbor,
                                                     MPI_iCommData & commData,
                                                     stdVector< array1d< indexType > > & listReceivedFromEachNeighbor )
{
  // Number of neighboring partitions
  unsigned int nn = m_neighbors.size();

  // Pack the outgoing lists of local indices
  stdVector< unsigned int > sizeOfPacked( nn );
  stdVector< buffer_type > sendBuffer( nn );
  for( size_t n=0; n<nn; n++ )
  {
    unsigned int sizeToBePacked = 0;                                                  // size of the outgoing data with packing=false (we
                                                                                      // need to run through it first without packing so we
                                                                                      // can size the buffer)
    sizeOfPacked[n] = 0;                                                              // size of the outgoing data with packing=true
    buffer_unit_type * junk;                                                          // junk buffer, stores nothing since we're just
                                                                                      // getting the buffer size on the first pass
    sizeToBePacked += bufferOps::Pack< false >( junk,
                                                listSendingToEachNeighbor[n] );       // get the size of the list of local indices
    sendBuffer[n].resize( sizeToBePacked );                                           // the actual sized buffer that we pack into
    buffer_unit_type * sendBufferPtr = sendBuffer[n].data();                          // get a pointer to the buffer
    sizeOfPacked[n] += bufferOps::Pack< true >( sendBufferPtr,
                                                listSendingToEachNeighbor[n] );       // pack the list of local indices into the buffer
    GEOS_ERROR_IF_NE( sizeToBePacked, sizeOfPacked[n] );                             // make sure the packer is self-consistent
  }

  // Declare the receive buffers
  stdVector< unsigned int > sizeOfReceived( nn ); // TODO: decide if these number-of-neighbor-sized arrays should be array1d, stdVector
                                                  // or std::array
  stdVector< buffer_type > receiveBuffer( nn );

  // send the list of local indices to each neighbor using an asynchronous send
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    // Send/receive the size of the packed buffer
    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      m_neighbors[n].mpiISendReceive( &(sizeOfPacked[n]),
                                      1,
                                      sendRequest[n],
                                      &(sizeOfReceived[n]),
                                      1,
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Send/receive the buffer containing the list of local indices
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      receiveBuffer[n].resize( sizeOfReceived[n] );
      m_neighbors[n].mpiISendReceive( sendBuffer[n].data(),   // TODO: This can't be sendBufferPtr, why not? I guess cuz sendBufferPtr gets
                                                              // incremented (moved) during packing.
                                      sizeOfPacked[n],
                                      sendRequest[n],
                                      receiveBuffer[n].data(),
                                      sizeOfReceived[n],
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Unpack the received list of local indices from each neighbor
  for( size_t n=0; n<nn; n++ )
  {
    // Unpack the buffer to an array of coordinates.
    const buffer_unit_type * receiveBufferPtr = receiveBuffer[n].data();  // needed for const cast
    bufferOps::Unpack( receiveBufferPtr, listReceivedFromEachNeighbor[n] );
  }
}

void SpatialPartition::sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                                stdVector< int > const & newParticleStartingIndices,
                                                stdVector< int > const & numberOfIncomingParticles,
                                                MPI_iCommData & commData,
                                                stdVector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor )
{
  unsigned int nn = m_neighbors.size();

  // Pack the send buffer for the particles being sent to each neighbor
  stdVector< buffer_type > sendBuffer( nn );
  stdVector< unsigned int > sizeOfPacked( nn );

  for( size_t n=0; n<nn; n++ )
  {
    sizeOfPacked[n] = subRegion.particlePack( sendBuffer[n], particleLocalIndicesToSendToEachNeighbor[n].toView(), false );
    sendBuffer[n].resize( sizeOfPacked[n] );
    unsigned int sizeCheck = subRegion.particlePack( sendBuffer[n], particleLocalIndicesToSendToEachNeighbor[n].toView(), true );
    GEOS_ERROR_IF_NE( sizeCheck, sizeOfPacked[n] );
  }

  // Declare the receive buffers
  stdVector< unsigned int > sizeOfReceived( nn ); // TODO: decide if these number-of-neighbor-sized arrays should be array1d, stdVector
                                                  // or std::array
  stdVector< buffer_type > receiveBuffer( nn );

  // send/receive the size of the packed particle data to each neighbor using an asynchronous send
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      m_neighbors[n].mpiISendReceive( &(sizeOfPacked[n]),
                                      1,
                                      sendRequest[n],
                                      &(sizeOfReceived[n]),
                                      1,
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data() );
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data() );
  }

  // Send/receive the buffer containing the list of local indices
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      receiveBuffer[n].resize( sizeOfReceived[n] );
      m_neighbors[n].mpiISendReceive( sendBuffer[n].data(),
                                      sizeOfPacked[n],
                                      sendRequest[n],
                                      receiveBuffer[n].data(),
                                      sizeOfReceived[n],
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Unpack the received particle data.
  for( size_t n=0; n<nn; n++ )
  {
    // Unpack the buffer to an array of coordinates.
    subRegion.particleUnpack( receiveBuffer[n], newParticleStartingIndices[n], numberOfIncomingParticles[n] );
  }

}

REGISTER_CATALOG_ENTRY( PartitionBase, SpatialPartition, string const &, dataRepository::Group * const )
}
