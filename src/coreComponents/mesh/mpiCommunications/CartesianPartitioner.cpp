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
 * @file CartesianPartitioner.cpp
 */

#include "CartesianPartitioner.hpp"
#include "common/MpiWrapper.hpp"
#include "codingUtilities/Utilities.hpp"   // for isOdd
#include "LvArray/src/genericTensorOps.hpp"
#include <numeric>


namespace geos
{
using namespace dataRepository;


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

// Helper class to manage the MPI_Comm lifetime
class ScopedMpiComm
{
public:
  // Constructor is implicit
  // Destructor
  ~ScopedMpiComm()
  {
    if( m_comm != MPI_COMM_NULL )
    {
      MPI_Comm_free( &m_comm );
    }
  }

  // Allow the class to be used as an MPI_Comm handle
  operator MPI_Comm &() { return m_comm; }

private:
  MPI_Comm m_comm = MPI_COMM_NULL;
};



CartesianPartitioner::CartesianPartitioner( string const & name,
                                            Group * const parent ):
  PartitionerBase( name, parent ),
  m_periodic( m_ndim ),
  m_coords( m_ndim ),
  m_partitionCounts(),
  m_localMin{ 0.0, 0.0, 0.0 },
  m_localMax{ 0.0, 0.0, 0.0 },
  m_localSize{ 0.0, 0.0, 0.0 },
  m_globalGridMin{ 0.0, 0.0, 0.0 },
  m_globalGridMax{ 0.0, 0.0, 0.0 },
  m_globalGridSize{ 0.0, 0.0, 0.0 },
  m_cartRank( 0 )
{
  registerWrapper( viewKeyStruct::partitionCountsString(), &m_partitionCounts ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Method (library) used to partition the mesh" );
}

void CartesianPartitioner::partition( real64 const ( &globalGridMin )[ 3 ],
                                      real64 const ( &globalGridMax )[ 3 ] )
{
  ScopedMpiComm cartComm;
  initializeCartesianCommunicator( cartComm );
  determineNeighborsRank( cartComm );
  color();
  setGlobalGridValues( globalGridMin, globalGridMax );
  setLocalPartitionValues( globalGridMin );
}

void CartesianPartitioner::setNeighborsRank( const std::vector< int > & GEOS_UNUSED_PARAM( neighborsRank ) )
{}

void CartesianPartitioner::initializeCartesianCommunicator( MPI_Comm & cartComm )
{
  GEOS_ERROR_IF( m_partitionCounts.size() != 3, "Partition counts were not set before initializeCartesianCommunicator." );

  int size = MpiWrapper::commSize( MPI_COMM_GEOS );
  validatePartitionSize( size );

  int reorder = 0;
  MpiWrapper::cartCreate( MPI_COMM_GEOS, m_ndim, m_partitionCounts.data(),
                          m_periodic.data(), reorder, &cartComm );

  m_cartRank = MpiWrapper::commRank( cartComm );
  MpiWrapper::cartCoords( cartComm, m_cartRank, m_ndim, m_coords.data());
}

void CartesianPartitioner::validatePartitionSize( int size ) const
{
  string_view partitionsLogMessage =
    "The total number of processes = {} does not correspond to the total number of partitions = {}.\n"
    "The number of cells in an axis cannot be lower than the partition count of this axis\n";
  GEOS_ERROR_IF_NE_MSG( m_numPartitions, size, GEOS_FMT( partitionsLogMessage, size, m_numPartitions ));
}

void CartesianPartitioner::color()
{
  int color = 0;
  if( isOdd( m_coords[0] ))
    color += 1;
  if( isOdd( m_coords[1] ))
    color += 2;
  if( isOdd( m_coords[2] ))
    color += 4;

  // This color numbering may have gaps, so number of colors is just an upper bound, not a strict definition
  m_numColors = MpiWrapper::max( color ) + 1;
  m_color = color;
}

void CartesianPartitioner::processCommandLineOverrides( unsigned int xparCL, unsigned int yparCL, unsigned int zparCL )
{
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );

  // Case 1: User provided command-line overrides
  if( xparCL != 0 || yparCL != 0 || zparCL != 0 )
  {
    integer const xpar = (xparCL == 0) ? 1 : xparCL;
    integer const ypar = (yparCL == 0) ? 1 : yparCL;
    integer const zpar = (zparCL == 0) ? 1 : zparCL;

    integer const totalPartitions = xpar * ypar * zpar;

    // Validate that partition counts match MPI size
    GEOS_ERROR_IF( totalPartitions != mpiSize,
                   GEOS_FMT( "Partition count mismatch: -x {} -y {} -z {} = {} total partitions, "
                             "but running with {} MPI ranks. Product must equal MPI ranks.",
                             xpar, ypar, zpar, totalPartitions, mpiSize ));

    setPartitionCounts( xpar, ypar, zpar );
  }
  // Case 2: No command-line overrides
  else
  {
    GEOS_LOG_RANK_0 ( "CartesianPartitioner calculating default partition layout along z-axis." );
    // By default: a 1D decomposition along the z-axis.
    setPartitionCounts( 1, 1, mpiSize );
  }
}

void CartesianPartitioner::setPartitionCounts( unsigned int xPartitions, unsigned int yPartitions, unsigned int zPartitions )
{
  m_partitionCounts.clear();
  m_partitionCounts.emplace_back( xPartitions );
  m_partitionCounts.emplace_back( yPartitions );
  m_partitionCounts.emplace_back( zPartitions );

  m_numPartitions = xPartitions * yPartitions * zPartitions;
}

void CartesianPartitioner::addNeighbors( const unsigned int idim,
                                         int * ncoords,
                                         MPI_Comm cartComm )
{
  if( idim == m_ndim )
  {
    bool me = true;
    for( int i = 0; i < m_ndim; i++ )
    {
      if( ncoords[i] != this->m_coords( i ))
      {
        me = false;
        break;
      }
    }
    if( !me )
    {
      int const rank = MpiWrapper::cartRank( cartComm, ncoords );
      m_neighborsRank.push_back( rank );
    }
  }
  else
  {
    const int dim = this->m_partitionCounts( LvArray::integerConversion< localIndex >( idim ) );
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
        addNeighbors( idim + 1, ncoords, cartComm );
      }
    }
  }
}

void CartesianPartitioner::determineNeighborsRank( MPI_Comm cartComm )
{
  m_neighborsRank.clear();

  int ncoords[3];
  addNeighbors( 0, ncoords, cartComm );
}


void CartesianPartitioner::setGlobalGridValues( const real64 (& gridMin)[m_ndim], const real64 (& gridMax)[m_ndim] )
{
  LvArray::tensorOps::copy< 3 >( m_globalGridMin, gridMin );
  LvArray::tensorOps::copy< 3 >( m_globalGridMax, gridMax );
  LvArray::tensorOps::copy< 3 >( m_globalGridSize, gridMax );
  LvArray::tensorOps::subtract< 3 >( m_globalGridSize, gridMin );
}


void CartesianPartitioner::setLocalPartitionValues( const real64 (& gridMin)[m_ndim] )
{
  LvArray::tensorOps::copy< 3 >( m_localSize, m_globalGridSize );
  LvArray::tensorOps::copy< 3 >( m_localMin, gridMin );
  for( int i = 0; i < m_ndim; ++i )
  {
    m_localSize[i] /= m_partitionCounts( i );
    m_localMin[i] += m_coords( i ) * m_localSize[i];
    m_localMax[i] = gridMin[i] + (m_coords( i ) + 1) * m_localSize[i];
  }
}


bool CartesianPartitioner::isCoordInPartition( const real64 & coord, const int dir ) const
{
  bool rval = true;
  const int i = dir;
  if( m_periodic( i ))
  {
    if( m_partitionCounts( i ) != 1 )
    {
      real64 localCenter = MapValueToRange( coord, m_globalGridMin[ i ], m_globalGridMax[ i ] );
      rval = rval && localCenter >= m_localMin[ i ] && localCenter < m_localMax[ i ];
    }
  }
  else
  {
    rval = rval && (m_partitionCounts[ i ] == 1 || (coord >= m_localMin[ i ] && coord < m_localMax[ i ]));
  }

  return rval;
}

REGISTER_CATALOG_ENTRY( PartitionerBase, CartesianPartitioner, string const &, Group * const )

} // namespace geos
