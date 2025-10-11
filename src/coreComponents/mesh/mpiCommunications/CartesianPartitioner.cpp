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
real64 Mod( real64 const num, real64 const denom )
{
  if( fabs( denom ) < fabs( num ) * 1.0e-14 )
  {
    return num;
  }

  return num - denom * std::floor( num / denom );
}

// MapValueToRange
// returns a periodic value in the range [min, max)
real64 MapValueToRange( real64 const value, real64 const min, real64 const max )
{
  return Mod( value - min, max - min ) + min;
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
                                            Group * const parent )
  : GeometricPartitioner( name, parent ),
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
  registerWrapper( viewKeyStruct::partitionCountsString(), &m_partitionCounts )
    .setInputFlag( InputFlags::OPTIONAL )
    .setDescription( "Number of partitions in each direction [nx, ny, nz]" );
}

void CartesianPartitioner::initializeDomain( R1Tensor const & globalGridMin,
                                             R1Tensor const & globalGridMax )
{
  ScopedMpiComm cartComm;
  initializeCartesianCommunicator( cartComm );

  // Compute neighbors from Cartesian topology
  computeNeighborsFromTopology( cartComm );

  // Compute checkerboard coloring
  color();

  // Set global and local domain values
  setGlobalGridValues( globalGridMin, globalGridMax );
  setLocalPartitionValues( globalGridMin );
}

void CartesianPartitioner::initializeCartesianCommunicator( MPI_Comm & cartComm )
{
  GEOS_ERROR_IF( m_partitionCounts.size() != 3,
                 "Partition counts were not set before initializeCartesianCommunicator." );

  int const size = MpiWrapper::commSize( MPI_COMM_GEOS );
  validatePartitionSize( size );


  int const reorder = 0;
  MpiWrapper::cartCreate( MPI_COMM_GEOS, m_ndim, m_partitionCounts.data(),
                          m_periodic.data(), reorder, &cartComm );

  m_cartRank = MpiWrapper::commRank( cartComm );
  MpiWrapper::cartCoords( cartComm, m_cartRank, m_ndim, m_coords.data() );
}

void CartesianPartitioner::validatePartitionSize( int const size ) const
{
  GEOS_ERROR_IF( m_numPartitions != size,
                 GEOS_FMT( "The total number of processes = {} does not correspond to the total number of partitions = {}.\n"
                           "The number of cells in an axis cannot be lower than the partition count of this axis",
                           size, m_numPartitions ) );
}

void CartesianPartitioner::color()
{
  // Specialized checkerboard coloring for Cartesian grid
  // Much more efficient than generic graph coloring!

  int color = 0;
  if( isOdd( m_coords[0] ) )
    color += 1;
  if( isOdd( m_coords[1] ) )
    color += 2;
  if( isOdd( m_coords[2] ) )
    color += 4;

  // This color numbering may have gaps, so number of colors is just an upper bound
  m_numColors = MpiWrapper::max( color, MPI_COMM_GEOS ) + 1;
  m_color = color;

  GEOS_LOG_RANK_0( GEOS_FMT( "CartesianPartitioner: Checkerboard coloring complete. "
                             "Using {} colors (upper bound)", m_numColors ) );
}

void CartesianPartitioner::processCommandLineOverrides( unsigned int const xparCL,
                                                        unsigned int const yparCL,
                                                        unsigned int const zparCL )
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
    GEOS_THROW_IF( totalPartitions != mpiSize,
                   GEOS_FMT( "Partition count mismatch: -x {} -y {} -z {} = {} total partitions, "
                             "but running with {} MPI ranks. Product must equal MPI ranks.",
                             xpar, ypar, zpar, totalPartitions, mpiSize ), InputError );

    setPartitionCounts( xpar, ypar, zpar );
  }
  // Case 2: No command-line overrides - check if XML values are set
  else if( m_partitionCounts.empty() )
  {
    GEOS_LOG_RANK_0( "CartesianPartitioner: No partition counts specified. "
                     "Using default 1D decomposition along z-axis." );
    setPartitionCounts( 1, 1, mpiSize );
  }
  // Case 3: XML values already set, validate them
  else
  {
    integer const totalPartitions = m_partitionCounts[0] * m_partitionCounts[1] * m_partitionCounts[2];
    GEOS_THROW_IF( totalPartitions != mpiSize,
                   GEOS_FMT( "Partition count mismatch: XML specifies {} x {} x {} = {} total partitions, "
                             "but running with {} MPI ranks. Product must equal MPI ranks.",
                             m_partitionCounts[0], m_partitionCounts[1], m_partitionCounts[2],
                             totalPartitions, mpiSize ), InputError );
  }
}

void CartesianPartitioner::setPartitionCounts( unsigned int const xPartitions,
                                               unsigned int const yPartitions,
                                               unsigned int const zPartitions )
{
  m_partitionCounts.resize( 3 );
  m_partitionCounts[0] = xPartitions;
  m_partitionCounts[1] = yPartitions;
  m_partitionCounts[2] = zPartitions;

  m_numPartitions = xPartitions * yPartitions * zPartitions;
}

void CartesianPartitioner::addNeighbors( unsigned int const idim,
                                         int * ncoords,
                                         MPI_Comm const cartComm )
{
  if( idim == m_ndim )
  {
    bool me = true;
    for( int i = 0; i < m_ndim; i++ )
    {
      if( ncoords[i] != this->m_coords( i ) )
      {
        me = false;
        break;
      }
    }
    if( !me )
    {
      int const rank = MpiWrapper::cartRank( cartComm, ncoords );

      if( std::find( m_neighborsRank.begin(), m_neighborsRank.end(), rank )
          == m_neighborsRank.end() )
      {
        m_neighborsRank.push_back( rank );
      }
    }
  }
  else
  {
    int const dim = this->m_partitionCounts( LvArray::integerConversion< localIndex >( idim ) );
    bool const periodic = this->m_periodic( LvArray::integerConversion< localIndex >( idim ) );
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

void CartesianPartitioner::computeNeighborsFromTopology( MPI_Comm & cartComm )
{
  m_neighborsRank.clear();

  int ncoords[3];

  // Recursively find all neighbors in 3×3×3 stencil
  addNeighbors( 0, ncoords, cartComm );
}

void CartesianPartitioner::setGlobalGridValues( R1Tensor const & gridMin,
                                                R1Tensor const & gridMax )
{
  LvArray::tensorOps::copy< 3 >( m_globalGridMin, gridMin );
  LvArray::tensorOps::copy< 3 >( m_globalGridMax, gridMax );
  LvArray::tensorOps::copy< 3 >( m_globalGridSize, gridMax );
  LvArray::tensorOps::subtract< 3 >( m_globalGridSize, gridMin );
}

void CartesianPartitioner::setLocalPartitionValues( R1Tensor const & gridMin )
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

bool CartesianPartitioner::isCoordInPartition( R1Tensor const & coords ) const
{
  return isCoordInPartition( coords[0], 0 ) &&
         isCoordInPartition( coords[1], 1 ) &&
         isCoordInPartition( coords[2], 2 );
}

bool CartesianPartitioner::isCoordInPartition( real64 const & coord, int const dir ) const
{
  bool rval = true;
  if( m_periodic( dir ) )
  {
    if( m_partitionCounts( dir ) != 1 )
    {
      real64 const localCenter = MapValueToRange( coord, m_globalGridMin[dir], m_globalGridMax[dir] );
      rval = rval && localCenter >= m_localMin[dir] && localCenter < m_localMax[dir];
    }
  }
  else
  {
    rval = rval && (m_partitionCounts[dir] == 1 || (coord >= m_localMin[dir] && coord < m_localMax[dir]));
  }

  return rval;
}

// Register in DomainPartitioner catalog
REGISTER_CATALOG_ENTRY( DomainPartitioner, CartesianPartitioner, string const &, Group * const )

} // namespace geos
