
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


//#include <mpi.h>
//#include "LvArray/genericTensorOps.hpp"


#include "CartesianPartitioner.hpp"




namespace geos {


void CartesianPartitioner::partition() {
}

std::vector<int> CartesianPartitioner::getNeighbors() const {

#if 0
    // Return neighbor partitions
    int rank;
    MPI_Comm_rank(m_cartComm, &rank);

    int coords;
    MPI_Cart_coords(m_cartComm, rank, 2, coords);

    std::vector<int> neighbors;

    for (int i = 0; i < 2; ++i) {
        int neighborCoords = {coords, coords};
        neighborCoords[i] = (coords[i] + 1) % dims[i];
        int neighborRank;
        MPI_Cart_rank(m_cartComm, neighborCoords, &neighborRank);
        neighbors.push_back(neighborRank);

        neighborCoords[i] = (coords[i] - 1 + dims[i]) % dims[i];
        MPI_Cart_rank(m_cartComm, neighborCoords, &neighborRank);
        neighbors.push_back(neighborRank);
        
    }

    return neighbors;

  
#endif    
    std::vector<int> u;
    return u;
}

int CartesianPartitioner::getColor() const {

return -1;
}



#if 0
void CartesianPartitioner::addNeighbors( const unsigned int idim,
                                         MPI_Comm & cartcomm,
                                         int * ncoords )
{


#if 0

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

     std::cout << "Adding to: "<< MpiWrapper::commRank(MPI_COMM_GEOS) << " neighbor: "<< rank <<std::endl;
     //m_neighbors.push_back( NeighborCommunicator( rank ) );
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

#endif  
}

#endif



} // namespace geos












#if 0
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
  m_Periodic( nsdof ),
  m_coords( nsdof ),
  m_min{ 0.0 },
  m_max{ 0.0 },
  m_blockSize{ 1.0 },
  m_gridSize{ 0.0 },
  m_gridMin{ 0.0 },
  m_gridMax{ 0.0 },
  m_Partitions()
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

std::cout << "Adding to: "<< MpiWrapper::commRank(MPI_COMM_GEOS) << " neighbor: "<< rank <<std::endl;
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

void SpatialPartition::updateSizes( arrayView1d< real64 > const domainL,
                                    real64 const dt )
{
  for( int i=0; i<3; i++ )
  {
    real64 ratio = 1.0 + domainL[i] * dt;
    m_min[i] *= ratio;
    m_max[i] *= ratio;
    //m_PartitionLocations[i] *= ratio; ?
    m_blockSize[i] *= ratio;
    m_gridSize[i] *= ratio;
    m_gridMin[i] *= ratio;
    m_gridMax[i] *= ratio;
    m_contactGhostMin[i] *= ratio;
    m_contactGhostMax[i] *= ratio;
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
      int check = 1;
      for( int i = 0; i < nsdof; i++ )
      {
        check *= this->m_Partitions( i );
      }
      GEOS_ERROR_IF_NE( check, m_size );
    }

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOS, nsdof, m_Partitions.data(), m_Periodic.data(), reorder, &cartcomm );
    }
    m_rank = MpiWrapper::commRank( cartcomm );
    MpiWrapper::cartCoords( cartcomm, m_rank, nsdof, m_coords.data());

    //add neighbors
    // {
    //   int ncoords[nsdof];
    //   m_neighbors.clear();
    //   addNeighbors( 0, cartcomm, ncoords );
    // }

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
      GEOS_ERROR( "SpatialPartition::setSizes(): number of partition locations does not equal number of partitions - 1\n" );
    }
  }
}

bool SpatialPartition::isCoordInPartition( const real64 & coord, const int dir ) const
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

bool SpatialPartition::isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                                      const real64 & boundaryRadius ) const
// test a point relative to a boundary box. If non-zero buffer specified, expand the box.
{
  for( int i = 0; i < nsdof; i++ )
  {
    // Is particle already in bounds of partition?
    if( !(m_Partitions( i )==1 || ( elemCenter[i] >= (m_min[i] - boundaryRadius) && elemCenter[i] <= (m_max[i] + boundaryRadius) ) ) )
    {
      // Particle not in bounds, check if direction has a periodic boundary
      if( m_Periodic( i ) && (m_coords[i] == 0 || m_coords[i] == m_Partitions[i] - 1) )
      {
        // Partition minimum boundary is periodic
        if( m_coords[i] == 0 && ( (elemCenter[i] - m_gridSize[i]) < (m_min[i] - boundaryRadius) ) )
        {
          return false;
        }
        // Partition maximum boundary is periodic
        if( m_coords[i] == m_Partitions[i] - 1 && ( (elemCenter[i] + m_gridSize[i]) > (m_max[i] + boundaryRadius) ) )
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






#endif