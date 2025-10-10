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
 * @file HypreExport.cpp
 */

#include "HypreExport.hpp"

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "LvArray/src/sortedArrayManipulation.hpp"

#include <_hypre_parcsr_mv.h>
#include <_hypre_IJ_mv.h>
#include <algorithm>

namespace geos
{

HypreExport::HypreExport() = default;

HypreExport::HypreExport( HypreMatrix const & mat,
                          integer const targetRank )
  : m_targetRank( targetRank )
{
  // make a sub-communicator for scatter and ensure target rank is mapped to 0 in new comm
  int const rank = MpiWrapper::commRank( mat.comm() );
  int const color = ( mat.numLocalRows() > 0 ) ? 0 : MPI_UNDEFINED;
  int const key = ( rank == m_targetRank ) ? 0 : ( rank < m_targetRank ) ? rank + 1 : rank;
  m_subComm = MpiWrapper::commSplit( mat.comm(), color, key );
}

HypreExport::~HypreExport()
{
  if( m_subComm != MPI_COMM_NULL )
  {
    MpiWrapper::commFree( m_subComm );
  }
}

namespace
{

template< typename T >
void exportArray( HYPRE_MemoryLocation const location,
                  T const * const src,
                  arrayView1d< T > const & dst )
{
  // hypre does not maintain const-correctness of its APIs, hence the const_cast<>
  dst.move( hypre::getLvArrayMemorySpace( location ), true );
  hypre_TMemcpy( dst.data(), const_cast< T * >( src ), T, dst.size(), location, location );
}

template< typename T >
void exportArray( HYPRE_MemoryLocation const location,
                  arrayView1d< T const > const & src,
                  T * const dst )
{
  // hypre does not maintain const-correctness of its APIs, hence the const_cast<>
  src.move( hypre::getLvArrayMemorySpace( location ), false );
  hypre_TMemcpy( dst, const_cast< T * >( src.data() ), T, src.size(), location, location );
}

template< typename T, typename U >
void exportArray( HYPRE_MemoryLocation const location,
                  T const * const src,
                  arrayView1d< U > const & dst )
{
  if( location == HYPRE_MEMORY_HOST )
  {
    dst.move( hostMemorySpace, true );
    std::transform( src, src + dst.size(), dst.begin(),
                    []( T const v ) { return static_cast< U >( v ); } );
  }
  else // src is on device
  {
    forAll< hypre::execPolicy >( dst.size(), [dst, src] GEOS_HYPRE_DEVICE ( localIndex const i )
    {
      dst[i] = static_cast< U >( src[i] );
    } );
  }
}

template< typename T, typename U >
void exportArray( HYPRE_MemoryLocation const location,
                  arrayView1d< T const > const & src,
                  U * const dst )
{
  if( location == HYPRE_MEMORY_HOST )
  {
    src.move( hostMemorySpace, false );
    std::transform( src.begin(), src.end(), dst,
                    []( T const v ) { return static_cast< U >( v ); } );
  }
  else // dst is on device
  {
    forAll< hypre::execPolicy >( src.size(), [dst, src] GEOS_HYPRE_DEVICE ( localIndex const i )
    {
      dst[i] = static_cast< U >( src[i] );
    } );
  }
}

} // namespace

template< typename OFFSET_TYPE, typename COLUMN_TYPE >
void HypreExport::exportCRS( HypreMatrix const & mat,
                             arrayView1d< OFFSET_TYPE > const & rowOffsets,
                             arrayView1d< COLUMN_TYPE > const & colIndices,
                             arrayView1d< real64 > const & values ) const
{
  int const rank = MpiWrapper::commRank( mat.comm() );

  // import on target rank if needed, or extract diag+offd part in each rank
  hypre_CSRMatrix * const localMatrix = m_targetRank < 0
                                      ? hypre_MergeDiagAndOffd( mat.unwrapped() )
                                      : hypre_ParCSRMatrixToCSRMatrixAll( mat.unwrapped() );
  GEOS_ERROR_IF( rank == m_targetRank && !localMatrix, "HypreExport: matrix is empty on target rank" );

  if( m_targetRank < 0 || m_targetRank == rank )
  {
    HYPRE_Int const numRow = hypre_CSRMatrixNumRows( localMatrix );
    HYPRE_Int const numNz  = hypre_CSRMatrixNumNonzeros( localMatrix );

    GEOS_LAI_ASSERT_EQ( rowOffsets.size(), numRow + 1 );
    GEOS_LAI_ASSERT_EQ( colIndices.size(), numNz );
    GEOS_LAI_ASSERT_EQ( values.size(), numNz );

    HYPRE_MemoryLocation const location = hypre_CSRMatrixMemoryLocation( localMatrix );

    exportArray( location, hypre_CSRMatrixI( localMatrix ), rowOffsets );
    exportArray( location, hypre_CSRMatrixData( localMatrix ), values );

    // We have to handle two cases differently because hypre uses two different struct members
    // (j/big_j) to store the column indices depending on how we obtained the local matrix.
    if( m_targetRank < 0 )
    {
      exportArray( location, hypre_CSRMatrixBigJ( localMatrix ), colIndices );
    }
    else
    {
      exportArray( location, hypre_CSRMatrixJ( localMatrix ), colIndices );
    }

    // Sort the values by column index after copying (some solvers expect this)
    forAll< hypre::execPolicy >( numRow, [rowOffsets, colIndices, values] GEOS_HYPRE_DEVICE ( HYPRE_Int const i )
    {
      LvArray::sortedArrayManipulation::dualSort( colIndices.data() + rowOffsets[i],
                                                  colIndices.data() + rowOffsets[i + 1],
                                                  values.data() + rowOffsets[i] );
    } );
  }

  GEOS_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( localMatrix ) );
}

void HypreExport::exportVector( HypreVector const & vec,
                                arrayView1d< real64 > const & values ) const
{
  int const rank = MpiWrapper::commRank( vec.comm() );

  // Gather vector on target rank, or just get the local part
  hypre_Vector * const localVector = m_targetRank < 0
                                   ? hypre_ParVectorLocalVector( vec.unwrapped() )
                                   : (hypre_Vector *)hypre::parVectorToVector( vec.unwrapped(), m_targetRank );
  GEOS_ERROR_IF( rank == m_targetRank && !localVector, "HypreExport: vector is empty on target rank" );

  if( m_targetRank < 0 || m_targetRank == rank )
  {
    GEOS_LAI_ASSERT_EQ( values.size(), hypre_VectorSize( localVector ) );
    exportArray( hypre_VectorMemoryLocation( localVector ),
                 hypre_VectorData( localVector ),
                 values );
  }

  if( m_targetRank >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( localVector ) );
  }
}

void HypreExport::importVector( arrayView1d< real64 const > const & values,
                                HypreVector & vec ) const
{
  hypre_Vector * const localVector = hypre_ParVectorLocalVector( vec.unwrapped() );

  // Single-rank import: scatter the global vector from target rank to all ranks with rows
  if( m_targetRank >= 0 && m_subComm != MPI_COMM_NULL )
  {
    int const subRank = MpiWrapper::commRank( m_subComm );
    int const subSize = MpiWrapper::commSize( m_subComm );
    int const parentRank = MpiWrapper::commRank( vec.comm() );

    // Root of the sub-communicator must be the target rank (the constructor maps it to 0)
    if( subRank == 0 )
    {
      GEOS_ERROR_IF( parentRank != m_targetRank,
                     "HypreExport::importVector: target rank has no rows and is not in the sub-communicator" );
    }

    // Build counts and displacements from actual local sizes
    int const myLocal = LvArray::integerConversion< int >( vec.localSize() );
    array1d< int > counts( subSize );
    MPI_CHECK_ERROR( MpiWrapper::allgather( &myLocal, 1, counts.data(), 1, m_subComm ) );

    array1d< int > displs( subSize );
    displs[0] = 0;
    for( int i = 1; i < subSize; ++i )
    {
      displs[i] = displs[i-1] + counts[i-1];
    }

    // Root validates buffer size and prepares send pointer
    real64 const * sendBuf = nullptr;
    if( subRank == 0 )
    {
      int const total = displs[subSize-1] + counts[subSize-1];
      GEOS_ERROR_IF_NE( vec.globalSize(), total );
      GEOS_ERROR_IF_NE( values.size(), vec.globalSize() );
      values.move( hostMemorySpace, false );
      sendBuf = values.data();
    }

    // Receive local chunk into host buffer and then copy to hypre local vector (host/device)
    array1d< real64 > recvBuf( myLocal );

    MPI_CHECK_ERROR( MpiWrapper::scatterv( sendBuf,
                                           counts.data(),
                                           displs.data(),
                                           recvBuf.data(),
                                           myLocal,
                                           /*root*/ 0,
                                           m_subComm ) );

    hypre_TMemcpy( hypre_VectorData( localVector ),
                   recvBuf.data(),
                   HYPRE_Real,
                   vec.localSize(),
                   hypre_VectorMemoryLocation( localVector ),
                   HYPRE_MEMORY_HOST );
  }
  else
  {
    exportArray( hypre_VectorMemoryLocation( localVector ),
                 values,
                 hypre_VectorData( localVector ) );
  }

  vec.touch();
}

/**
 * Explicit template instantiation macro for HypreExport::exportCRS.
 * We need to explicitly instantiate this function template because:
 * - we want CRS consumers to specify their own destination buffer types;
 * - we're "hiding" Hypre headers from being included by consumer code.
 */
#define INST_HYPRE_EXPORT_CRS( OFFSET_TYPE, COLUMN_TYPE ) \
  template void \
  HypreExport::exportCRS< OFFSET_TYPE, COLUMN_TYPE >( HypreMatrix const &, \
                                                      arrayView1d< OFFSET_TYPE > const &, \
                                                      arrayView1d< COLUMN_TYPE > const &, \
                                                      arrayView1d< real64 > const & ) const

// Add other instantiations as needed (only use built-in types)
INST_HYPRE_EXPORT_CRS( int, int );
INST_HYPRE_EXPORT_CRS( long, long );
INST_HYPRE_EXPORT_CRS( long long, long long );

#undef INST_HYPRE_EXPORT_CRS

} // namespace geos
