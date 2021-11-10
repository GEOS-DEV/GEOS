/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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

namespace geosx
{

HypreExport::HypreExport() = default;

HypreExport::HypreExport( HypreMatrix const & mat,
                          integer const targetRank )
  : m_targetRank( targetRank )
{
  // make a sub-communicator for scatter and ensure target rank is mapped to 0 in new comm
  int const rank = MpiWrapper::commRank( mat.getComm() );
  int const color = ( mat.numLocalRows() > 0 ) ? 0 : MPI_UNDEFINED;
  int const key = ( rank == m_targetRank ) ? 0 : ( rank < m_targetRank ) ? rank + 1 : rank;
  m_subComm = MpiWrapper::commSplit( mat.getComm(), color, key );
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

template< typename T, typename U >
std::enable_if_t< std::is_same< T, U >::value >
copyOrTransform( T const * const first,
                 T const * const last,
                 U * const out )
{
  std::copy( first, last, out );
}

template< typename T, typename U >
std::enable_if_t< !std::is_same< T, U >::value >
copyOrTransform( T const * const first,
                 T const * const last,
                 U * const out )
{
  std::transform( first, last, out, []( T const v ) { return static_cast< U >( v ); } );
}

template< typename HYPRE_TYPE, typename GEOSX_TYPE >
void exportArray( HYPRE_MemoryLocation const memorySpace,
                  HYPRE_TYPE const * const hypreArray,
                  arrayView1d< GEOSX_TYPE > const & geosxArray )
{
  if( memorySpace == HYPRE_MEMORY_HOST )
  {
    geosxArray.move( LvArray::MemorySpace::host, true );
    copyOrTransform( hypreArray, hypreArray + geosxArray.size(), geosxArray.data() );
  }
  else // hypreArray is on device
  {
    forAll< hypre::execPolicy >( geosxArray.size(),
                                 [geosxArray, hypreArray] GEOSX_HYPRE_DEVICE ( localIndex const i )
    {
      geosxArray[i] = static_cast< GEOSX_TYPE >( hypreArray[i] );
    } );
  }
}

template< typename HYPRE_TYPE, typename GEOSX_TYPE >
void importArray( HYPRE_MemoryLocation const memorySpace,
                  arrayView1d< GEOSX_TYPE const > const & geosxArray,
                  HYPRE_TYPE * const hypreArray )
{
  if( memorySpace == HYPRE_MEMORY_HOST )
  {
    geosxArray.move( LvArray::MemorySpace::host, false );
    copyOrTransform( geosxArray.data(), geosxArray.data() + geosxArray.size(), hypreArray );
  }
  else // hypreArray is on device
  {
    forAll< hypre::execPolicy >( geosxArray.size(),
                                 [geosxArray, hypreArray] GEOSX_HYPRE_DEVICE ( localIndex const i )
    {
      hypreArray[i] = static_cast< HYPRE_TYPE >( geosxArray[i] );
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
  int const rank = MpiWrapper::commRank( mat.getComm() );

  // import on target rank if needed, or extract diag+offd part in each rank
  hypre_CSRMatrix * const localMatrix = m_targetRank < 0
                                      ? hypre_MergeDiagAndOffd( mat.unwrapped() )
                                      : hypre_ParCSRMatrixToCSRMatrixAll( mat.unwrapped() );
  GEOSX_ERROR_IF( rank == m_targetRank && !localMatrix, "HypreExport: matrix is empty on target rank" );

  if( m_targetRank < 0 || m_targetRank == rank )
  {
    HYPRE_Int const numRow = hypre_CSRMatrixNumRows( localMatrix );
    HYPRE_Int const numNz  = hypre_CSRMatrixNumNonzeros( localMatrix );

    GEOSX_LAI_ASSERT_EQ( rowOffsets.size(), numRow + 1 );
    GEOSX_LAI_ASSERT_EQ( colIndices.size(), numNz );
    GEOSX_LAI_ASSERT_EQ( values.size(), numNz );

    HYPRE_MemoryLocation const memorySpace = hypre_CSRMatrixMemoryLocation( localMatrix );

    exportArray( memorySpace, hypre_CSRMatrixI( localMatrix ), rowOffsets );
    exportArray( memorySpace, hypre_CSRMatrixData( localMatrix ), values );

    // We have to handle two cases differently because hypre uses two different struct members
    // (j/big_j) to store the column indices depending on how we obtained the local matrix.
    if( m_targetRank < 0 )
    {
      exportArray( memorySpace, hypre_CSRMatrixBigJ( localMatrix ), colIndices );
    }
    else
    {
      exportArray( memorySpace, hypre_CSRMatrixJ( localMatrix ), colIndices );
    }

    // Sort the values by column index after copying (some solvers expect this)
    forAll< hypre::execPolicy >( numRow, [rowOffsets, colIndices, values] GEOSX_HYPRE_DEVICE ( HYPRE_Int const i )
    {
      using LvArray::sortedArrayManipulation::dualSort;
      dualSort( colIndices.data() + rowOffsets[i],
                colIndices.data() + rowOffsets[i + 1],
                values.data() + rowOffsets[i] );
    } );
  }

  GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( localMatrix ) );
}

void HypreExport::exportVector( HypreVector const & vec,
                                arrayView1d< real64 > const & values ) const
{
  int const rank = MpiWrapper::commRank( vec.getComm() );

  // Gather vector on target rank, or just get the local part
  hypre_Vector * const localVector = m_targetRank < 0
                                   ? hypre_ParVectorLocalVector( vec.unwrapped() )
                                   : hypre_ParVectorToVectorAll( vec.unwrapped() );
  GEOSX_ERROR_IF( rank == m_targetRank && !localVector, "HypreExport: vector is empty on target rank" );

  if( m_targetRank < 0 || m_targetRank == rank )
  {
    GEOSX_LAI_ASSERT_EQ( values.size(), hypre_VectorSize( localVector ) );
    exportArray( hypre_VectorMemoryLocation( localVector ),
                 hypre_VectorData( localVector ),
                 values );
  }

  if( m_targetRank >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( localVector ) );
  }
}

void HypreExport::importVector( arrayView1d< real64 const > const & values,
                                HypreVector & vec ) const
{
  if( m_subComm != MPI_COMM_NULL )
  {
    hypre_Vector * localVector{};
    if( MpiWrapper::commRank( vec.getComm() ) == m_targetRank )
    {
      GEOSX_LAI_ASSERT_EQ( values.size(), vec.globalSize() );
      values.move( LvArray::MemorySpace::host, false );

      // HACK: create a hypre vector that points to local data; we have to use const_cast,
      //       but this is ok because we don't modify the values, only scatter the vector.
      localVector = hypre_SeqVectorCreate( LvArray::integerConversion< HYPRE_Int >( values.size() ) );
      hypre_VectorOwnsData( localVector ) = false;
      hypre_VectorData( localVector ) = const_cast< real64 * >( values.data() );
      hypre_SeqVectorInitialize_v2( localVector, HYPRE_MEMORY_HOST );
    }

    // scatter the data
    hypre_ParVector * const parVector = hypre_VectorToParVector( m_subComm,
                                                                 localVector,
                                                                 hypre_ParVectorPartitioning( vec.unwrapped() ) );
    // copy local part of the data over to the output vector
    HYPRE_Real const * const parVectorData = hypre_VectorData( hypre_ParVectorLocalVector( parVector ) );
    if( hypre_ParVectorMemoryLocation( vec.unwrapped() ) == HYPRE_MEMORY_HOST )
    {
      std::copy( parVectorData, parVectorData + vec.localSize(), vec.extractLocalVector() );
    }
    else
    {
#ifdef GEOSX_USE_HYPRE_CUDA
      cudaMemcpy( vec.extractLocalVector(),
                  parVectorData,
                  vec.localSize() * sizeof( HYPRE_Real ),
                  cudaMemcpyHostToDevice );
#else
      GEOSX_ERROR( "HypreExport: invalid memory space of target vector" );
#endif
    }

    GEOSX_LAI_CHECK_ERROR( hypre_ParVectorDestroy( parVector ) );
    GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( localVector ) );
  }
  else
  {
    hypre_Vector * const localVector = hypre_ParVectorLocalVector( vec.unwrapped() );
    importArray( hypre_VectorMemoryLocation( localVector ),
                 values,
                 hypre_VectorData( localVector ) );
  }
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

} // namespace geosx
