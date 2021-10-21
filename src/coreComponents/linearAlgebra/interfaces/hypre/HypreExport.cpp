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
#include "LvArray/src/sortedArrayManipulation.hpp"

#include <_hypre_parcsr_mv.h>
#include <_hypre_IJ_mv.h>

namespace geosx
{

HypreExport::HypreExport() = default;

HypreExport::HypreExport( HypreMatrix const & mat,
                          integer targetRank )
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

template< typename OFFSET_TYPE, typename COLUMN_TYPE >
void HypreExport::exportCRS( HypreMatrix const & mat,
                             OFFSET_TYPE * const rowOffsets,
                             COLUMN_TYPE * const colIndices,
                             real64 * const values ) const
{
  int const rank = MpiWrapper::commRank( mat.getComm() );

  // import on target rank if needed, or extract diag+offdiag part in each rank
  hypre_CSRMatrix * localMatrix;
  if( m_targetRank < 0 )
  {
    localMatrix = hypre_MergeDiagAndOffd( mat.unwrapped() );
  }
  else
  {
    GEOSX_ERROR_IF( rank == m_targetRank && !( rowOffsets && colIndices && values ),
                    "HypreExport: must pass non-null pointers on the target rank" );
    localMatrix = hypre_ParCSRMatrixToCSRMatrixAll( mat.unwrapped() );
    GEOSX_ERROR_IF( rank == m_targetRank && !localMatrix, "HypreExport: matrix is empty on target rank" );
  }

  // export the raw CRS data
  if( m_targetRank < 0 || m_targetRank == rank )
  {
    HYPRE_Int const numRows = hypre_CSRMatrixNumRows( localMatrix );
    HYPRE_Int const numNz = hypre_CSRMatrixNumNonzeros( localMatrix );
    HYPRE_Int const * const ia = hypre_CSRMatrixI( localMatrix );
    HYPRE_Real const * const va = hypre_CSRMatrixData( localMatrix );

    std::transform( ia, ia + numRows + 1, rowOffsets, LvArray::integerConversion< OFFSET_TYPE, HYPRE_Int > );
    std::copy( va, va + numNz, values );

    // We have to handle two cases differently because hypre uses two different struct members
    // (j/big_j) to store the column indices depending on how we obtained the local matrix.
    if( m_targetRank < 0 )
    {
      HYPRE_BigInt const * const ja = hypre_CSRMatrixBigJ( localMatrix );
      std::transform( ja, ja + numNz, colIndices, LvArray::integerConversion< COLUMN_TYPE, HYPRE_BigInt > );
    }
    else
    {
      HYPRE_Int const * const ja = hypre_CSRMatrixJ( localMatrix );
      std::transform( ja, ja + numNz, colIndices, LvArray::integerConversion< COLUMN_TYPE, HYPRE_Int > );
    }

    // Sort the values by column index after copying (some solvers expect this)
    for( localIndex i = 0; i < numRows; ++i )
    {
      using LvArray::sortedArrayManipulation::dualSort;
      dualSort( colIndices + rowOffsets[i], colIndices + rowOffsets[i + 1], values + rowOffsets[i] );
    }
  }

  GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( localMatrix ) );
}

void HypreExport::exportVector( HypreVector const & vec,
                                real64 * values ) const
{
  if( m_targetRank >= 0 )
  {
    int const rank = MpiWrapper::commRank( vec.getComm() );
    GEOSX_ERROR_IF( rank == m_targetRank && !values, "HypreExport: must pass non-null pointers on the target rank" );

    hypre_Vector * const localVector = hypre_ParVectorToVectorAll( vec.unwrapped() );
    if( rank == m_targetRank )
    {
      HYPRE_Real const * const data = hypre_VectorData( localVector );
      std::copy( data, data + vec.globalSize(), values );
    }
    GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( localVector ) );
  }
  else
  {
    real64 const * const data = vec.extractLocalVector();
    std::copy( data, data + vec.localSize(), values );
  }
}

void HypreExport::importVector( real64 const * values,
                                HypreVector & vec ) const
{
  if( m_targetRank >= 0 )
  {
    int const rank = MpiWrapper::commRank( vec.getComm() );
    GEOSX_ERROR_IF( rank == m_targetRank && !values, "HypreExport: must pass non-null pointers on the target rank" );

    if( m_subComm != MPI_COMM_NULL )
    {
      HYPRE_Vector localVector = nullptr;
      if( rank == m_targetRank )
      {
        // HACK: create a hypre vector that points to local data; we have to use const_cast,
        //       but this is ok because we don't modify the values, only scatter the vector.
        localVector = HYPRE_VectorCreate( LvArray::integerConversion< HYPRE_Int >( vec.globalSize() ) );
        hypre_VectorOwnsData( ( hypre_Vector * ) localVector ) = false;
        hypre_VectorData( ( hypre_Vector * ) localVector ) = const_cast< real64 * >( values );
        hypre_SeqVectorInitialize_v2( ( hypre_Vector * ) localVector, HYPRE_MEMORY_HOST );
      }

      // output vector partitioning array
      HYPRE_BigInt * const partitioning = hypre_ParVectorPartitioning( vec.unwrapped() );

      // scatter the data
      HYPRE_ParVector parVector;
      GEOSX_LAI_CHECK_ERROR( HYPRE_VectorToParVector( m_subComm,
                                                      localVector,
                                                      partitioning,
                                                      &parVector ) );


      // copy parVector local part over to the output vector
      HYPRE_Real * const parVectorData = hypre_VectorData( hypre_ParVectorLocalVector( parVector ) );
      std::copy( parVectorData, parVectorData + vec.localSize(), vec.extractLocalVector() );

      GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorDestroy( parVector ) );
      GEOSX_LAI_CHECK_ERROR( HYPRE_VectorDestroy( localVector ) );
    }
  }
  else
  {
    real64 * const data = vec.extractLocalVector();
    std::copy( values, values + vec.localSize(), data );
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
                                                      OFFSET_TYPE * const, \
                                                      COLUMN_TYPE * const, \
                                                      real64 * const ) const

// Add other instantiations as needed (only use built-in types)
INST_HYPRE_EXPORT_CRS( int, int );
INST_HYPRE_EXPORT_CRS( long, long );
INST_HYPRE_EXPORT_CRS( long long, long long );

#undef INST_HYPRE_EXPORT_CRS

} // namespace geosx
