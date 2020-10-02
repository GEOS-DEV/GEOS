/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreSuiteSparse.cpp
 */

#include <numeric>

#include "HypreSuiteSparse.hpp"
#include "common/Stopwatch.hpp"

#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"
#include "seq_mv.h"

namespace geosx
{

namespace
{
void sortIntReal64( Int * arrayInt, real64 * arrayReal64, localIndex const size )
{
  array1d< localIndex > indices( size );
  std::iota( indices.begin(), indices.end(), 0 );
  std::sort( indices.begin(), indices.end(), [&]( localIndex const & a, localIndex const & b ) { return arrayInt[a] < arrayInt[b]; } );
  array1d< localIndex > permutation( size );
  for( localIndex i = 0; i < size; ++i )
  {
    permutation[indices[i]] = i;
  }
  for( localIndex s = 1, d; s < size; ++s )
  {
    for( d = permutation[s]; d < s; d = permutation[d] )
      ;
    if( d == s )
    {
      while( d = permutation[d], d != s )
      {
        std::swap( arrayInt[s], arrayInt[d] );
        std::swap( arrayReal64[s], arrayReal64[d] );
      }
    }
  }
}
}

void ConvertHypreToSuiteSparseMatrix( HypreMatrix const & matrix,
                                      SuiteSparseData & SSData )
{
  // Copy distributed parcsr matrix in a local CSR matrix on every process
  // with at least one row
  // Warning: works for a parcsr matrix that is smaller than 2^31-1
  hypre_CSRMatrix * CSRmatrix = hypre_ParCSRMatrixToCSRMatrixAll( matrix.unwrapped() );

  // Identify the smallest process where CSRmatrix exists
  int rank = MpiWrapper::Comm_rank( matrix.getComm() );
  if( CSRmatrix == 0 )
  {
    rank = MpiWrapper::Comm_size( matrix.getComm() );
  }
  SSData.workingRank = MpiWrapper::Min( rank, matrix.getComm() );

  rank = MpiWrapper::Comm_rank( matrix.getComm() );
  if( rank == SSData.workingRank )
  {
    SSData.numRows = toSuiteSparse_Int( matrix.numGlobalRows() );
    SSData.numCols = toSuiteSparse_Int( matrix.numGlobalCols() );
    SSData.nonZeros = toSuiteSparse_Int( hypre_CSRMatrixNumNonzeros( CSRmatrix ) );
    HYPRE_Int const * const hypreI = hypre_CSRMatrixI( CSRmatrix );
    SSData.rowPtr = new Int[SSData.numRows+1];
    for( localIndex i = 0; i <= SSData.numRows; ++i )
    {
      SSData.rowPtr[i] = LvArray::integerConversion< Int >( hypreI[i] );
    }
    SSData.colIndices = new Int[SSData.nonZeros];
    SSData.data = new real64[SSData.nonZeros];
    HYPRE_Int const * const hypreJ = hypre_CSRMatrixJ( CSRmatrix );
    for( localIndex i = 0; i < SSData.nonZeros; ++i )
    {
      SSData.colIndices[i] = LvArray::integerConversion< Int >( hypreJ[i] );
    }
    std::copy( hypre_CSRMatrixData( CSRmatrix ),
               hypre_CSRMatrixData( CSRmatrix ) + SSData.nonZeros,
               SSData.data );
    for( localIndex i = 0; i < SSData.numRows; ++i )
    {
      localIndex rowLength = LvArray::integerConversion< localIndex >( SSData.rowPtr[i+1] - SSData.rowPtr[i] );
      sortIntReal64( &( SSData.colIndices[SSData.rowPtr[i]] ), &( SSData.data[SSData.rowPtr[i]] ), rowLength );
    }
  }

  // Destroy CSRmatrix
  if( CSRmatrix )
  {
    GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( CSRmatrix ) );
  }

  // Save communicator
  SSData.comm = matrix.getComm();
}

int SuiteSparseSolve( SuiteSparseData & SSData,
                      HypreVector const & b,
                      HypreVector & x,
                      real64 & time )
{
  int status = 0;

  hypre_Vector * b_vector = hypre_ParVectorToVectorAll( b.unwrapped() );
  HYPRE_Vector sol_Vector = nullptr;

  int const rank = MpiWrapper::Comm_rank( SSData.comm );
  if( rank == SSData.workingRank )
  {
    // Create local vector to store the solution
    hypre_Vector * sol_vector = hypre_SeqVectorCreate( SSData.numRows );
    hypre_VectorMemoryLocation( sol_vector ) = HYPRE_MEMORY_HOST;
    hypre_SeqVectorInitialize( sol_vector );

    // solve Ax=b
    status = SuiteSparseSolveWorkingRank( SSData,
                                          hypre_VectorData( sol_vector ),
                                          hypre_VectorData( b_vector ),
                                          time );

    sol_Vector = ( HYPRE_Vector ) sol_vector;
  }

  // Copy partitioning otherwise hypre would point to data in x (issues with deallocation!!!)
  HYPRE_BigInt * partitioning = hypre_CTAlloc( HYPRE_BigInt, 2, HYPRE_MEMORY_HOST );
  partitioning[0] = hypre_ParVectorPartitioning( x.unwrapped() )[0];
  partitioning[1] = hypre_ParVectorPartitioning( x.unwrapped() )[1];

  HYPRE_ParVector sol_ParVector;
  GEOSX_LAI_CHECK_ERROR( HYPRE_VectorToParVector( x.getComm(),
                                                  sol_Vector,
                                                  partitioning,
                                                  &sol_ParVector ) );
  std::copy( hypre_VectorData( hypre_ParVectorLocalVector( sol_ParVector ) ),
             hypre_VectorData( hypre_ParVectorLocalVector( sol_ParVector ) ) + x.localSize(),
             x.extractLocalVector() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_VectorDestroy( sol_Vector ) );
  GEOSX_LAI_CHECK_ERROR( hypre_ParVectorDestroy( sol_ParVector ) );
  GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( b_vector ) );

  MpiWrapper::bcast( &time, 1, SSData.workingRank, SSData.comm );
  MpiWrapper::bcast( &status, 1, SSData.workingRank, SSData.comm );

  return status;
}

}
