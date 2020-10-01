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
 * @file HypreSuperlu.cpp
 */

#if !defined(DLONG)
#define DLONG
#endif

#include <numeric>

#include "HypreSuiteSparse.hpp"
#include "common/Stopwatch.hpp"

#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"
#include "seq_mv.h"

namespace geosx
{

// Check matching requirements on index/value types between GEOSX and SuiteSparse

static_assert( sizeof( Int ) == sizeof( globalIndex ),
               "SuiteSparse Int and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< Int >::value == std::is_signed< globalIndex >::value,
               "SuiteSparse Int and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuiteSparse real and geosx::real64 must be the same type" );

namespace
{

/**
 * @brief Convert GEOSX globalIndex value to SuiteSparse Int
 * @param index the input value
 * @return the converted value
 */
inline Int toSuiteSparse_Int( globalIndex const index )
{
  return LvArray::integerConversion< Int >( index );
}

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

/**
 * @brief Converts a matrix from Hypre to SuiteSparse format
 * @param[in] matrix the HypreMatrix object
 * @param[out] SSData the structure containing the matrix in SuiteSparse format
 */
void ConvertToSuiteSparseMatrix( HypreMatrix const & matrix,
                                 SuiteSparseData & ssData )
{
  // Copy distributed parcsr matrix in a local CSR matrix on every process
  // with at least one row
  // Warning: works for a parcsr matrix that is smaller than 2^31-1
  ssData.CSRmatrix = hypre_ParCSRMatrixToCSRMatrixAll( matrix.unwrapped() );

  // Identify the smallest process where CSRmatrix exists
  int rank = MpiWrapper::Comm_rank( matrix.getComm() );
  if( ssData.CSRmatrix == 0 )
  {
    rank = MpiWrapper::Comm_size( matrix.getComm() );
  }
  ssData.workingRank = MpiWrapper::Min( rank, matrix.getComm() );

  rank = MpiWrapper::Comm_rank( matrix.getComm() );
  if( rank == ssData.workingRank )
  {
    ssData.numRows = toSuiteSparse_Int( matrix.numGlobalRows() );
    ssData.numCols = toSuiteSparse_Int( matrix.numGlobalCols() );
    ssData.nonZeros = toSuiteSparse_Int( hypre_CSRMatrixNumNonzeros( ssData.CSRmatrix ) );
    HYPRE_Int const * const hypreI = hypre_CSRMatrixI( ssData.CSRmatrix );
    ssData.rowPtr = new Int[ssData.numRows+1];
    for( localIndex i = 0; i <= ssData.numRows; ++i )
    {
      ssData.rowPtr[i] = LvArray::integerConversion< Int >( hypreI[i] );
    }
    ssData.colIndices = new Int[ssData.nonZeros];
    ssData.data = new real64[ssData.nonZeros];
    HYPRE_Int const * const hypreJ = hypre_CSRMatrixJ( ssData.CSRmatrix );
    for( localIndex i = 0; i < ssData.nonZeros; ++i )
    {
      ssData.colIndices[i] = LvArray::integerConversion< Int >( hypreJ[i] );
    }
    std::copy( hypre_CSRMatrixData( ssData.CSRmatrix ),
               hypre_CSRMatrixData( ssData.CSRmatrix ) + ssData.nonZeros,
               ssData.data );
    for( localIndex i = 0; i < ssData.numRows; ++i )
    {
      localIndex rowLength = LvArray::integerConversion< localIndex >( ssData.rowPtr[i+1] - ssData.rowPtr[i] );
      sortIntReal64( &( ssData.colIndices[ssData.rowPtr[i]] ), &( ssData.data[ssData.rowPtr[i]] ), rowLength );
    }
  }
  else
  {
    // Destroy CSRmatrix
    if( ssData.CSRmatrix )
    {
      GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( ssData.CSRmatrix ) );
    }
  }

  // Save communicator
  ssData.comm = matrix.getComm();
}
}

void SuiteSparseCreate( HypreMatrix const & matrix,
                        LinearSolverParameters const & params,
                        SuiteSparseData & ssData )
{
  // Get the default control parameters
  umfpack_dl_defaults( ssData.Control );
  ssData.Control[UMFPACK_PRL] = params.logLevel > 1 ? 6 : 1;
  ssData.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;

  // Convert matrix from Hypre to SuiteSparse format
  ConvertToSuiteSparseMatrix( matrix, ssData );
}

int SuiteSparseSetup( SuiteSparseData & ssData,
                      real64 & time )
{
  Stopwatch watch;

  int status = 0;

  int const rank = MpiWrapper::Comm_rank( ssData.comm );
  if( rank == ssData.workingRank )
  {
    // symbolic factorization
    status = umfpack_dl_symbolic( ssData.numCols,
                                  ssData.numRows,
                                  ssData.rowPtr,
                                  ssData.colIndices,
                                  ssData.data,
                                  &ssData.Symbolic,
                                  ssData.Control,
                                  ssData.Info );
    if( status < 0 )
    {
      umfpack_dl_report_info( ssData.Control, ssData.Info );
      umfpack_dl_report_status( ssData.Control, status );
      GEOSX_ERROR( "Hypre SuiteSparse interface: umfpack_dl_symbolic failed." );
    }

    // print the symbolic factorization
    if( ssData.logLevel > 1 )
    {
      umfpack_dl_report_symbolic( ssData.Symbolic, ssData.Control );
    }

    // numeric factorization
    status = umfpack_dl_numeric( ssData.rowPtr,
                                 ssData.colIndices,
                                 ssData.data,
                                 ssData.Symbolic,
                                 &ssData.Numeric,
                                 ssData.Control,
                                 ssData.Info );
    if( status < 0 )
    {
      umfpack_dl_report_info( ssData.Control, ssData.Info );
      umfpack_dl_report_status( ssData.Control, status );
      GEOSX_ERROR( "Hypre SuiteSparse interface: umfpack_dl_numeric failed." );
    }

    // print the numeric factorization
    if( ssData.logLevel > 1 )
    {
      umfpack_dl_report_numeric( ssData.Symbolic, ssData.Control );
    }
  }

  time = watch.elapsedTime();
  MpiWrapper::bcast( &time, 1, ssData.workingRank, ssData.comm );

  return status;
}

int SuiteSparseSolve( SuiteSparseData & ssData,
                      HypreVector const & b,
                      HypreVector & x,
                      real64 & time )
{
  Stopwatch watch;

  int status = 0;

  hypre_Vector * b_vector = hypre_ParVectorToVectorAll( b.unwrapped() );
  HYPRE_Vector sol_Vector = nullptr;

  int const rank = MpiWrapper::Comm_rank( ssData.comm );
  if( rank == ssData.workingRank )
  {
    // Create local vector to store the solution
    hypre_Vector * sol_vector = hypre_SeqVectorCreate( ssData.numRows );
    hypre_VectorMemoryLocation( sol_vector ) = HYPRE_MEMORY_HOST;
    hypre_SeqVectorInitialize( sol_vector );

    // solve Ax=b
    status = umfpack_dl_solve( UMFPACK_At,
                               ssData.rowPtr,
                               ssData.colIndices,
                               ssData.data,
                               hypre_VectorData( sol_vector ),
                               hypre_VectorData( b_vector ),
                               ssData.Numeric,
                               ssData.Control,
                               ssData.Info );

    if( status < 0 )
    {
      umfpack_dl_report_info( ssData.Control, ssData.Info );
      umfpack_dl_report_status( ssData.Control, status );
      GEOSX_ERROR( "Hypre SuiteSparse interface: umfpack_dl_solve failed." );
    }
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

  time = watch.elapsedTime();
  MpiWrapper::bcast( &time, 1, ssData.workingRank, ssData.comm );

  return status;
}

void SuiteSparseDestroy( SuiteSparseData & ssData )
{
  int const rank = MpiWrapper::Comm_rank( ssData.comm );
  if( rank == ssData.workingRank )
  {
    umfpack_dl_free_symbolic( &ssData.Symbolic );
    umfpack_dl_free_numeric( &ssData.Numeric );
    delete [] ssData.rowPtr;
    delete [] ssData.colIndices;
    delete [] ssData.data;
    GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( ssData.CSRmatrix ) );
  }
}

}
