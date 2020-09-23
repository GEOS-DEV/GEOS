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
                                 SuiteSparseData & SSData )
{
  // Copy distributed parcsr matrix in a local CSR matrix on every process
  // with at least one row
  // Warning: works for a parcsr matrix that is smaller than 2^31-1
  SSData.CSRmatrix = hypre_ParCSRMatrixToCSRMatrixAll( matrix.unwrapped() );

  // Identify the smallest process where CSRmatrix exists
  int rank = MpiWrapper::Comm_rank( matrix.getComm() );
  if( SSData.CSRmatrix == 0 )
  {
    rank = MpiWrapper::Comm_size( matrix.getComm() );
  }
  SSData.workingRank = MpiWrapper::Min( rank, matrix.getComm() );

  rank = MpiWrapper::Comm_rank( matrix.getComm() );
  if( rank == SSData.workingRank )
  {
    SSData.numRows = toSuiteSparse_Int( matrix.numGlobalRows() );
    SSData.numCols = toSuiteSparse_Int( matrix.numGlobalCols() );
    SSData.nonZeros = toSuiteSparse_Int( hypre_CSRMatrixNumNonzeros( SSData.CSRmatrix ) );
    HYPRE_Int const * const hypreI = hypre_CSRMatrixI( SSData.CSRmatrix );
    SSData.rowPtr = new Int[SSData.numRows+1];
    for( localIndex i = 0; i <= SSData.numRows; ++i )
    {
      SSData.rowPtr[i] = LvArray::integerConversion< Int >( hypreI[i] );
    }
    SSData.colIndices = new Int[SSData.nonZeros];
    SSData.data = new real64[SSData.nonZeros];
    HYPRE_Int const * const hypreJ = hypre_CSRMatrixJ( SSData.CSRmatrix );
    for( localIndex i = 0; i < SSData.nonZeros; ++i )
    {
      SSData.colIndices[i] = LvArray::integerConversion< Int >( hypreJ[i] );
    }
    std::copy( hypre_CSRMatrixData( SSData.CSRmatrix ),
               hypre_CSRMatrixData( SSData.CSRmatrix ) + SSData.nonZeros,
               SSData.data );
    for( localIndex i = 0; i < SSData.numRows; ++i )
    {
      localIndex rowLength = LvArray::integerConversion< localIndex >( SSData.rowPtr[i+1] - SSData.rowPtr[i] );
      sortIntReal64( &( SSData.colIndices[SSData.rowPtr[i]] ), &( SSData.data[SSData.rowPtr[i]] ), rowLength );
    }
  }
  else
  {
    // Destroy CSRmatrix
    if( SSData.CSRmatrix )
    {
      GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( SSData.CSRmatrix ) );
    }
  }

  // Save communicator
  SSData.comm = matrix.getComm();
}

/*--------------------------------------------------------------------------
 * hypre_VectorToParVector:
 * generates a ParVector from a Vector on proc 0 and distributes the pieces
 * to the other procs in comm
 *
 * this is not being optimized to use HYPRE_NO_GLOBAL_PARTITION
 *--------------------------------------------------------------------------*/

hypre_ParVector * HYPRE_VectorToParVectorGeneral ( MPI_Comm comm,
                                                   int proc,
                                                   HYPRE_Vector v,
                                                   HYPRE_BigInt *vec_starts )
{
  HYPRE_BigInt global_size;
  HYPRE_BigInt *global_vec_starts = NULL;
  HYPRE_BigInt first_index;
  HYPRE_BigInt last_index;
  HYPRE_Int local_size;
  HYPRE_Int num_vectors;
  HYPRE_Int num_procs, my_id;
  HYPRE_Int global_vecstride, vecstride, idxstride;
  hypre_ParVector *par_vector;
  hypre_Vector *local_vector;
  HYPRE_Complex *v_data;
  HYPRE_Complex *local_data;
  hypre_MPI_Request *requests;
  hypre_MPI_Status *status, status0;
  HYPRE_Int i, j, k, p;

  hypre_MPI_Comm_size( comm, &num_procs );
  hypre_MPI_Comm_rank( comm, &my_id );

  if( my_id == proc )
  {
    global_size = (HYPRE_BigInt)hypre_VectorSize( (hypre_Vector *) v );
    v_data = hypre_VectorData( (hypre_Vector *) v );
    num_vectors = hypre_VectorNumVectors( (hypre_Vector *) v ); /* for multivectors */
    global_vecstride = hypre_VectorVectorStride( (hypre_Vector *) v );
  }

  hypre_MPI_Bcast( &global_size, 1, HYPRE_MPI_INT, proc, comm );
  hypre_MPI_Bcast( &num_vectors, 1, HYPRE_MPI_INT, proc, comm );
  hypre_MPI_Bcast( &global_vecstride, 1, HYPRE_MPI_INT, proc, comm );

  if( num_vectors == 1 )
  {
    par_vector = hypre_ParVectorCreate( comm, global_size, vec_starts );
  }
  else
  {
    par_vector = hypre_ParMultiVectorCreate( comm, global_size, vec_starts, num_vectors );
  }

  vec_starts  = hypre_ParVectorPartitioning( par_vector );
  first_index = hypre_ParVectorFirstIndex( par_vector );
  last_index  = hypre_ParVectorLastIndex( par_vector );
  local_size  = (HYPRE_Int)(last_index - first_index) + 1;

#ifdef HYPRE_NO_GLOBAL_PARTITION
  if( my_id == proc )
  {
    global_vec_starts = hypre_CTAlloc( HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST );
  }
  hypre_MPI_Gather( &first_index, 1, HYPRE_MPI_BIG_INT, global_vec_starts,
                    1, HYPRE_MPI_BIG_INT, proc, comm );
  if( my_id == proc )
  {
    global_vec_starts[num_procs] = hypre_ParVectorGlobalSize( par_vector );
  }
#else
  global_vec_starts = vec_starts;
#endif

  hypre_ParVectorInitialize( par_vector );
  local_vector = hypre_ParVectorLocalVector( par_vector );
  local_data = hypre_VectorData( local_vector );
  vecstride = hypre_VectorVectorStride( local_vector );
  idxstride = hypre_VectorIndexStride( local_vector );
  /* so far the only implemented multivector StorageMethod is 0 */
  hypre_assert( idxstride==1 );

  if( my_id == proc )
  {
    requests = hypre_CTAlloc( hypre_MPI_Request, num_vectors*(num_procs-1), HYPRE_MEMORY_HOST );
    status = hypre_CTAlloc( hypre_MPI_Status, num_vectors*(num_procs-1), HYPRE_MEMORY_HOST );
    k = 0;
    for( p = 0; p<num_procs; p++ )
    {
      if( p != proc )
      {
        for( j = 0; j<num_vectors; ++j )
        {
          hypre_MPI_Isend( &v_data[(HYPRE_Int) global_vec_starts[p]] + j*global_vecstride,
                           (HYPRE_Int)(global_vec_starts[p+1] - global_vec_starts[p]),
                           HYPRE_MPI_COMPLEX, p, 0, comm, &requests[k++] );
        }
      }
    }
    if( num_vectors == 1 )
    {
      for( i = 0; i < local_size; i++ )
        local_data[i] = v_data[i];
    }
    else
    {
      for( j = 0; j<num_vectors; ++j )
      {
        for( i = 0; i < local_size; i++ )
        {
          local_data[i+j*vecstride] = v_data[i+j*global_vecstride];
        }
      }
    }
    hypre_MPI_Waitall( num_procs-1, requests, status );
    hypre_TFree( requests, HYPRE_MEMORY_HOST );
    hypre_TFree( status, HYPRE_MEMORY_HOST );
  }
  else
  {
    for( j=0; j<num_vectors; ++j )
    {
      hypre_MPI_Recv( local_data+j*vecstride, local_size, HYPRE_MPI_COMPLEX,
                      proc, 0, comm, &status0 );
    }
  }

#ifdef HYPRE_NO_GLOBAL_PARTITION
  if( global_vec_starts )
  {
    hypre_TFree( global_vec_starts, HYPRE_MEMORY_HOST );
  }
#endif

  return par_vector;
}
}

void SuiteSparseCreate( HypreMatrix const & matrix,
                        LinearSolverParameters const & params,
                        SuiteSparseData & SSData )
{
  // Get the default control parameters
  umfpack_dl_defaults( SSData.Control );
  SSData.Control[UMFPACK_PRL] = params.logLevel > 1 ? 6 : 1;
  SSData.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;

  // Convert matrix from Hypre to SuiteSparse format
  ConvertToSuiteSparseMatrix( matrix, SSData );
}

int SuiteSparseSetup( SuiteSparseData & SSData,
                      real64 & time )
{
  Stopwatch watch;

  int status = 0;

  int const rank = MpiWrapper::Comm_rank( SSData.comm );
  if( rank == SSData.workingRank )
  {
    // symbolic factorization
    status = umfpack_dl_symbolic( SSData.numCols,
                                  SSData.numRows,
                                  SSData.rowPtr,
                                  SSData.colIndices,
                                  SSData.data,
                                  &SSData.Symbolic,
                                  SSData.Control,
                                  SSData.Info );
    if( status < 0 )
    {
      umfpack_dl_report_info( SSData.Control, SSData.Info );
      umfpack_dl_report_status( SSData.Control, status );
      GEOSX_ERROR( "Hypre SuiteSparse interface: umfpack_dl_symbolic failed." );
    }

    // print the symbolic factorization
    if( SSData.logLevel > 1 )
    {
      umfpack_dl_report_symbolic( SSData.Symbolic, SSData.Control );
    }

    // numeric factorization
    status = umfpack_dl_numeric( SSData.rowPtr,
                                 SSData.colIndices,
                                 SSData.data,
                                 SSData.Symbolic,
                                 &SSData.Numeric,
                                 SSData.Control,
                                 SSData.Info );
    if( status < 0 )
    {
      umfpack_dl_report_info( SSData.Control, SSData.Info );
      umfpack_dl_report_status( SSData.Control, status );
      GEOSX_ERROR( "Hypre SuiteSparse interface: umfpack_dl_numeric failed." );
    }

    // print the numeric factorization
    if( SSData.logLevel > 1 )
    {
      umfpack_dl_report_numeric( SSData.Symbolic, SSData.Control );
    }
  }
  time = watch.elapsedTime();
  MpiWrapper::bcast( &time, 1, SSData.workingRank, SSData.comm );

  return status;
}

int SuiteSparseSolve( SuiteSparseData & SSData,
                      HypreVector const & b,
                      HypreVector & x,
                      real64 & time )
{
  Stopwatch watch;

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
    status = umfpack_dl_solve( UMFPACK_At,
                               SSData.rowPtr,
                               SSData.colIndices,
                               SSData.data,
                               hypre_VectorData( sol_vector ),
                               hypre_VectorData( b_vector ),
                               SSData.Numeric,
                               SSData.Control,
                               SSData.Info );

    if( status < 0 )
    {
      umfpack_dl_report_info( SSData.Control, SSData.Info );
      umfpack_dl_report_status( SSData.Control, status );
      GEOSX_ERROR( "Hypre SuiteSparse interface: umfpack_dl_solve failed." );
    }
    sol_Vector = ( HYPRE_Vector ) sol_vector;
  }

  // Copy partitioning otherwise hypre would point to data in x (issues with deallocation!!!)
  HYPRE_BigInt * partitioning = hypre_CTAlloc( HYPRE_BigInt, 2, HYPRE_MEMORY_HOST );
  partitioning[0] = hypre_ParVectorPartitioning( x.unwrapped() )[0];
  partitioning[1] = hypre_ParVectorPartitioning( x.unwrapped() )[1];

  HYPRE_ParVector sol_ParVector;
  // hypre assumes that the vector to scatter is ALWAYS on processor 0
  //GEOSX_LAI_CHECK_ERROR( HYPRE_VectorToParVector( x.getComm(),
  //                                                sol_Vector,
  //                                                partitioning,
  //                                                &sol_ParVector ) );

  // Same routine, with possibility to specify the root processor
  sol_ParVector = HYPRE_VectorToParVectorGeneral( x.getComm(),
                                                  SSData.workingRank,
                                                  sol_Vector,
                                                  partitioning );
  std::copy( hypre_VectorData( hypre_ParVectorLocalVector( sol_ParVector ) ),
             hypre_VectorData( hypre_ParVectorLocalVector( sol_ParVector ) ) + x.localSize(),
             x.extractLocalVector() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_VectorDestroy( sol_Vector ) );
  GEOSX_LAI_CHECK_ERROR( hypre_ParVectorDestroy( sol_ParVector ) );
  GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( b_vector ) );

  time = watch.elapsedTime();
  MpiWrapper::bcast( &time, 1, SSData.workingRank, SSData.comm );

  return status;
}

real64 SuiteSparseCondEst( SuiteSparseData const & SSData )
{
  return 1.0 / SSData.Info[UMFPACK_RCOND];
}

void SuiteSparseDestroy( SuiteSparseData & SSData )
{
  int const rank = MpiWrapper::Comm_rank( SSData.comm );
  if( rank == SSData.workingRank )
  {
    umfpack_dl_free_symbolic( &SSData.Symbolic );
    umfpack_dl_free_numeric( &SSData.Numeric );
    delete [] SSData.rowPtr;
    delete [] SSData.colIndices;
    delete [] SSData.data;
    GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( SSData.CSRmatrix ) );
  }
}

}
