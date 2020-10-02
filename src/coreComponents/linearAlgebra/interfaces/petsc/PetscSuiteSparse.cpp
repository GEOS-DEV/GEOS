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
 * @file PetscSuiteSparse.cpp
 */

#include "PetscSuiteSparse.hpp"
#include "common/Stopwatch.hpp"

#include "petsc.h"

namespace geosx
{

void ConvertPetscToSuiteSparseMatrix( PetscMatrix const & matrix,
                                      SuiteSparseData & SSData )
{
  // Perform everything on rank 0
  SSData.workingRank = 0;

  MPI_Comm const comm = matrix.getComm();

  int const rank = MpiWrapper::Comm_rank( comm );

  Mat * localMatrix;
  IS set;
  if( rank == SSData.workingRank )
  {
    GEOSX_LAI_CHECK_ERROR( ISCreateStride( PETSC_COMM_SELF, matrix.numGlobalRows(), 0, 1, &set ) );
  }
  GEOSX_LAI_CHECK_ERROR( MatCreateSubMatrices( matrix.unwrapped(), rank==SSData.workingRank, &set, &set, MAT_INITIAL_MATRIX, &localMatrix ) );

  if( rank == SSData.workingRank )
  {
    PetscInt numRows;
    const PetscInt * ia;
    const PetscInt * ja;
    PetscBool status;
    GEOSX_LAI_CHECK_ERROR( MatGetRowIJ( localMatrix[0], 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
    if( !status )
    {
      GEOSX_ERROR( "ConvertPetscToSuiteSparseMatrix: MatGetRowIJ reported an error." );
    }

    real64 * array;
    GEOSX_LAI_CHECK_ERROR( MatSeqAIJGetArray( localMatrix[0], &array ) );

    MatInfo info;
    GEOSX_LAI_CHECK_ERROR( MatGetInfo( localMatrix[0], MAT_LOCAL, &info ) );

    SSData.numRows = toSuiteSparse_Int( matrix.numGlobalRows() );
    SSData.numCols = toSuiteSparse_Int( matrix.numGlobalCols() );
    SSData.nonZeros = toSuiteSparse_Int( info.nz_used );
    SSData.rowPtr = new Int[SSData.numRows+1];
    for( localIndex i = 0; i <= SSData.numRows; ++i )
    {
      SSData.rowPtr[i] = LvArray::integerConversion< Int >( ia[i] );
    }
    SSData.colIndices = new Int[SSData.nonZeros];
    SSData.data = new real64[SSData.nonZeros];
    for( localIndex i = 0; i < SSData.nonZeros; ++i )
    {
      SSData.colIndices[i] = LvArray::integerConversion< Int >( ja[i] );
    }
    std::copy( array, array + SSData.nonZeros, SSData.data );

    GEOSX_LAI_CHECK_ERROR( MatRestoreRowIJ( localMatrix[0], 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
    if( !status )
    {
      GEOSX_ERROR( "ConvertPetscToSuiteSparseMatrix: MatGetRowIJ reported an error." );
    }
    GEOSX_LAI_CHECK_ERROR( MatSeqAIJRestoreArray( localMatrix[0], &array ) );
  }

  // Destroy localMatrix
  GEOSX_LAI_CHECK_ERROR( MatDestroySubMatrices( rank==0, &localMatrix ) );

  // Save communicator
  SSData.comm = matrix.getComm();
}

int SuiteSparseSolve( SuiteSparseData & SSData,
                      PetscVector const & b,
                      PetscVector & x,
                      real64 & time )
{
  int status = 0;

  VecScatter rhsScatter;
  Vec localRhs;
  GEOSX_LAI_CHECK_ERROR( VecScatterCreateToZero( b.unwrapped(), &rhsScatter, &localRhs ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterBegin( rhsScatter, b.unwrapped(), localRhs, INSERT_VALUES, SCATTER_FORWARD ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterEnd( rhsScatter, b.unwrapped(), localRhs, INSERT_VALUES, SCATTER_FORWARD ) );

  VecScatter solScatter;
  Vec localSol;
  GEOSX_LAI_CHECK_ERROR( VecScatterCreateToZero( x.unwrapped(), &solScatter, &localSol ) );

  int const rank = MpiWrapper::Comm_rank( SSData.comm );
  if( rank == SSData.workingRank )
  {
    real64 * dataRhs;
    real64 * dataSol;
    GEOSX_LAI_CHECK_ERROR( VecGetArray( localRhs, &dataRhs ) );
    GEOSX_LAI_CHECK_ERROR( VecGetArray( localSol, &dataSol ) );

    // solve Ax=b
    status = SuiteSparseSolveWorkingRank( SSData,
                                          dataSol,
                                          dataRhs,
                                          time );

    GEOSX_LAI_CHECK_ERROR( VecRestoreArray( localRhs, &dataRhs ) );
    GEOSX_LAI_CHECK_ERROR( VecRestoreArray( localSol, &dataSol ) );
  }

  GEOSX_LAI_CHECK_ERROR( VecScatterBegin( solScatter, localSol, x.unwrapped(), INSERT_VALUES, SCATTER_REVERSE ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterEnd( solScatter, localSol, x.unwrapped(), INSERT_VALUES, SCATTER_REVERSE ) );
  GEOSX_LAI_CHECK_ERROR( VecDestroy( &localRhs ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterDestroy( &rhsScatter ) );
  GEOSX_LAI_CHECK_ERROR( VecDestroy( &localSol ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterDestroy( &solScatter ) );

  // broadcast status and times
  MpiWrapper::bcast( &time, 1, SSData.workingRank, SSData.comm );
  MpiWrapper::bcast( &status, 1, SSData.workingRank, SSData.comm );

  return status;
}

}
