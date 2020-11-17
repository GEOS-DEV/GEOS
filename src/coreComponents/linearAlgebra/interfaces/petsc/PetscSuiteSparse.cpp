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
#include "linearAlgebra/interfaces/direct/Arnoldi.hpp"

#include "petsc.h"

namespace geosx
{

void ConvertPetscToSuiteSparseMatrix( PetscMatrix const & matrix,
                                      SuiteSparse & SSData )
{
  // Perform everything on rank 0
  int const workingRank = 0;
  SSData.setWorkingRank( workingRank );

  MPI_Comm const comm = matrix.getComm();

  int const rank = MpiWrapper::Comm_rank( comm );

  Mat * localMatrix;
  IS set;
  if( rank == workingRank )
  {
    GEOSX_LAI_CHECK_ERROR( ISCreateStride( PETSC_COMM_SELF, matrix.numGlobalRows(), 0, 1, &set ) );
  }
  GEOSX_LAI_CHECK_ERROR( MatCreateSubMatrices( matrix.unwrapped(), rank==workingRank, &set, &set, MAT_INITIAL_MATRIX, &localMatrix ) );

  if( rank == workingRank )
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

    SSData.resize( toSuiteSparse_Int( matrix.numGlobalRows() ),
                   toSuiteSparse_Int( matrix.numGlobalCols() ),
                   toSuiteSparse_Int( info.nz_used ) );

    for( localIndex i = 0; i <= SSData.numRows(); ++i )
    {
      SSData.rowPtr()[i] = LvArray::integerConversion< SSInt >( ia[i] );
    }
    for( localIndex i = 0; i < SSData.nonZeros(); ++i )
    {
      SSData.colIndices()[i] = LvArray::integerConversion< SSInt >( ja[i] );
    }
    std::copy( array, array + SSData.nonZeros(), SSData.values().data() );

    GEOSX_LAI_CHECK_ERROR( MatRestoreRowIJ( localMatrix[0], 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
    if( !status )
    {
      GEOSX_ERROR( "ConvertPetscToSuiteSparseMatrix: MatGetRowIJ reported an error." );
    }
    GEOSX_LAI_CHECK_ERROR( MatSeqAIJRestoreArray( localMatrix[0], &array ) );
  }

  // Destroy localMatrix
  GEOSX_LAI_CHECK_ERROR( MatDestroySubMatrices( rank==workingRank, &localMatrix ) );

  // Save communicator
  SSData.setComm( matrix.getComm() );
}

int SuiteSparseSolve( SuiteSparse & SSData,
                      PetscVector const & b,
                      PetscVector & x,
                      bool transpose )
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

  int const rank = MpiWrapper::Comm_rank( SSData.getComm() );
  if( rank == SSData.workingRank() )
  {
    real64 * dataRhs;
    real64 * dataSol;
    GEOSX_LAI_CHECK_ERROR( VecGetArray( localRhs, &dataRhs ) );
    GEOSX_LAI_CHECK_ERROR( VecGetArray( localSol, &dataSol ) );

    // solve Ax=b
    status = SSData.solveWorkingRank( dataSol, dataRhs, transpose );

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
  SSData.syncTimes();
  MpiWrapper::bcast( &status, 1, SSData.workingRank(), SSData.getComm() );

  return status;
}

namespace
{
class InverseOperator : public LinearOperator< PetscVector >
{
public:

  void set( PetscMatrix const & matrix, SuiteSparse & SSData )
  {
    m_SSData = &SSData;
    m_comm = SSData.getComm();
    m_numGlobalRows = matrix.numGlobalRows();
    m_numGlobalCols = matrix.numGlobalCols();
    m_numLocalRows = matrix.numLocalRows();
  }

  globalIndex numGlobalRows() const override
  {
    return m_numGlobalRows;
  }

  globalIndex numGlobalCols() const override
  {
    return m_numGlobalCols;
  }

  localIndex numLocalRows() const
  {
    return m_numLocalRows;
  }

  MPI_Comm const & getComm() const
  {
    return m_comm;
  }

  void apply( PetscVector const & x, PetscVector & y ) const override
  {
    SuiteSparseSolve( *m_SSData, x, y, false );
    SuiteSparseSolve( *m_SSData, y, y, true );
  }

private:

  SuiteSparse * m_SSData;

  MPI_Comm m_comm;

  globalIndex m_numGlobalRows;

  globalIndex m_numGlobalCols;

  localIndex m_numLocalRows;
};
}

real64 PetscSuiteSparseCond( PetscMatrix const & matrix, SuiteSparse & SSData )
{
  localIndex const numIterations = 4;

  using DirectOperator = DirectOperator< PetscMatrix, PetscVector >;
  DirectOperator directOperator;
  directOperator.set( matrix, matrix.getComm() );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue( directOperator, numIterations );

  InverseOperator inverseOperator;
  inverseOperator.set( matrix, SSData );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue( inverseOperator, numIterations );

  return sqrt( lambdaDirect * lambdaInverse );
}

}
