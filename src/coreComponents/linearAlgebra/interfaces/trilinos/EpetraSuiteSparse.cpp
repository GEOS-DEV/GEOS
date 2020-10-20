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
 * @file EpetraSuiteSparse.cpp
 */

#include "EpetraSuiteSparse.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/direct/Arnoldi.hpp"

#include <Epetra_Import.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>

#ifdef GEOSX_USE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
using Epetra_MpiComm = Epetra_SerialComm;
#endif

namespace geosx
{

void ConvertEpetraToSuiteSparseMatrix( EpetraMatrix const & matrix,
                                       SuiteSparse & SSData,
                                       Epetra_Map * & serialMap,
                                       Epetra_Import * & importToSerial )
{
  // Perform everything on rank 0
  int const workingRank = 0;
  SSData.setWorkingRank( workingRank );

  MPI_Comm const comm = matrix.getComm();

  const Epetra_Map & originalMap = matrix.unwrapped().DomainMap();
  globalIndex const numGlobalElements = matrix.numGlobalRows();

  int const rank = MpiWrapper::Comm_rank( matrix.getComm() );

  globalIndex numMyElements = 0;
  if( rank == workingRank )
  {
    numMyElements = numGlobalElements;
  }

  //  Convert Original Matrix to Serial
  serialMap = new Epetra_Map( numGlobalElements, numMyElements, 0, Epetra_MpiComm( MPI_PARAM( comm ) ));
  importToSerial = new Epetra_Import( *serialMap, originalMap );
  Epetra_CrsMatrix serialCrsMatrix( Copy, *serialMap, 0 );
  serialCrsMatrix.Import( matrix.unwrapped(), *importToSerial, Insert );
  serialCrsMatrix.FillComplete();

  if( rank == workingRank )
  {
    globalIndex const nonZeros = matrix.numGlobalNonzeros();

    SSData.setNumRows( toSuiteSparse_Int( matrix.numGlobalRows() ) );
    SSData.setNumCols( toSuiteSparse_Int( matrix.numGlobalCols() ) );
    SSData.setNonZeros( toSuiteSparse_Int( nonZeros ) );

    SSData.createInternalStorage();

    int const numEntries = serialCrsMatrix.MaxNumEntries();
    int numEntriesThisRow;
    int aiIndex = 0;
    array1d< int > ai( SSData.nonZeros() );
    for( localIndex i = 0; i < SSData.numRows(); ++i )
    {
      SSData.rowPtr()[i] = LvArray::integerConversion< Int >( aiIndex );
      GEOSX_LAI_CHECK_ERROR( serialCrsMatrix.ExtractMyRowCopy( i, numEntries, numEntriesThisRow, &SSData.values()[aiIndex], &ai[aiIndex] ) );
      aiIndex += numEntriesThisRow;
    }
    SSData.rowPtr()[SSData.numRows()] = LvArray::integerConversion< Int >( aiIndex );

    for( localIndex i = 0; i < SSData.nonZeros(); ++i )
    {
      SSData.colIndices()[i] = LvArray::integerConversion< Int >( ai[i] );
    }
  }

  // Save communicator
  SSData.setComm( matrix.getComm() );
}

int SuiteSparseSolve( SuiteSparse & SSData,
                      Epetra_Map const * serialMap,
                      Epetra_Import const * importToSerial,
                      EpetraVector const & b,
                      EpetraVector & x,
                      bool transpose )
{
  int status = 0;

  Epetra_MultiVector serialX( *serialMap, 1 );
  Epetra_MultiVector serialB( *serialMap, 1 );

  serialB.Import( b.unwrapped(), *importToSerial, Insert );

  int const rank = MpiWrapper::Comm_rank( SSData.getComm() );
  if( rank == SSData.workingRank() )
  {
    // Extract Serial versions of X and B
    real64 *serialXvalues;
    real64 *serialBvalues;

    int serialBlda, serialXlda;
    GEOSX_LAI_CHECK_ERROR( serialB.ExtractView( &serialBvalues, &serialBlda ) );
    GEOSX_LAI_CHECK_ERROR( serialX.ExtractView( &serialXvalues, &serialXlda ) );
    GEOSX_ASSERT( serialBlda == LvArray::integerConversion< int >( SSData.numRows() ) );
    GEOSX_ASSERT( serialXlda == LvArray::integerConversion< int >( SSData.numCols() ) );

    // solve Ax=b
    status = SSData.solveWorkingRank( serialXvalues, serialBvalues, transpose );
  }

  x.unwrapped().Export( serialX, *importToSerial, Insert );

  SSData.syncTimes();
  MpiWrapper::bcast( &status, 1, SSData.workingRank(), SSData.getComm() );

  return status;
}

namespace
{
class InverseOperator
{
public:

  void set( EpetraMatrix const & matrix,
            Epetra_Map const * serialMap,
            Epetra_Import const * importToSerial,
            SuiteSparse & SSData )
  {
    m_SSData = &SSData;
    m_comm = SSData.getComm();
    m_serialMap = serialMap;
    m_importToSerial = importToSerial;
    m_numGlobalRows = matrix.numGlobalRows();
    m_numLocalRows = matrix.numLocalRows();
  }

  globalIndex globalSize() const
  {
    return m_numGlobalRows;
  }

  localIndex localSize() const
  {
    return m_numLocalRows;
  }

  MPI_Comm const & getComm() const
  {
    return m_comm;
  }

  void apply( EpetraVector const & x, EpetraVector & y ) const
  {
    SuiteSparseSolve( *m_SSData, m_serialMap, m_importToSerial, x, y, false );
    SuiteSparseSolve( *m_SSData, m_serialMap, m_importToSerial, y, y, true );
  }

private:

  SuiteSparse * m_SSData;

  Epetra_Map const * m_serialMap;

  Epetra_Import const * m_importToSerial;

  MPI_Comm m_comm;

  globalIndex m_numGlobalRows;

  localIndex m_numLocalRows;
};
}

real64 EpetraSuiteSparseCond( EpetraMatrix const & matrix,
                              Epetra_Map const * serialMap,
                              Epetra_Import const * importToSerial,
                              SuiteSparse & SSData )
{
  localIndex const numIterations = 4;

  using DirectOperator = DirectOperator< EpetraMatrix, EpetraVector >;
  DirectOperator directOperator;
  directOperator.set( matrix );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue< EpetraVector, DirectOperator >( directOperator, numIterations );

  InverseOperator inverseOperator;
  inverseOperator.set( matrix, serialMap, importToSerial, SSData );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue< EpetraVector, InverseOperator >( inverseOperator, numIterations );

  return sqrt( lambdaDirect * lambdaInverse );
}

}
