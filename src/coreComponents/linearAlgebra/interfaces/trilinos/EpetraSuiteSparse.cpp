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

#include <Epetra_Import.h>

#ifdef GEOSX_USE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
using Epetra_MpiComm = Epetra_SerialComm;
#endif

namespace geosx
{

void ConvertEpetraToSuiteSparseMatrix( EpetraMatrix const & matrix,
                                       SuiteSparseData & SSData,
                                       Epetra_Map * & serialMap,
                                       Epetra_Import * & importToSerial )
{
  // Perform everything on rank 0
  SSData.workingRank = 0;

  MPI_Comm const comm = matrix.getComm();

  const Epetra_Map & originalMap = matrix.unwrapped().DomainMap();
  globalIndex const numGlobalElements = matrix.numGlobalRows();

  int const rank = MpiWrapper::Comm_rank( matrix.getComm() );

  globalIndex numMyElements = 0;
  if( rank == SSData.workingRank )
  {
    numMyElements = numGlobalElements;
  }

  //  Convert Original Matrix to Serial
  serialMap = new Epetra_Map( numGlobalElements, numMyElements, 0, Epetra_MpiComm( MPI_PARAM( comm ) ));
  importToSerial = new Epetra_Import( *serialMap, originalMap );
  Epetra_CrsMatrix serialCrsMatrix( Copy, *serialMap, 0 );
  serialCrsMatrix.Import( matrix.unwrapped(), *importToSerial, Insert );
  serialCrsMatrix.FillComplete();

  if( rank == SSData.workingRank )
  {
    globalIndex const nonZeros = matrix.numGlobalNonzeros();

    SSData.numRows = toSuiteSparse_Int( matrix.numGlobalRows() );
    SSData.numCols = toSuiteSparse_Int( matrix.numGlobalCols() );
    SSData.nonZeros = toSuiteSparse_Int( nonZeros );

    SSData.rowPtr = new Int[SSData.numRows+1];
    SSData.colIndices = new Int[SSData.nonZeros];
    SSData.data = new real64[SSData.nonZeros];

    array1d< int > ai( SSData.nonZeros );

    int const numEntries = serialCrsMatrix.MaxNumEntries();
    int numEntriesThisRow;
    int aiIndex = 0;
    for( localIndex i = 0; i < SSData.numRows; ++i )
    {
      SSData.rowPtr[i] = LvArray::integerConversion< Int >( aiIndex );
      GEOSX_LAI_CHECK_ERROR( serialCrsMatrix.ExtractMyRowCopy( i, numEntries, numEntriesThisRow, &SSData.data[aiIndex], &ai[aiIndex] ) );
      aiIndex += numEntriesThisRow;
    }
    SSData.rowPtr[SSData.numRows] = LvArray::integerConversion< Int >( aiIndex );

    for( localIndex i = 0; i < SSData.nonZeros; ++i )
    {
      SSData.colIndices[i] = LvArray::integerConversion< Int >( ai[i] );
    }
  }

  // Save communicator
  SSData.comm = matrix.getComm();
}

int SuiteSparseSolve( SuiteSparseData & SSData,
                      Epetra_Map const * serialMap,
                      Epetra_Import const * importToSerial,
                      EpetraVector const & b,
                      EpetraVector & x,
                      real64 & time )
{
  int status = 0;

  Epetra_MultiVector serialX( *serialMap, 1 );
  Epetra_MultiVector serialB( *serialMap, 1 );

  serialB.Import( b.unwrapped(), *importToSerial, Insert );

  int const rank = MpiWrapper::Comm_rank( SSData.comm );
  if( rank == SSData.workingRank )
  {
    // Extract Serial versions of X and B
    real64 *serialXvalues;
    real64 *serialBvalues;

    int serialBlda, serialXlda;
    GEOSX_LAI_CHECK_ERROR( serialB.ExtractView( &serialBvalues, &serialBlda ) );
    GEOSX_LAI_CHECK_ERROR( serialX.ExtractView( &serialXvalues, &serialXlda ) );
    GEOSX_ASSERT( serialBlda == LvArray::integerConversion< int >( SSData.numRows ) );
    GEOSX_ASSERT( serialXlda == LvArray::integerConversion< int >( SSData.numCols ) );

    // solve Ax=b
    status = SuiteSparseSolveWorkingRank( SSData,
                                          serialXvalues,
                                          serialBvalues,
                                          time );

  }

  x.unwrapped().Export( serialX, *importToSerial, Insert );

  MpiWrapper::bcast( &time, 1, SSData.workingRank, SSData.comm );
  MpiWrapper::bcast( &status, 1, SSData.workingRank, SSData.comm );

  return status;
}

}
