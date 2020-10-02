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
                                       Epetra_Map * SerialMap,
                                       Epetra_Import * ImportToSerial )
{
  // Perform everything on rank 0
  SSData.workingRank = 0;

  MPI_Comm const comm = matrix.getComm();

  //const Epetra_Map &OriginalMap = matrix.unwrapped().RowMatrixRowMap();
  const Epetra_Map & OriginalMap = matrix.unwrapped().DomainMap();
  int const NumGlobalElements_ = LvArray::integerConversion< int >( matrix.numGlobalRows() );

  int const rank = MpiWrapper::Comm_rank( matrix.getComm() );

  int NumMyElements_ = 0;
  if( rank == SSData.workingRank )
    NumMyElements_ = NumGlobalElements_;

  //  Convert Original Matrix to Serial
  SerialMap = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Epetra_MpiComm( MPI_PARAM( comm ) ));
  std::cout << "DOES NOT WORK!!!\n";
  ImportToSerial = new Epetra_Import( *SerialMap, OriginalMap );
  std::cout << "DOES NOT WORK!!! (int excpetion)\n";
  Epetra_CrsMatrix SerialCrsMatrix( Copy, *SerialMap, 0 );
  SerialCrsMatrix.Import( matrix.unwrapped(), *ImportToSerial, Insert );
  SerialCrsMatrix.FillComplete();

  if( rank == SSData.workingRank )
  {
    globalIndex const nonZeros = matrix.numGlobalNonzeros();

    SSData.numRows = toSuiteSparse_Int( matrix.numGlobalRows() );
    SSData.numCols = toSuiteSparse_Int( matrix.numGlobalCols() );
    SSData.nonZeros = toSuiteSparse_Int( nonZeros );

    SSData.rowPtr = new Int[SSData.numRows+1];
    SSData.colIndices = new Int[SSData.nonZeros];
    SSData.data = new real64[SSData.nonZeros];

    array1d< int > Ai( SSData.nonZeros );

    int const NumEntries = SerialCrsMatrix.MaxNumEntries();
    int NumEntriesThisRow;
    int Ai_index = 0;
    for( localIndex i = 0; i < SSData.numRows; ++i )
    {
      SSData.rowPtr[i] = LvArray::integerConversion< Int >( Ai_index );
      GEOSX_LAI_CHECK_ERROR( SerialCrsMatrix.ExtractMyRowCopy( i, NumEntries, NumEntriesThisRow, &SSData.data[Ai_index], &Ai[Ai_index] ) );
      Ai_index += NumEntriesThisRow;
    }
    SSData.rowPtr[SSData.numRows] = LvArray::integerConversion< Int >( Ai_index );

    for( localIndex i = 0; i < SSData.nonZeros; ++i )
    {
      SSData.colIndices[i] = LvArray::integerConversion< Int >( Ai[i] );
    }
  }

  // Save communicator
  SSData.comm = matrix.getComm();
}

int SuiteSparseSolve( SuiteSparseData & SSData,
                      Epetra_Map const * SerialMap,
                      Epetra_Import const * ImportToSerial,
                      EpetraVector const & b,
                      EpetraVector & x,
                      real64 & time )
{
  int status = 0;

  Epetra_MultiVector SerialX( *SerialMap, 1 );
  Epetra_MultiVector SerialB( *SerialMap, 1 );

  SerialB.Import( b.unwrapped(), *ImportToSerial, Insert );

  int SerialBlda, SerialXlda;
  std::cout << "HERE\n";

  int const rank = MpiWrapper::Comm_rank( SSData.comm );
  if( rank == SSData.workingRank )
  {
    //  Extract Serial versions of X and B
    real64 *SerialXvalues;
    real64 *SerialBvalues;

    GEOSX_LAI_CHECK_ERROR( SerialB.ExtractView( &SerialBvalues, &SerialBlda ) );
    GEOSX_LAI_CHECK_ERROR( SerialX.ExtractView( &SerialXvalues, &SerialXlda ) );
    GEOSX_ASSERT( SerialBlda == LvArray::integerConversion< int >( SSData.nonZeros ) );
    GEOSX_ASSERT( SerialXlda == LvArray::integerConversion< int >( SSData.nonZeros ) );

    // solve Ax=b
    status = SuiteSparseSolveWorkingRank( SSData,
                                          SerialXvalues,
                                          SerialBvalues,
                                          time );

  }

  x.unwrapped().Export( SerialX, *ImportToSerial, Insert );

  MpiWrapper::bcast( &time, 1, SSData.workingRank, SSData.comm );
  MpiWrapper::bcast( &status, 1, SSData.workingRank, SSData.comm );

  return status;
}

}
