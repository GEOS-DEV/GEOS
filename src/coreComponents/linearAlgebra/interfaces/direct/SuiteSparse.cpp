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
 * @file SuiteSparse.cpp
 */

#if !defined(DLONG)
#define DLONG
#endif

#include "SuiteSparse.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "common/Stopwatch.hpp"

namespace geosx
{

// Check matching requirements on index/value types between GEOSX and SuiteSparse

static_assert( sizeof( Int ) == sizeof( globalIndex ),
               "SuiteSparse Int and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< Int >::value == std::is_signed< globalIndex >::value,
               "SuiteSparse Int and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuiteSparse real and geosx::real64 must be the same type" );

void SuiteSparseCreate( LinearSolverParameters const & params,
                        SuiteSparseData & SSData )
{
  // Get the default control parameters
  umfpack_dl_defaults( SSData.Control );
  SSData.Control[UMFPACK_PRL] = params.logLevel > 1 ? 6 : 1;
  SSData.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
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
      GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_symbolic failed." );
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
      GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_numeric failed." );
    }

    // print the numeric factorization
    if( SSData.logLevel > 1 )
    {
      umfpack_dl_report_numeric( SSData.Symbolic, SSData.Control );
    }
  }

  time = watch.elapsedTime();

  return status;
}

int SuiteSparseSolveWorkingRank( SuiteSparseData & SSData,
                                 real64 * b,
                                 real64 * x,
                                 real64 & time )
{
  Stopwatch watch;

  int status = 0;
  // solve Ax=b
  status = umfpack_dl_solve( UMFPACK_At,
                             SSData.rowPtr,
                             SSData.colIndices,
                             SSData.data,
                             b,
                             x,
                             SSData.Numeric,
                             SSData.Control,
                             SSData.Info );

  if( status < 0 )
  {
    umfpack_dl_report_info( SSData.Control, SSData.Info );
    umfpack_dl_report_status( SSData.Control, status );
    GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_solve failed." );
  }

  time = watch.elapsedTime();

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
  }
}

}
