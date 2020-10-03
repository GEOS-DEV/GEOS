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

// Add two orders of magnitude to allow small error in condition number estimate
static real64 const machinePrecision = 100.0 * std::numeric_limits< real64 >::epsilon();

// Check matching requirements on index/value types between GEOSX and SuiteSparse

static_assert( sizeof( Int ) == sizeof( globalIndex ),
               "SuiteSparse Int and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< Int >::value == std::is_signed< globalIndex >::value,
               "SuiteSparse Int and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuiteSparse real and geosx::real64 must be the same type" );

SuiteSparse::SuiteSparse():
  m_numRows( 0 ),
  m_numCols( 0 ),
  m_nonZeros( 0 ),
  m_setupTime( 0 ),
  m_solveTime( 0 )
{}

SuiteSparse::SuiteSparse( LinearSolverParameters const & params ):
  m_numRows( 0 ),
  m_numCols( 0 ),
  m_nonZeros( 0 ),
  m_setupTime( 0 ),
  m_solveTime( 0 )
{
  create( params );
}

SuiteSparse::~SuiteSparse()
{
  destroy();
}

void SuiteSparse::create( LinearSolverParameters const & params )
{
  // Get the default control parameters
  umfpack_dl_defaults( m_Control );
  m_Control[UMFPACK_PRL] = params.logLevel > 1 ? 6 : 1;
  m_Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
}

int SuiteSparse::setup()
{
  Stopwatch watch;

  int status = 0;

  int const rank = MpiWrapper::Comm_rank( m_comm );
  if( rank == m_workingRank )
  {
    // symbolic factorization
    status = umfpack_dl_symbolic( m_numCols,
                                  m_numRows,
                                  m_rowPtr.data(),
                                  m_colIndices.data(),
                                  m_values.data(),
                                  &m_Symbolic,
                                  m_Control,
                                  m_Info );
    if( status < 0 )
    {
      umfpack_dl_report_info( m_Control, m_Info );
      umfpack_dl_report_status( m_Control, status );
      GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_symbolic failed." );
    }

    // print the symbolic factorization
    if( m_logLevel > 1 )
    {
      umfpack_dl_report_symbolic( m_Symbolic, m_Control );
    }

    // numeric factorization
    status = umfpack_dl_numeric( m_rowPtr.data(),
                                 m_colIndices.data(),
                                 m_values.data(),
                                 m_Symbolic,
                                 &m_Numeric,
                                 m_Control,
                                 m_Info );

    if( status < 0 )
    {
      umfpack_dl_report_info( m_Control, m_Info );
      umfpack_dl_report_status( m_Control, status );
      GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_numeric failed." );
    }

    // print the numeric factorization
    if( m_logLevel > 1 )
    {
      umfpack_dl_report_numeric( m_Symbolic, m_Control );
    }

    // save condition number
    m_condEst = 1.0 / m_Info[UMFPACK_RCOND];
  }
  MpiWrapper::bcast( &m_condEst, 1, m_workingRank, m_comm );

  m_setupTime = watch.elapsedTime();

  return status;
}

int SuiteSparse::solveWorkingRank( real64 * b, real64 * x )
{
  Stopwatch watch;

  int status = 0;
  // solve Ax=b
  status = umfpack_dl_solve( UMFPACK_At,
                             m_rowPtr.data(),
                             m_colIndices.data(),
                             m_values.data(),
                             b,
                             x,
                             m_Numeric,
                             m_Control,
                             m_Info );

  if( status < 0 )
  {
    umfpack_dl_report_info( m_Control, m_Info );
    umfpack_dl_report_status( m_Control, status );
    GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_solve failed." );
  }

  m_solveTime = watch.elapsedTime();

  return status;
}

void SuiteSparse::syncTimes()
{
  MpiWrapper::bcast( &m_solveTime, 1, m_workingRank, m_comm );
}

real64 SuiteSparse::condEst() const
{
  return m_condEst;
}

real64 SuiteSparse::relativeTolerance() const
{
  return m_condEst * machinePrecision;
}

void SuiteSparse::destroy()
{
  int const rank = MpiWrapper::Comm_rank( m_comm );
  if( rank == m_workingRank )
  {
    umfpack_dl_free_symbolic( &m_Symbolic );
    umfpack_dl_free_numeric( &m_Numeric );
  }
}

void SuiteSparse::setWorkingRank( int const workingRank )
{
  m_workingRank = workingRank;
}

int SuiteSparse::workingRank() const
{
  return m_workingRank;
}

void SuiteSparse::setComm( MPI_Comm const comm )
{
  m_comm = comm;
}

MPI_Comm SuiteSparse::getComm() const
{
  return m_comm;
}

void SuiteSparse::setNumRows( Int const numRows )
{
  m_numRows = numRows;
}

Int SuiteSparse::numRows() const
{
  return m_numRows;
}

void SuiteSparse::setNumCols( Int const numCols )
{
  m_numCols = numCols;
}

Int SuiteSparse::numCols() const
{
  return m_numCols;
}

void SuiteSparse::setNonZeros( Int const nonZeros )
{
  m_nonZeros = nonZeros;
}

Int SuiteSparse::nonZeros() const
{
  return m_nonZeros;
}

void SuiteSparse::createInternalStorage()
{
  m_rowPtr.resize( m_numRows + 1 );
  m_colIndices.resize( m_nonZeros );
  m_values.resize( m_nonZeros );
}

array1d< Int > & SuiteSparse::rowPtr()
{
  return m_rowPtr;
}

array1d< Int > & SuiteSparse::colIndices()
{
  return m_colIndices;
}

array1d< real64 > & SuiteSparse::values()
{
  return m_values;
}

real64 SuiteSparse::setupTime() const
{
  return m_setupTime;
}

real64 SuiteSparse::solveTime() const
{
  return m_solveTime;
}

}
