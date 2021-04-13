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
#include "common/MpiWrapper.hpp"
#include "common/Stopwatch.hpp"

namespace geosx
{

// Check matching requirements on index/value types between GEOSX and SuiteSparse

static_assert( sizeof( SSInt ) == sizeof( globalIndex ),
               "SuiteSparse Int and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< SSInt >::value == std::is_signed< globalIndex >::value,
               "SuiteSparse Int and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuiteSparse real and geosx::real64 must be the same type" );

SuiteSparse::SuiteSparse():
  m_numRows( 0 ),
  m_numCols( 0 ),
  m_nonZeros( 0 ),
  m_subComm( MPI_COMM_NULL ),
  m_usingSubComm( false ),
  m_setupTime( 0 ),
  m_solveTime( 0 )
{}

SuiteSparse::SuiteSparse( LinearSolverParameters const & params ):
  m_numRows( 0 ),
  m_numCols( 0 ),
  m_nonZeros( 0 ),
  m_subComm( MPI_COMM_NULL ),
  m_usingSubComm( false ),
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
  // No need of subcommunicator: just point to the global communicator
  if( !m_usingSubComm )
  {
    m_subComm = m_comm;
  }

  Stopwatch watch;

  int status = 0;

  if( m_subComm != MPI_COMM_NULL )
  {
    int const rank = MpiWrapper::commRank( m_subComm );
    if( rank == m_subCommWorkingRank )
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
  }
  MpiWrapper::bcast( &m_condEst, 1, m_workingRank, m_comm );

  m_setupTime = watch.elapsedTime();

  return status;
}

int SuiteSparse::solveWorkingRank( real64 * b, real64 * x, bool transpose )
{
  Stopwatch watch;

  int status = 0;
  // solve Ax=b
  if( !transpose )
  {
    status = umfpack_dl_solve( UMFPACK_At,
                               m_rowPtr.data(),
                               m_colIndices.data(),
                               m_values.data(),
                               b,
                               x,
                               m_Numeric,
                               m_Control,
                               m_Info );
  }
  else
  {
    status = umfpack_dl_solve( UMFPACK_A,
                               m_rowPtr.data(),
                               m_colIndices.data(),
                               m_values.data(),
                               b,
                               x,
                               m_Numeric,
                               m_Control,
                               m_Info );
  }

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
  return m_condEst * m_precisionTolerance;
}

void SuiteSparse::destroy()
{
  if( m_subComm != MPI_COMM_NULL )
  {
    int const rank = MpiWrapper::commRank( m_subComm );
    if( rank == m_subCommWorkingRank )
    {
      umfpack_dl_free_symbolic( &m_Symbolic );
      umfpack_dl_free_numeric( &m_Numeric );
    }
    if( m_subComm != m_comm )
    {
      MpiWrapper::commFree( m_subComm );
    }
  }
}

void SuiteSparse::setWorkingRank( int const workingRank )
{
  m_workingRank = workingRank;
  if( !m_usingSubComm )
  {
    m_subCommWorkingRank = workingRank;
  }
}

int SuiteSparse::workingRank() const
{
  return m_workingRank;
}

void SuiteSparse::setSubCommWorkingRank( int const subCommWorkingRank )
{
  m_subCommWorkingRank = subCommWorkingRank;
}

int SuiteSparse::subCommWorkingRank() const
{
  return m_subCommWorkingRank;
}

void SuiteSparse::setComm( MPI_Comm const comm )
{
  m_comm = comm;
}

MPI_Comm SuiteSparse::getComm() const
{
  return m_comm;
}

void SuiteSparse::setSubComm( MPI_Comm const subComm )
{
  m_subComm = subComm;
  m_usingSubComm = true;
}

MPI_Comm SuiteSparse::getSubComm() const
{
  return m_subComm;
}

SSInt SuiteSparse::numRows() const
{
  return m_numRows;
}

SSInt SuiteSparse::numCols() const
{
  return m_numCols;
}

SSInt SuiteSparse::nonZeros() const
{
  return m_nonZeros;
}

void SuiteSparse::resize( SSInt const numRows, SSInt const numCols, SSInt const nonZeros )
{
  m_numRows = numRows;
  m_numCols = numCols;
  m_nonZeros = nonZeros;
  m_rowPtr.resize( m_numRows + 1 );
  m_colIndices.resize( m_nonZeros );
  m_values.resize( m_nonZeros );
}

arrayView1d< SSInt > SuiteSparse::rowPtr()
{
  return m_rowPtr;
}

arrayView1d< SSInt > SuiteSparse::colIndices()
{
  return m_colIndices;
}

arrayView1d< real64 > SuiteSparse::values()
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

real64 SuiteSparse::precisionTolerance() const
{
  return m_precisionTolerance;
}

}
