/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SuiteSparse.cpp
 */

#include "SuiteSparse.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/common/common.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/Arnoldi.hpp"
#include "linearAlgebra/utilities/NormalOperator.hpp"
#include "linearAlgebra/utilities/InverseNormalOperator.hpp"
#include "common/MpiWrapper.hpp"

// Pre-define some suitesparse variables since they are not properly defined
// in the header for alternate index types.
//#if GEOSX_GLOBALINDEX_TYPE_FLAG==0
#if 0
/// Set alias for SuiteSparse_long
#define SuiteSparse_long int

/// Set value for SuiteSparse_long_max
#define SuiteSparse_long_max 2147483647

/// Set string key for SuiteSparse_long
#define SuiteSparse_long_idd "d"

/// Set printf key for SuiteSparse_long
#define SuiteSparse_long_id "%d"
#endif

#if !defined(DLONG)
#define DLONG
#endif

#include <umfpack.h>

namespace geosx
{

/**
 * SuiteSparse integer definition
 */
using SSlong = SuiteSparse_long;

// Check matching requirements on index/value types between GEOSX and SuiteSparse

//static_assert( sizeof( SSlong ) == sizeof( globalIndex ),
//               "SuiteSparse Int and geosx::globalIndex must have the same size" );
//
//static_assert( std::is_signed< SSlong >::value == std::is_signed< globalIndex >::value,
//               "SuiteSparse Int and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuiteSparse real and geosx::real64 must be the same type" );

struct SuiteSparseData
{
  array1d< SSlong > rowPtr{};          /// row pointers
  array1d< SSlong > colIndices{};      /// column indices
  array1d< double > values{};          /// values
  array1d< double > rhs{};             /// right-hand side vector
  array1d< double > sol{};             /// solution vector
  real64 info[UMFPACK_INFO]{};       /// data structure to gather various info
  real64 control[UMFPACK_CONTROL]{}; /// SuiteSparse options
  void * symbolic{};                 /// pointer to the symbolic factorization
  void * numeric{};                  /// pointer to the numeric factorization

  SuiteSparseData( SSlong const numRows,
                   SSlong const numNonzeros )
    : rowPtr( numRows + 1 ),
    colIndices( numNonzeros ),
    values( numNonzeros ),
    rhs( numRows ),
    sol( numRows )
  {}

  ~SuiteSparseData()
  {
    if( symbolic )
    {
      umfpack_dl_free_symbolic( &symbolic );
    }
    if( numeric )
    {
      umfpack_dl_free_numeric( &numeric );
    }
  }
};

template< typename LAI >
SuiteSparse< LAI >::SuiteSparse( LinearSolverParameters params )
  : Base( std::move( params ) ),
  m_workingRank( -1 ),
  m_condEst( -1.0 )
{}

template< typename LAI >
SuiteSparse< LAI >::~SuiteSparse() = default;

template< typename LAI >
void SuiteSparse< LAI >::setup( Matrix const & mat )
{
  clear();
  PreconditionerBase< LAI >::setup( mat );

  // Choose working rank that will carry out the solve
  int const rank = MpiWrapper::commRank( mat.getComm() );
  m_workingRank = MpiWrapper::min( mat.numLocalRows() > 0 ? rank : std::numeric_limits< int >::max(), mat.getComm() );

  SSlong const numGR = LvArray::integerConversion< SSlong >( mat.numGlobalRows() );
  SSlong const numNZ = LvArray::integerConversion< SSlong >( mat.numGlobalNonzeros() );

  // Allocate memory and control structures on working rank only
  if( m_workingRank == rank )
  {
    m_data = std::make_unique< SuiteSparseData >( numGR, numNZ );
    setOptions();
  }

  // Export needs to be carried collectively on all ranks
  m_export = std::make_unique< typename Matrix::Export >( mat, m_workingRank );

  arrayView1d< SSlong > rowPtrView;
  arrayView1d< SSlong > colIndicesView;
  arrayView1d< double > valuesView;
  if( rank == m_workingRank )
  {
    rowPtrView = m_data->rowPtr.toView();
    colIndicesView = m_data->colIndices.toView();
    valuesView = m_data->values.toView();
  }
  else
  {
    array1d< SSlong > dummySSlongArray1d;
    array1d< double > dummyDoubleArray1d;
    rowPtrView = dummySSlongArray1d.toView();
    colIndicesView = dummySSlongArray1d.toView();
    valuesView = dummyDoubleArray1d.toView();
  }

  m_export->exportCRS( mat,
                       rowPtrView,
                       colIndicesView,
                       valuesView );

  if( rank == m_workingRank )
  {
    m_data->rowPtr.move( LvArray::MemorySpace::host, false );
    m_data->colIndices.move( LvArray::MemorySpace::host, false );
    m_data->values.move( LvArray::MemorySpace::host, false );
  }

  // Perform matrix factorization on working rank and sync timer
  {
    Stopwatch timer( m_result.setupTime );
    factorize();
  }
  MpiWrapper::bcast( &m_result.setupTime, 1, m_workingRank, mat.getComm() );
}

template< typename LAI >
void SuiteSparse< LAI >::apply( Vector const & src,
                                Vector & dst ) const
{
  doSolve( src, dst, false );
}

template< typename LAI >
void SuiteSparse< LAI >::applyTranspose( Vector const & src,
                                         Vector & dst ) const
{
  doSolve( src, dst, true );
}

template< typename LAI >
void SuiteSparse< LAI >::clear()
{
  PreconditionerBase< LAI >::clear();
  m_data.reset();
  m_condEst = -1.0;
}

template< typename LAI >
void SuiteSparse< LAI >::setOptions()
{
  // Get the default control parameters
  umfpack_dl_defaults( m_data->control );
  m_data->control[UMFPACK_PRL] = m_params.logLevel > 1 ? 6 : 1;
  m_data->control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
}

template< typename LAI >
void SuiteSparse< LAI >::factorize()
{
  int const rank = MpiWrapper::commRank( matrix().getComm() );

  if( rank == m_workingRank )
  {
    // To be able to use UMFPACK direct solver we need to disable floating point exceptions
    LvArray::system::FloatingPointExceptionGuard guard;

    SSlong status;
    SSlong const numRows = m_data->rowPtr.size() - 1;

    // symbolic factorization
    status = umfpack_dl_symbolic( numRows,
                                  numRows,
                                  m_data->rowPtr.data(),
                                  m_data->colIndices.data(),
                                  m_data->values.data(),
                                  &m_data->symbolic,
                                  m_data->control,
                                  m_data->info );
    if( status < 0 )
    {
      umfpack_dl_report_info( m_data->control, m_data->info );
      umfpack_dl_report_status( m_data->control, status );
      GEOSX_ERROR( "SuiteSparse: umfpack_dl_symbolic failed." );
    }

    // print the symbolic factorization
    if( m_params.logLevel > 1 )
    {
      umfpack_dl_report_symbolic( m_data->symbolic, m_data->control );
    }

    // numeric factorization
    status = umfpack_dl_numeric( m_data->rowPtr.data(),
                                 m_data->colIndices.data(),
                                 m_data->values.data(),
                                 m_data->symbolic,
                                 &m_data->numeric,
                                 m_data->control,
                                 m_data->info );

    if( status < 0 )
    {
      umfpack_dl_report_info( m_data->control, m_data->info );
      umfpack_dl_report_status( m_data->control, status );
      GEOSX_ERROR( "SuiteSparse: umfpack_dl_numeric failed." );
    }

    // print the numeric factorization
    if( m_params.logLevel > 1 )
    {
      umfpack_dl_report_numeric( m_data->symbolic, m_data->control );
    }
  }
}

template< typename LAI >
void SuiteSparse< LAI >::solve( Vector const & rhs,
                                Vector & sol ) const
{
  real64 const bnorm = rhs.norm2();
  if( isZero( bnorm, 0.0 ) )
  {
    sol.zero();
    m_result.numIterations = 0;
    m_result.residualReduction = 0.0;
    m_result.solveTime = 0.0;
    m_result.status = LinearSolverResult::Status::Success;
    return;
  }

  {
    Stopwatch timer( m_result.solveTime );
    apply( rhs, sol );
  }
  MpiWrapper::bcast( &m_result.solveTime, 1, m_workingRank, rhs.getComm() );

  Vector r( rhs );
  matrix().residual( sol, r, r );
  m_result.residualReduction = r.norm2() / bnorm;

  m_result.status = LinearSolverResult::Status::Success;
  m_result.numIterations = 1;

  if( m_params.direct.checkResidual )
  {
    real64 constexpr precTol = 100.0 * std::numeric_limits< real64 >::epsilon();
    real64 condEst = estimateConditionNumberBasic();
    if( m_result.residualReduction > condEst * precTol )
    {
      condEst = estimateConditionNumberAdvanced();
      if( m_result.residualReduction > condEst * precTol )
      {
        if( m_params.logLevel > 0 )
        {
          GEOSX_WARNING( "SuiteSparse: failed to reduce residual below tolerance.\n"
                         "Condition number estimate: " << condEst );
        }
        m_result.status = LinearSolverResult::Status::Breakdown;
      }
    }
  }

  if( m_params.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( "\t\tLinear Solver | " << m_result.status <<
                      " | Iterations: " << m_result.numIterations <<
                      " | Final Rel Res: " << m_result.residualReduction <<
                      " | Setup Time: " << m_result.setupTime << " s" <<
                      " | Solve Time: " << m_result.solveTime << " s" );
  }
}

template< typename LAI >
void SuiteSparse< LAI >::doSolve( Vector const & b, Vector & x, bool transpose ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( b.ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( b.localSize(), x.localSize() );
  GEOSX_LAI_ASSERT_EQ( b.localSize(), matrix().numLocalRows() );

  int const rank = MpiWrapper::commRank( b.getComm() );

  arrayView1d< double > rhsView;
  if( rank == m_workingRank )
  {
    rhsView = m_data->rhs.toView();
  }
  else
  {
    array1d< double > dummyDoubleArray1d;
    rhsView = dummyDoubleArray1d.toView();
  }

  m_export->exportVector( b, rhsView );

  if( rank == m_workingRank )
  {
    m_data->rhs.move( LvArray::MemorySpace::host, false );

    // To be able to use UMFPACK direct solver we need to disable floating point exceptions
    LvArray::system::FloatingPointExceptionGuard guard;

    // Note: UMFPACK expects column-sparse matrix, but we have row-sparse, so we flip the transpose flag
    SSlong const status = umfpack_dl_solve( transpose ? UMFPACK_A : UMFPACK_At,
                                            m_data->rowPtr.data(),
                                            m_data->colIndices.data(),
                                            m_data->values.data(),
                                            m_data->sol.data(),
                                            m_data->rhs.data(),
                                            m_data->numeric,
                                            m_data->control,
                                            m_data->info );

    if( status < 0 )
    {
      umfpack_dl_report_info( m_data->control, m_data->info );
      umfpack_dl_report_status( m_data->control, status );
      GEOSX_ERROR( "SuiteSparse interface: umfpack_dl_solve failed." );
    }
  }

  arrayView1d< double > solView;
  if( rank == m_workingRank )
  {
    solView = m_data->sol.toView();
  }
  else
  {
    array1d< double > dummyDoubleArray1d;
    solView = dummyDoubleArray1d.toView();
  }
  m_export->importVector( solView, x );
}

template< typename LAI >
real64 SuiteSparse< LAI >::estimateConditionNumberBasic() const
{
  GEOSX_LAI_ASSERT( ready() );
  if( m_condEst >= 0.0 )
  {
    return m_condEst; // used cached result, possibly more accurate
  }

  int const rank = MpiWrapper::commRank( matrix().getComm() );
  if( rank == m_workingRank )
  {
    m_condEst = 1.0 / m_data->info[UMFPACK_RCOND];
  }
  MpiWrapper::bcast( &m_condEst, 1, m_workingRank, matrix().getComm() );

  return m_condEst;
}

template< typename LAI >
real64 SuiteSparse< LAI >::estimateConditionNumberAdvanced() const
{
  GEOSX_LAI_ASSERT( ready() );
  localIndex constexpr numIterations = 4;

  NormalOperator< LAI > const normalOperator( matrix() );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue( normalOperator, numIterations );

  InverseNormalOperator< LAI, SuiteSparse > const inverseNormalOperator( matrix(), *this );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue( inverseNormalOperator, numIterations );

  m_condEst = sqrt( lambdaDirect * lambdaInverse );
  return m_condEst;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class SuiteSparse< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class SuiteSparse< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class SuiteSparse< PetscInterface >;
#endif

}
