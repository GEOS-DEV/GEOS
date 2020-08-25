/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DirectSolver.cpp
 */

#include "common/DataTypes.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "DirectSolver.hpp"
#include "DirectUtils.hpp"

#ifdef GEOSX_USE_SUPERLU_DIST
#include <superlu_ddefs.h>
#endif

namespace geosx
{

template< typename LAI >
DirectSolver< LAI >::DirectSolver( LinearSolverParameters const params )
  : m_params( params )
{}

#ifdef GEOSX_USE_SUPERLU_DIST

template< typename LAI >
void ConvertToSuperMatrix( typename LAI::ParallelMatrix const & matrix,
                           array1d< globalIndex > & rowPtr,
                           array1d< globalIndex > & cols,
                           array1d< real64 > & vals,
                           SuperMatrix & SLUDMat )
{
  localIndex const numLocalRows = matrix.numLocalRows();
  rowPtr.resize( numLocalRows+1 );
  cols.reserve( matrix.numLocalNonzeros() );
  vals.reserve( matrix.numLocalNonzeros() );
  rowPtr[0] = 0;
  for( localIndex i = 0; i < numLocalRows; ++i )
  {
    localIndex const nonZeros = matrix.localRowLength( i );
    array1d< globalIndex > colIndices( nonZeros );
    array1d< real64 > values( nonZeros );
    matrix.getRowCopy( matrix.getGlobalRowID( i ), colIndices, values );
    cols.insert( rowPtr[i], colIndices.begin(), colIndices.end() );
    vals.insert( rowPtr[i], values.begin(), values.end() );
    rowPtr[i+1] = rowPtr[i] + nonZeros;
  }

  dCreate_CompRowLoc_Matrix_dist( &SLUDMat,
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numLocalNonzeros() ),
                                  toSuperlu_intT( numLocalRows ),
                                  toSuperlu_intT( matrix.ilower() ),
                                  vals.data(),
                                  toSuperlu_intT( cols.data() ),
                                  toSuperlu_intT( rowPtr.data() ),
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );
}

template< typename LAI >
int SolveSuperMatrix( SuperMatrix & SLUDMat,
                      typename LAI::ParallelVector const & b,
                      typename LAI::ParallelVector & x,
                      MPI_Comm const & comm,
                      superlu_dist_options_t & options,
                      integer const & logLevel,
                      real64 & timeFact,
                      real64 & timeSolve )
{
  Stopwatch watch;

  int_t const m = SLUDMat.nrow;
  int_t const n = SLUDMat.ncol;

  // Initialize ScalePermstruct.
  dScalePermstruct_t ScalePermstruct;
  dScalePermstructInit( m, n, &ScalePermstruct );

  // Initialize LUstruct.
  dLUstruct_t LUstruct;
  dLUstructInit( n, &LUstruct );

  // Initialize the statistics variables.
  SuperLUStat_t stat;
  PStatInit( &stat );

  // Create process grid.
  int const num_procs = MpiWrapper::Comm_size( comm );
  int pcols = 1;
  int prows = 1;
  while( prows*pcols <= num_procs )
  {
    ++prows;
  }
  --prows;
  pcols = num_procs/prows;
  while( prows*pcols != num_procs )
  {
    prows -= 1;
    pcols = num_procs/prows;
  }
  gridinfo_t grid;
  superlu_gridinit( comm, prows, pcols, &grid );

  // Call the linear equation solver.
  int nrhs = 0;
  int const ldb = b.localSize();
  dSOLVEstruct_t SOLVEstruct;
  array1d< real64 > berr( b.localSize() );
  int info = 0;

  options.Fact = DOFACT;
  pdgssvx( &options,
           &SLUDMat,
           &ScalePermstruct,
           NULL,
           ldb,
           nrhs,
           &grid,
           &LUstruct,
           &SOLVEstruct,
           berr.data(),
           &stat,
           &info );

  timeFact = watch.elapsedTime();
  watch.zero();

  if( info == 0 )
  {
    options.Fact = FACTORED;
    nrhs = 1;
    x.copy( b );
    pdgssvx( &options,
             &SLUDMat,
             &ScalePermstruct,
             x.extractLocalVector(),
             ldb,
             nrhs,
             &grid,
             &LUstruct,
             &SOLVEstruct,
             berr.data(),
             &stat,
             &info );
  }

  timeSolve = watch.elapsedTime();

  if( logLevel > 0 )
  {
    // Print the statistics.
    PStatPrint( &options, &stat, &grid );
  }

  // Deallocate other SuperLU data structures
  dScalePermstructFree( &ScalePermstruct );
  dDestroy_LU( n, &grid, &LUstruct );
  dLUstructFree( &LUstruct );
  if( options.SolveInitialized )
  {
    dSolveFinalize( &options, &SOLVEstruct );
  }

  return info;
}

template< typename LAI >
void DirectSolver< LAI >::solve( Matrix & matrix,
                                 Vector & b,
                                 Vector & x ) const
{
  GEOSX_LOG_RANK_0( "===========> Using direct solver from SuperLU_Dist package. <===========" );
  // Convert matrix from Matrix to SuperMatrix format
  array1d< globalIndex > rowPtr;
  array1d< globalIndex > cols;
  array1d< real64 > vals;
  SuperMatrix SLUDMat;
  ConvertToSuperMatrix< LAI >( matrix, rowPtr, cols, vals, SLUDMat );

  MPI_Comm const comm = matrix.getComm();

  // Initialize options.
  superlu_dist_options_t options;
  set_default_options_dist( &options );
  options.ReplaceTinyPivot = YES;
  if( m_params.logLevel > 1 )
  {
    options.PrintStat = YES;
  }
  else
  {
    options.PrintStat = NO;
  }

  if( m_params.logLevel > 0 )
  {
    print_sp_ienv_dist( &options );
    print_options_dist( &options );
  }

  real64 const bnorm2 = b.norm2();

  real64 timeFact, timeSolve;
  int const info = SolveSuperMatrix< LAI >( SLUDMat, b, x, comm, options, m_params.logLevel, timeFact, timeSolve );

  Vector res( b );
  matrix.gemv( -1.0, x, 1.0, res );
  m_result.residualReduction = res.norm2() / bnorm2;

  if( info == 0 && m_result.residualReduction < m_params.krylov.relTolerance )
  {
    m_result.status = LinearSolverResult::Status::Success;
    m_result.setupTime = timeFact;
    m_result.solveTime = timeSolve;
  }
  else
  {
    m_result.status = LinearSolverResult::Status::NotConverged;
  }

  // Second and last trial
  if( m_result.status == LinearSolverResult::Status::NotConverged )
  {
    options.ParSymbFact = YES;
    options.ColPerm = PARMETIS;
    if( m_params.logLevel > 0 )
    {
      GEOSX_LOG_RANK_0( "Second attempt with different options." );
      print_sp_ienv_dist( &options );
      print_options_dist( &options );
    }

    real64 timeFact2, timeSolve2;
    int const info2 = SolveSuperMatrix< LAI >( SLUDMat, b, x, comm, options, m_params.logLevel, timeFact2, timeSolve2 );

    res.copy( b );
    matrix.gemv( -1.0, x, 1.0, res );
    m_result.residualReduction = res.norm2() / bnorm2;

    if( info2 == 0 && m_result.residualReduction < m_params.krylov.relTolerance )
    {
      m_result.status = LinearSolverResult::Status::Success;
      m_result.setupTime = timeFact + timeFact2;
      m_result.solveTime = timeSolve + timeSolve2;
    }
    else
    {
      m_result.status = LinearSolverResult::Status::NotConverged;
    }
  }

  // rowPtr, cols, vals will be deallocated when out of scope
  SUPERLU_FREE( SLUDMat.Store );
}

#else

template< typename LAI >
void DirectSolver< LAI >::solve( Matrix & matrix,
                                 Vector & b,
                                 Vector & x ) const
{
  GEOSX_LOG_RANK_0( "===========> Using direct solver from specific LAI. <===========" );
  typename LAI::LinearSolver solver( m_params );
  solver.solve( matrix, x, b );
  m_result = solver.result();
}

#endif

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class DirectSolver< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class DirectSolver< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class DirectSolver< PetscInterface >;
#endif

}
