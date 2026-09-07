/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreSolver.cpp
 */

#include "HypreSolver.hpp"

#include "common/MpiWrapper.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_IJ_mv.h>
#include <HYPRE_krylov.h>

#include <algorithm>

namespace geos
{

/**
 * @brief Container for hypre Krylov solver function pointers.
 */
struct HypreSolverWrapper : public HyprePrecWrapper
{
  /// Alias for set preconditioner function type
  using SetPrecondFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_Solver );

  /// Alias for get number of iterations function type
  using GetNumIter = HYPRE_Int ( * )( HYPRE_Solver solver,
                                      HYPRE_Int * num_iterations );

  /// Alias for get final residual norm function type
  using GetFinalNorm = HYPRE_Int ( * )( HYPRE_Solver solver,
                                        HYPRE_Real * norm );

  SetPrecondFunc setPrecond{}; ///< pointer to set preconditioner function
  GetNumIter getNumIter{};     ///< pointer to get number of iterations function
  GetFinalNorm getFinalNorm{}; ///< pointer to get final residual norm function

  /// DoF-component labels used by hypre's tagged Krylov reductions.
  array1d< HYPRE_Int > krylovDofTags;
  HYPRE_Int numKrylovDofTags = 1;
  HypreVector dummy;
};

HypreSolver::HypreSolver( LinearSolverParameters parameters )
  : Base( std::move( parameters ) ),
  m_precond( m_params )
{ }

HypreSolver::~HypreSolver()
{
  HypreSolver::clear();
}

namespace
{

void createHypreGMRES( LinearSolverParameters const & params,
                       MPI_Comm const comm,
                       HypreSolverWrapper & solver )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESCreate( comm, &solver.ptr ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetKDim( solver.ptr, params.krylov.maxRestart ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  HYPRE_Int logLevel = (params.logLevel >= 3) ? 3 : 0;

  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrintLevel( solver.ptr, logLevel ) ); // print iteration info
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetLogging( solver.ptr, 1 ) ); /* needed to get run info later */

  solver.setPrecond = HYPRE_ParCSRGMRESSetPrecond;
  solver.setup = HYPRE_ParCSRGMRESSetup;
  solver.solve = HYPRE_ParCSRGMRESSolve;
  solver.getNumIter = HYPRE_GMRESGetNumIterations;
  solver.getFinalNorm = HYPRE_GMRESGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRGMRESDestroy;
}

void createHypreFlexGMRES( LinearSolverParameters const & params,
                           MPI_Comm const comm,
                           HypreSolverWrapper & solver )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESCreate( comm, &solver.ptr ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetKDim( solver.ptr, params.krylov.maxRestart ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  HYPRE_Int logLevel = (params.logLevel >= 3) ? 2 : 0;

  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetPrintLevel( solver.ptr, logLevel ) ); // print iteration info
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetLogging( solver.ptr, 1 ) ); /* needed to get run info later */

  solver.setPrecond = HYPRE_ParCSRFlexGMRESSetPrecond;
  solver.setup = HYPRE_ParCSRFlexGMRESSetup;
  solver.solve = HYPRE_ParCSRFlexGMRESSolve;
  solver.getNumIter = HYPRE_FlexGMRESGetNumIterations;
  solver.getFinalNorm = HYPRE_FlexGMRESGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRFlexGMRESDestroy;
}

void createHypreBiCGSTAB( LinearSolverParameters const & params,
                          MPI_Comm const comm,
                          HypreSolverWrapper & solver )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABCreate( comm, &solver.ptr ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  HYPRE_Int logLevel = (params.logLevel >= 3) ? 1 : 0;

  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetPrintLevel( solver.ptr, logLevel ) ); // print iteration info
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetLogging( solver.ptr, 1 ) ); // needed to get run info later

  solver.setPrecond = HYPRE_ParCSRBiCGSTABSetPrecond;
  solver.setup = HYPRE_ParCSRBiCGSTABSetup;
  solver.solve = HYPRE_ParCSRBiCGSTABSolve;
  solver.getNumIter = HYPRE_BiCGSTABGetNumIterations;
  solver.getFinalNorm = HYPRE_BiCGSTABGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRBiCGSTABDestroy;
}

void createHypreCG( LinearSolverParameters const & params,
                    MPI_Comm const comm,
                    HypreSolverWrapper & solver )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_ParCSRPCGCreate( comm, &solver.ptr ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  HYPRE_Int logLevel = (params.logLevel >= 3) ? 1 : 0;

  GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetPrintLevel( solver.ptr, logLevel ) ); /* print the iteration info */
  GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetLogging( solver.ptr, 1 ) );    /* needed to get run info later */
  GEOS_LAI_CHECK_ERROR( HYPRE_PCGSetTwoNorm( solver.ptr, 1 ) );    /* use the two norm as the stopping criteria */

  solver.setPrecond = HYPRE_ParCSRPCGSetPrecond;
  solver.setup = HYPRE_ParCSRPCGSetup;
  solver.solve = HYPRE_ParCSRPCGSolve;
  solver.getNumIter = HYPRE_PCGGetNumIterations;
  solver.getFinalNorm = HYPRE_PCGGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRPCGDestroy;
}

void createHypreKrylovSolver( LinearSolverParameters const & params,
                              MPI_Comm const comm,
                              HypreSolverWrapper & solver )
{
  switch( params.solverType )
  {
    case LinearSolverParameters::SolverType::gmres:
    {
      createHypreGMRES( params, comm, solver );
      break;
    }
    case LinearSolverParameters::SolverType::fgmres:
    {
      createHypreFlexGMRES( params, comm, solver );
      break;
    }
    case LinearSolverParameters::SolverType::bicgstab:
    {
      createHypreBiCGSTAB( params, comm, solver );
      break;
    }
    case LinearSolverParameters::SolverType::cg:
    {
      createHypreCG( params, comm, solver );
      break;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Solver type not supported in hypre interface: {}", params.solverType ) );
    }
  }
}

#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CUDA && GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
bool hasMatchingKrylovDofTags( HYPRE_IJVector const vector,
                               HYPRE_Int const numTags,
                               HYPRE_Int const * const tags )
{
  if( vector == nullptr || hypre_IJVectorObjectType( vector ) != HYPRE_PARCSR )
  {
    return false;
  }

  hypre_ParVector * const parVector = static_cast< hypre_ParVector * >( hypre_IJVectorObject( vector ) );
  if( parVector == nullptr )
  {
    return false;
  }

  hypre_Vector * const localVector = hypre_ParVectorLocalVector( parVector );
  if( localVector == nullptr || hypre_VectorNumTags( localVector ) != numTags )
  {
    return false;
  }

  HYPRE_Int const localSize = hypre_VectorSize( localVector );
  return localSize == 0 || (tags != nullptr &&
                            std::equal( tags, tags + localSize, hypre_VectorTags( localVector ) ));
}

void tagKrylovDofVector( HypreSolverWrapper & solver,
                         HypreVector const & vec )
{
  if( solver.numKrylovDofTags <= 1 )
  {
    return;
  }

  HYPRE_IJVector const ijVector = vec.unwrappedIJ();
  bool const tagsChanged = !hasMatchingKrylovDofTags( ijVector,
                                                      solver.numKrylovDofTags,
                                                      solver.krylovDofTags.data() );
  if( tagsChanged )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorSetTags( ijVector, 1, solver.numKrylovDofTags,
                                                 solver.krylovDofTags.data() ) );
  }
  if( MpiWrapper::max( tagsChanged ? 1 : 0, vec.comm() ) > 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( ijVector ) );
  }
}
#endif

} // namespace

void HypreSolver::setup( HypreMatrix const & mat )
{
  GEOS_MARK_FUNCTION;

  clear();
  Base::setup( mat );
  Stopwatch timer( m_result.setupTime );
  m_precond.setup( mat );

  m_solver = std::make_unique< HypreSolverWrapper >();
  createHypreKrylovSolver( m_params, mat.comm(), *m_solver );

  // Set the preconditioner
  GEOS_LAI_CHECK_ERROR( m_solver->setPrecond( m_solver->ptr,
                                              m_precond.unwrapped().solve,
                                              hypre::dummySetup,
                                              m_precond.unwrapped().ptr ) );

  // Setup the solver (need a dummy vector for rhs/sol to avoid hypre segfaulting in setup).
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CUDA && GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
  array1d< int > labels;
  hypre::fillKrylovDofLabels( mat, labels );
  hypre::assignKrylovDofTags( labels.toViewConst(),
                              mat.comm(),
                              m_solver->krylovDofTags,
                              m_solver->numKrylovDofTags );
#endif
  // hypre's Krylov work vectors borrow the setup vector's tag array, so keep
  // the tagged dummy alive for the lifetime of the solver.
  m_solver->dummy.create( mat.numLocalRows(), mat.comm() );
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CUDA && GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
  if( m_solver->numKrylovDofTags > 1 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorSetTags( m_solver->dummy.unwrappedIJ(),
                                                 0,
                                                 m_solver->numKrylovDofTags,
                                                 m_solver->krylovDofTags.data() ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( m_solver->dummy.unwrappedIJ() ) );
  }
#endif
  GEOS_LAI_CHECK_ERROR( m_solver->setup( m_solver->ptr,
                                         mat.unwrapped(),
                                         m_solver->dummy.unwrapped(),
                                         m_solver->dummy.unwrapped() ) );
}

int HypreSolver::doSolve( HypreVector const & rhs,
                          HypreVector & sol ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( sol.ready() );
  GEOS_LAI_ASSERT( rhs.ready() );
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CUDA && GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
  tagKrylovDofVector( *m_solver, rhs );
  tagKrylovDofVector( *m_solver, sol );
#endif
  HYPRE_Int const result = m_solver->solve( m_solver->ptr, matrix().unwrapped(), rhs.unwrapped(), sol.unwrapped() );
  sol.touch();
  return result;
}

void HypreSolver::apply( HypreVector const & rhs,
                         HypreVector & sol ) const
{
  doSolve( rhs, sol );
}

void HypreSolver::solve( HypreVector const & rhs,
                         HypreVector & sol ) const
{
  GEOS_MARK_FUNCTION;

  if( isZero( rhs.norm2(), 0.0 ) )
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
    int const result = doSolve( rhs, sol );
    m_result.status = result ? LinearSolverResult::Status::NotConverged : LinearSolverResult::Status::Success;
  }

  // Clear error code to avoid GEOS from crashing if Krylov method did not converge
  GEOS_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

  // Get final residual norm
  GEOS_LAI_CHECK_ERROR( m_solver->getFinalNorm( m_solver->ptr, &m_result.residualReduction ) );

  // Get number of iterations
  HYPRE_Int numIter;
  GEOS_LAI_CHECK_ERROR( m_solver->getNumIter( m_solver->ptr, &numIter ) );
  m_result.numIterations = numIter;

  if( m_params.logLevel >= 1 )
  {
    HYPRE_BigInt global_num_rows, global_num_nonzeros;

    // This involves an MPI collective call, and therefore we call it only when necessary
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixGetGlobalInfo( matrix().unwrappedIJ(),
                                                       &global_num_rows,
                                                       &global_num_rows, // This is intentional and assuming the matrix is square
                                                       &global_num_nonzeros ) );

    GEOS_LOG_RANK_0( GEOS_FMT( "        Linear Solver | {} | Unknowns: {} | Nonzeros: {} | Iterations: {} | Final Rel Res: {:.4e} | Setup Time: {:.3f} s | Solve Time: {:.3f} s",
                               m_result.status, stringutilities::addCommaSeparators( global_num_rows ),
                               stringutilities::addCommaSeparators( global_num_nonzeros ), m_result.numIterations,
                               m_result.residualReduction, m_result.setupTime, m_result.solveTime ) );
  }
}

void HypreSolver::clear()
{
  Base::clear();
  m_precond.clear();

  if( m_solver )
  {
    GEOS_LAI_CHECK_ERROR( m_solver->destroy( m_solver->ptr ) );
  }
  m_solver.reset();
}

} // end geos namespace
