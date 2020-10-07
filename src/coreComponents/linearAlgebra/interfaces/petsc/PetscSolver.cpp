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
 * @file PetscSolver.cpp
 */

#include "PetscSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/petsc/PetscVector.hpp"
#include "linearAlgebra/interfaces/petsc/PetscMatrix.hpp"
#include "linearAlgebra/interfaces/petsc/PetscPreconditioner.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include "PetscSuperlu.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#ifdef GEOSX_USE_SUITESPARSE
#include "PetscSuiteSparse.hpp"
#endif

// Put everything under the geosx namespace.
namespace geosx
{

PetscSolver::PetscSolver( LinearSolverParameters parameters )
  :
  m_parameters( std::move( parameters ) )
{}

void PetscSolver::solve( PetscMatrix & mat,
                         PetscVector & sol,
                         PetscVector & rhs,
                         DofManager const * const dofManager )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  GEOSX_UNUSED_VAR( dofManager );

  if( m_parameters.solverType == LinearSolverParameters::SolverType::direct )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs );
  }
}

namespace
{

void solve_parallelDirect( LinearSolverParameters const & parameters,
                           PetscMatrix & mat,
                           PetscVector & sol,
                           PetscVector & rhs,
                           LinearSolverResult & result )
{
  // To be able to use SuperLU_Dist solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;
  GEOSX_LAI_CHECK_ERROR( KSPCreate( comm, &ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, mat.unwrapped(), mat.unwrapped() ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPPREONLY ) );

  SuperLU_DistSetFromOptions( mat, parameters );

  // use direct solve preconditioner SUPERLU DIST
  Stopwatch watch;
  PC prec;
  GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp, &prec ) );
  GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCLU ) );
  GEOSX_LAI_CHECK_ERROR( PCFactorSetMatSolverType( prec, MATSOLVERSUPERLU_DIST ) );
  GEOSX_LAI_CHECK_ERROR( PCSetUp( prec ) );
  result.setupTime = watch.elapsedTime();

  // solve system
  watch.zero();
  GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, rhs.unwrapped(), sol.unwrapped() ) );
  result.solveTime = watch.elapsedTime();

  KSPConvergedReason reason;
  GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( ksp, &reason ) );

  result.status = reason >= 0 ? LinearSolverResult::Status::Success : LinearSolverResult::Status::Breakdown;

  if( result.status == LinearSolverResult::Status::Success )
  {
    PetscVector res( rhs );
    mat.gemv( -1.0, sol, 1.0, res );
    result.residualReduction = res.norm2() / rhs.norm2();

    // check for nan or inf
    if( std::isnan( result.residualReduction ) || std::isinf( result.residualReduction ) )
    {
      result.status = LinearSolverResult::Status::Breakdown;
    }
    else if( result.residualReduction < parameters.direct.checkResidualTolerance )
    {
      result.status = LinearSolverResult::Status::Success;
      result.numIterations = 1;
    }
    else
    {
      result.status = LinearSolverResult::Status::Breakdown;
    }
  }

  // destroy solver
  GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );
}

#ifdef GEOSX_USE_SUITESPARSE
void solve_serialDirect( LinearSolverParameters const & parameters,
                         PetscMatrix & mat,
                         PetscVector & sol,
                         PetscVector & rhs,
                         LinearSolverResult & result )
{
  // To be able to use UMFPACK direct solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  MPI_Comm const comm = mat.getComm();

  int const rank = MpiWrapper::Comm_rank( comm );

  SuiteSparseSetFromOptions( mat, parameters );

  Mat * localMatrix;
  IS set;
  if( rank == 0 )
  {
    GEOSX_LAI_CHECK_ERROR( ISCreateStride( PETSC_COMM_SELF, mat.numGlobalRows(), 0, 1, &set ) );
  }
  GEOSX_LAI_CHECK_ERROR( MatCreateSubMatrices( mat.unwrapped(), rank==0, &set, &set, MAT_INITIAL_MATRIX, &localMatrix ) );

  VecScatter rhsScatter;
  Vec localRhs;
  GEOSX_LAI_CHECK_ERROR( VecScatterCreateToZero( rhs.unwrapped(), &rhsScatter, &localRhs ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterBegin( rhsScatter, rhs.unwrapped(), localRhs, INSERT_VALUES, SCATTER_FORWARD ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterEnd( rhsScatter, rhs.unwrapped(), localRhs, INSERT_VALUES, SCATTER_FORWARD ) );

  VecScatter solScatter;
  Vec localSol;
  GEOSX_LAI_CHECK_ERROR( VecScatterCreateToZero( sol.unwrapped(), &solScatter, &localSol ) );

  localIndex status = 0;

  if( rank == 0 )
  {
    // create linear solver
    KSP ksp;
    GEOSX_LAI_CHECK_ERROR( KSPCreate( PETSC_COMM_SELF, &ksp ) );
    GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, localMatrix[0], localMatrix[0] ) );
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPPREONLY ) );

    // use direct solve preconditioner UMFPACK
    Stopwatch watch;
    PC prec;
    GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp, &prec ) );
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCLU ) );
    GEOSX_LAI_CHECK_ERROR( PCFactorSetMatSolverType( prec, MATSOLVERUMFPACK ) );
    GEOSX_LAI_CHECK_ERROR( PCSetUp( prec ) );
    result.setupTime = watch.elapsedTime();

    // solve system
    watch.zero();
    GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, localRhs, localSol ) );
    result.solveTime = watch.elapsedTime();

    KSPConvergedReason reason;
    GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( ksp, &reason ) );

    // save status
    status = !( reason >= 0 );

    // destroy solver
    GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );
  }

  // broadcast status and times
  MpiWrapper::bcast( &status, 1, 0, comm );
  MpiWrapper::bcast( &result.setupTime, 1, 0, comm );
  MpiWrapper::bcast( &result.solveTime, 1, 0, comm );

  GEOSX_LAI_CHECK_ERROR( VecScatterBegin( solScatter, localSol, sol.unwrapped(), INSERT_VALUES, SCATTER_REVERSE ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterEnd( solScatter, localSol, sol.unwrapped(), INSERT_VALUES, SCATTER_REVERSE ) );
  GEOSX_LAI_CHECK_ERROR( VecDestroy( &localRhs ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterDestroy( &rhsScatter ) );
  GEOSX_LAI_CHECK_ERROR( VecDestroy( &localSol ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterDestroy( &solScatter ) );
  GEOSX_LAI_CHECK_ERROR( MatDestroySubMatrices( rank==0, &localMatrix ) );

  result.status = status == 0 ? LinearSolverResult::Status::Success : LinearSolverResult::Status::Breakdown;

  if( result.status == LinearSolverResult::Status::Success )
  {
    PetscVector res( rhs );
    mat.gemv( -1.0, sol, 1.0, res );
    result.residualReduction = res.norm2() / rhs.norm2();

    // check for nan or inf
    if( std::isnan( result.residualReduction ) || std::isinf( result.residualReduction ) )
    {
      result.status = LinearSolverResult::Status::Breakdown;
    }
    else if( result.residualReduction < parameters.direct.checkResidualTolerance )
    {
      result.status = LinearSolverResult::Status::Success;
      result.numIterations = 1;
    }
    else
    {
      result.status = LinearSolverResult::Status::Breakdown;
    }
  }
}
#endif

void CreatePetscKrylovSolver( LinearSolverParameters const & params,
                              MPI_Comm const comm,
                              KSP & ksp )
{
  GEOSX_LAI_CHECK_ERROR( KSPCreate( comm, &ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetNormType( ksp, KSP_NORM_UNPRECONDITIONED ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetTolerances( ksp, params.krylov.relTolerance, PETSC_DEFAULT,
                                           PETSC_DEFAULT, params.krylov.maxIterations ) );

  // pick the solver type
  switch( params.solverType )
  {
    case LinearSolverParameters::SolverType::gmres:
    {
      GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPGMRES ) );
      GEOSX_LAI_CHECK_ERROR( KSPGMRESSetRestart( ksp, params.krylov.maxRestart ) );
      break;
    }
    case LinearSolverParameters::SolverType::bicgstab:
    {
      GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPBCGS ) );
      break;
    }
    case LinearSolverParameters::SolverType::cg:
    {
      GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPCG ) );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Solver type not supported in PETSc interface: " << params.solverType );
    }
  }
}

} // namespace

void PetscSolver::solve_direct( PetscMatrix & mat,
                                PetscVector & sol,
                                PetscVector & rhs )
{
  if( m_parameters.direct.parallel )
  {
    solve_parallelDirect( m_parameters, mat, sol, rhs, m_result );
  }
  else
  {
#ifdef GEOSX_USE_SUITESPARSE
    solve_serialDirect( m_parameters, mat, sol, rhs, m_result );
#else
    GEOSX_ERROR( "Petsc direct solver interface: serial direct solver not available (try to compile GEOSX TPLs with SuiteSparse)." );
#endif
  }
}

void PetscSolver::solve_krylov( PetscMatrix & mat,
                                PetscVector & sol,
                                PetscVector & rhs )
{
  Stopwatch watch;

  // create linear solver
  KSP ksp;
  CreatePetscKrylovSolver( m_parameters, mat.getComm(), ksp );

  // This can be used to extract residual norm history:
  //array1d< real64 > residualNorms( m_parameters.krylov.maxIterations + 1 );
  //GEOSX_LAI_CHECK_ERROR( KSPSetResidualHistory( ksp, residualNorms.data(), residualNorms.size(), PETSC_TRUE ) );

  // Deal with separate component approximation
  PetscMatrix separateComponentMatrix;
  if( m_parameters.amg.separateComponents )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat, separateComponentMatrix, m_parameters.dofsPerNode );
  }
  PetscMatrix & precondMat = m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // create and compute a preconditioner and set into KSP
  PetscPreconditioner precond( m_parameters );
  precond.compute( precondMat );
  KSPSetPC( ksp, precond.unwrapped() );

  GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, mat.unwrapped(), precondMat.unwrapped() ) );

  // display output
  if( m_parameters.logLevel > 0 )
  {
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-ksp_monitor", nullptr ) );
  }
  GEOSX_LAI_CHECK_ERROR( KSPSetFromOptions( ksp ) );

  m_result.setupTime = watch.elapsedTime();

  // Actually solve
  watch.zero();
  GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, rhs.unwrapped(), sol.unwrapped() ) );
  m_result.solveTime = watch.elapsedTime();

  // Get status indicator
  KSPConvergedReason result;
  GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( ksp, &result ) );
  GEOSX_WARNING_IF( result < 0, "PetscSolver: Krylov convergence not achieved" );

  switch( result )
  {
    case KSP_CONVERGED_RTOL:
    case KSP_CONVERGED_ATOL:
    case KSP_CONVERGED_ITS:
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }
    case KSP_DIVERGED_ITS:
    case KSP_DIVERGED_DTOL:
    {
      m_result.status = LinearSolverResult::Status::NotConverged;
      break;
    }
    default:
    {
      m_result.status = LinearSolverResult::Status::Breakdown;
    }
  }

  // Get number of iterations performed
  PetscInt numIter;
  GEOSX_LAI_CHECK_ERROR( KSPGetIterationNumber( ksp, &numIter ) );
  m_result.numIterations = numIter;

  // reset verbosity option
  GEOSX_LAI_CHECK_ERROR( PetscOptionsClearValue( nullptr, "-ksp_monitor" ) );

  GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );
}

} // end geosx namespace
