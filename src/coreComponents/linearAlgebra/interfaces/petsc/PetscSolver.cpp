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

#include "PetscSuperLU_Dist.hpp"
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

  if( rhs.norm2() > 0.0 )
  {
    if( m_parameters.solverType == LinearSolverParameters::SolverType::direct )
    {
      solve_direct( mat, sol, rhs );
    }
    else
    {
      solve_krylov( mat, sol, rhs );
    }
  }
  else
  {
    sol.zero();
    m_result.status = LinearSolverResult::Status::Success;
    m_result.setupTime = 0.0;
    m_result.solveTime = 0.0;
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

  SuperLU_Dist SLUDData( parameters );
  Mat localMatrix;
  PetscConvertToSuperMatrix( mat, localMatrix, SLUDData );

  GEOSX_LAI_CHECK_ERROR( SLUDData.setup() );
  GEOSX_LAI_CHECK_ERROR( SLUDData.solve( rhs.extractLocalVector(), sol.extractLocalVector() ) );

  // Save setup and solution times
  result.setupTime = SLUDData.setupTime();
  result.solveTime = SLUDData.solveTime();

  PetscVector res( rhs );
  mat.gemv( -1.0, sol, 1.0, res );
  result.residualReduction = res.norm2() / rhs.norm2();

  result.status = parameters.direct.checkResidual == 0 ? LinearSolverResult::Status::Success : LinearSolverResult::Status::Breakdown;
  result.numIterations = 1;
  if( parameters.direct.checkResidual )
  {
    if( result.residualReduction < SLUDData.relativeTolerance() )
    {
      result.status = LinearSolverResult::Status::Success;
    }
    else
    {
      real64 const cond = PetscSuperLU_DistCond( mat, SLUDData );
      if( parameters.logLevel > 0 )
      {
        GEOSX_LOG_RANK_0( "Using a more accurate estimate of the number of conditions" );
        GEOSX_LOG_RANK_0( "Condition number is " << cond );
      }
      if( result.residualReduction < SLUDData.machinePrecision() * cond )
      {
        result.status = LinearSolverResult::Status::Success;
      }
    }
  }

  PetscDestroyAdditionalData( localMatrix );
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

  SuiteSparse SSData( parameters );
  ConvertPetscToSuiteSparseMatrix( mat, SSData );

  GEOSX_LAI_CHECK_ERROR( SSData.setup() );
  GEOSX_LAI_CHECK_ERROR( SuiteSparseSolve( SSData, rhs, sol ) );

  // Save setup and solution times
  result.setupTime = SSData.setupTime();
  result.solveTime = SSData.solveTime();

  PetscVector res( rhs );
  mat.gemv( -1.0, sol, 1.0, res );
  result.residualReduction = res.norm2() / rhs.norm2();

  result.status = parameters.direct.checkResidual == 0 ? LinearSolverResult::Status::Success : LinearSolverResult::Status::Breakdown;
  result.numIterations = 1;
  if( parameters.direct.checkResidual )
  {
    if( result.residualReduction < SSData.relativeTolerance() )
    {
      result.status = LinearSolverResult::Status::Success;
    }
    else
    {
      real64 const cond = PetscSuiteSparseCond( mat, SSData );
      if( parameters.logLevel > 0 )
      {
        GEOSX_LOG_RANK_0( "Using a more accurate estimate of the number of conditions" );
        GEOSX_LOG_RANK_0( "Condition number is " << cond );
      }
      if( result.residualReduction < SSData.machinePrecision() * cond )
      {
        result.status = LinearSolverResult::Status::Success;
      }
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
