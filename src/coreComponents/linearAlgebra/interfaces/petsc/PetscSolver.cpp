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
 * @file PetscSolver.cpp
 */

#include "PetscSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <petscvec.h>
#include <petscksp.h>

// Put everything under the geosx namespace.
namespace geosx
{

PetscSolver::PetscSolver( LinearSolverParameters parameters )
  : Base( std::move( parameters ) ),
  m_precond( m_params )
{}

PetscSolver::~PetscSolver()
{
  PetscSolver::clear();
}

namespace
{

void createPetscKrylovSolver( LinearSolverParameters const & params,
                              MPI_Comm const & comm,
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

void PetscSolver::setup( PetscMatrix const & mat )
{
  clear();
  Base::setup( mat );
  Stopwatch timer( m_result.setupTime );

  m_precond.setup( mat );

  createPetscKrylovSolver( m_params, mat.getComm(), m_solver );
  GEOSX_LAI_CHECK_ERROR( KSPSetPC( m_solver, m_precond.unwrapped() ) );

  // display output
  if( m_params.logLevel >= 1 )
  {
    // cast needed because of "void *" vs "PetscViewerAndFormat *" in last parameter
    using MonitorFunc = PetscErrorCode ( * )( KSP, PetscInt, PetscReal, void * );
    GEOSX_LAI_CHECK_ERROR( KSPMonitorSet( m_solver, ( MonitorFunc ) KSPMonitorDefault, nullptr, nullptr ) );
  }

  // This can be used to extract residual norm history:
  //GEOSX_LAI_CHECK_ERROR( KSPSetResidualHistory( ksp, residualNorms.data(), residualNorms.size(), PETSC_TRUE ) );
}

void PetscSolver::apply( PetscVector const & rhs,
                         PetscVector & sol ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );
  GEOSX_LAI_CHECK_ERROR( KSPSolve( m_solver, rhs.unwrapped(), sol.unwrapped() ) );
}

void PetscSolver::solve( PetscVector const & rhs,
                         PetscVector & sol ) const
{
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
    Stopwatch watch( m_result.solveTime );
    apply( rhs, sol );
  }

  // Get status indicator
  KSPConvergedReason result;
  GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( m_solver, &result ) );
  m_result.status = result > 0
                    ? LinearSolverResult::Status::Success
                    : ( result == KSP_DIVERGED_ITS || result == KSP_DIVERGED_DTOL )
                      ? LinearSolverResult::Status::NotConverged
                      : LinearSolverResult::Status::Breakdown;

  // Get number of iterations performed
  PetscInt numIter;
  GEOSX_LAI_CHECK_ERROR( KSPGetIterationNumber( m_solver, &numIter ) );
  m_result.numIterations = numIter;

  GEOSX_LAI_CHECK_ERROR( KSPGetResidualNorm( m_solver, &m_result.residualReduction ) );
  m_result.residualReduction /= rhs.norm2(); // this assumes initial sol is zero

  if( m_params.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( "\t\tLinear Solver | " << m_result.status <<
                      " | Iterations: " << m_result.numIterations <<
                      " | Final Rel Res: " << m_result.residualReduction <<
                      " | Setup Time: " << m_result.setupTime << " s" <<
                      " | Solve Time: " << m_result.solveTime << " s" );
  }
}

void PetscSolver::clear()
{
  PreconditionerBase::clear();
  if( m_solver )
  {
    GEOSX_LAI_CHECK_ERROR( KSPDestroy( &m_solver ) );
  }
}

} // end geosx namespace
