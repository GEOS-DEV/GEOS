/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreSolver.cpp
 */

#include "HypreSolver.hpp"

#include "HypreMatrix.hpp"
#include "HypreVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include "_hypre_utilities.h"
#include "_hypre_parcsr_ls.h"
#include "_hypre_IJ_mv.h"
#include "krylov.h"

namespace geosx
{

typedef HYPRE_Int (* HYPRE_PtrToSolverDestroyFcn)( HYPRE_Solver );

HypreSolver::HypreSolver( LinearSolverParameters const & parameters )
  :
  m_parameters( parameters )
{ }


//// ----------------------------
//// Top-Level Solver
//// ----------------------------
//// We switch between different solverTypes here
//
void HypreSolver::solve( HypreMatrix & mat,
                         HypreVector & sol,
                         HypreVector & rhs )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  if( m_parameters.solverType == "direct" )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs );
  }
}

void HypreSolver::solve_direct( HypreMatrix & mat,
                                HypreVector & sol,
                                HypreVector & rhs )
{
  // Instantiate solver
  HYPRE_Solver solver;

  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistSetup( &solver,
                                             mat.unwrappedParCSR(),
                                             0 ) );
  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistSolve( solver,
                                             rhs.unwrapped(),
                                             sol.unwrapped() ) );
  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistDestroy( solver ) );

}

void HypreSolver::solve_krylov( HypreMatrix & mat,
                                HypreVector & sol,
                                HypreVector & rhs )
{

  // Instantiate preconditioner and solver
  HYPRE_Solver precond = nullptr;
  HYPRE_Solver solver;

  // Get MPI communicator
  MPI_Comm comm = mat.getComm();


  // Setup the preconditioner
  HYPRE_Int (* precondSetupFunction)( HYPRE_Solver,
                                      HYPRE_ParCSRMatrix,
                                      HYPRE_ParVector,
                                      HYPRE_ParVector ) = nullptr;
  HYPRE_Int (* precondApplyFunction)( HYPRE_Solver,
                                      HYPRE_ParCSRMatrix,
                                      HYPRE_ParVector,
                                      HYPRE_ParVector ) = nullptr;
  HYPRE_Int (* precondDestroyFunction)( HYPRE_Solver ) = nullptr;


  if( m_parameters.preconditionerType == "none" )
  {
    precondSetupFunction = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentitySetup;
    precondApplyFunction = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentity;
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
    precondSetupFunction = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentitySetup;
    precondApplyFunction = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentity;
  }
  else if( m_parameters.preconditionerType == "ilu" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &precond ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( precond, 0.0 ) );

    if( m_parameters.ilu.fill >= 0 )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( precond,
                                                      LvArray::integerConversion< HYPRE_Int >( m_parameters.ilu.fill ) ) );
    }
    precondSetupFunction = (HYPRE_PtrToParSolverFcn) HYPRE_ILUSetup;
    precondApplyFunction = (HYPRE_PtrToParSolverFcn) HYPRE_ILUSolve;
    precondDestroyFunction = (HYPRE_PtrToSolverDestroyFcn) HYPRE_ILUDestroy;
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &precond ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( precond, 0.0 ) );

    if( m_parameters.ilu.fill >= 0 )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( precond,
                                                      LvArray::integerConversion< HYPRE_Int >( m_parameters.ilu.fill ) ) );
    }
    if( m_parameters.ilu.threshold >= 0 )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetDropThreshold( precond,
                                                        m_parameters.ilu.threshold ) );
    }
    precondSetupFunction = (HYPRE_PtrToParSolverFcn) HYPRE_ILUSetup;
    precondApplyFunction = (HYPRE_PtrToParSolverFcn) HYPRE_ILUSolve;
    precondDestroyFunction = (HYPRE_PtrToSolverDestroyFcn) HYPRE_ILUDestroy;
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( precond, -1 ) ); /* print amg solution info */
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( precond, 6 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetOldDefault( precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 6 ) );   /* Sym G.S./Jacobi hybrid */
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( precond, 0.0 ) );       /* conv. tolerance zero */
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( precond, 1 ) );     /* do only one iteration! */

    precondSetupFunction = (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup;
    precondApplyFunction = (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve;
    precondDestroyFunction = (HYPRE_PtrToSolverDestroyFcn) HYPRE_BoomerAMGDestroy;
  }
  else
  {
    GEOSX_ERROR( "The requested preconditionerType doesn't seem to exist" );
  }

  HYPRE_Int result = 0;

  // Choose the solver type - set parameters - solve
  if( m_parameters.solverType == "gmres" )
  {

    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESCreate( comm, &solver ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetMaxIter( solver, m_parameters.krylov.maxIterations ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetKDim( solver, m_parameters.krylov.maxRestart ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetTol( solver, m_parameters.krylov.tolerance ) );

    // Default for now
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrintLevel( solver, m_parameters.logLevel ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetLogging( solver, 1 ) );

    // Set the preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrecond( solver,
                                                        precondApplyFunction,
                                                        precondSetupFunction,
                                                        precond ) );

    // Setup
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetup( solver,
                                                   mat.unwrappedParCSR(),
                                                   rhs.unwrapped(),
                                                   sol.unwrapped() ) );

    // Solve
    result = HYPRE_ParCSRGMRESSolve( solver,
                                     mat.unwrappedParCSR(),
                                     rhs.unwrapped(),
                                     sol.unwrapped() );

    // Clear error code to avoid GEOSX from crashing if Krylov method did not converge
    GEOSX_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

    // Destroy solver
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESDestroy( solver ) );
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABCreate( comm, &solver ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetMaxIter( solver, m_parameters.krylov.maxIterations ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetTol( solver, m_parameters.krylov.tolerance ) );

    // Default for now
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetPrintLevel( solver, m_parameters.logLevel ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetLogging( solver, 1 ) );

    // Set the preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetPrecond( solver,
                                                     precondApplyFunction,
                                                     precondSetupFunction,
                                                     precond ) );

    // Setup
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetup( solver,
                                                      mat.unwrappedParCSR(),
                                                      rhs.unwrapped(),
                                                      sol.unwrapped() ) );

    // Solve
    result = HYPRE_ParCSRBiCGSTABSolve( solver,
                                        mat.unwrappedParCSR(),
                                        rhs.unwrapped(),
                                        sol.unwrapped() );

    // Clear error code to avoid GEOSX from crashing if Krylov method did not converge
    GEOSX_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

    // Destroy solver
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABDestroy( solver ) );
  }
  else if( m_parameters.solverType == "cg" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGCreate( comm, &solver ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetMaxIter( solver, m_parameters.krylov.maxIterations ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTol( solver, m_parameters.krylov.tolerance ) );

    // Default for now
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetPrintLevel( solver, m_parameters.logLevel ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetLogging( solver, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTwoNorm( solver, 1 ) ); /* use the two norm as the stopping criteria */

    // Set the preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGSetPrecond( solver,
                                                      precondApplyFunction,
                                                      precondSetupFunction,
                                                      precond ) );

    // Setup
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGSetup( solver,
                                                 mat.unwrappedParCSR(),
                                                 rhs.unwrapped(),
                                                 sol.unwrapped() ) );

    // Solve
    result = HYPRE_ParCSRPCGSolve( solver,
                                   mat.unwrappedParCSR(),
                                   rhs.unwrapped(),
                                   sol.unwrapped() );

    // Clear error code to avoid GEOSX from crashing if Krylov method did not converge
    GEOSX_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

    // Destroy solver
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGDestroy( solver ) );
  }
  else
  {
    GEOSX_ERROR( "The requested linear solverType doesn't seem to exist" );
  }

  GEOSX_WARNING_IF( result, "HypreSolver: Krylov convergence not achieved" );

  // Destroy preconditioner
  if( precond != nullptr )
  {
    GEOSX_LAI_CHECK_ERROR( precondDestroyFunction( precond ) );
  }

}

} // end geosx namespace
