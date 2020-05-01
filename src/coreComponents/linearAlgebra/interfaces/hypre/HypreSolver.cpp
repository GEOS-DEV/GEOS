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

typedef HYPRE_Int (* HYPRE_PtrToDestroyFcn)( HYPRE_Solver );

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
                                      HYPRE_Matrix,
                                      HYPRE_Vector,
                                      HYPRE_Vector ) = nullptr;
  HYPRE_Int (* precondApplyFunction)( HYPRE_Solver,
                                      HYPRE_Matrix,
                                      HYPRE_Vector,
                                      HYPRE_Vector ) = nullptr;
  HYPRE_Int (* precondDestroyFunction)( HYPRE_Solver ) = nullptr;


  if( m_parameters.preconditionerType == "none" )
  {
    GEOSX_ERROR( "precond none: Not implemented yet" );
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
    GEOSX_ERROR( "precond Jacobi: Not implemented yet" );
  }
  else if( m_parameters.preconditionerType == "ilu" )
  {
    GEOSX_ERROR( "precond ilu: Not implemented yet" );
  }
  else if( m_parameters.preconditionerType == "icc" )
  {
    GEOSX_ERROR( "precond icc: Not implemented yet" );
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_EuclidCreate( comm, &precond ) );
    precondApplyFunction = (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve;
    precondSetupFunction = (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup;
    precondDestroyFunction = (HYPRE_PtrToDestroyFcn) HYPRE_EuclidDestroy;
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

    precondApplyFunction = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve;
    precondSetupFunction = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup;
    precondDestroyFunction = (HYPRE_PtrToDestroyFcn) HYPRE_BoomerAMGDestroy;
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_GMRESSetMaxIter( solver, m_parameters.krylov.maxIterations ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_GMRESSetKDim( solver, m_parameters.krylov.maxRestart ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_GMRESSetTol( solver, m_parameters.krylov.tolerance ) );

    // Default for now
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_GMRESSetPrintLevel( solver, m_parameters.logLevel ) ); /* prints out the iteration info */
    GEOSX_LAI_CHECK_ERROR( HYPRE_GMRESSetLogging( solver, 1 ) ); /* needed to get run info later */

    // Set the preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_GMRESSetPrecond( solver,
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_BiCGSTABSetMaxIter( solver, m_parameters.krylov.maxIterations ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BiCGSTABSetTol( solver, m_parameters.krylov.tolerance ) );

    // Default for now
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_BiCGSTABSetPrintLevel( solver, m_parameters.logLevel ) ); /* prints out the iteration info */
    GEOSX_LAI_CHECK_ERROR( HYPRE_BiCGSTABSetLogging( solver, 1 ) );    /* needed to get run info later */

    // Set the preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_BiCGSTABSetPrecond( solver,
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
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_PCGSetPrintLevel( solver, m_parameters.logLevel ) ); /* prints out the iteration info */
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetLogging( solver, 1 ) );    /* needed to get run info later */
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTwoNorm( solver, 1 ) );    /* use the two norm as the stopping criteria */

    // Set the preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetPrecond( solver,
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
  GEOSX_LAI_CHECK_ERROR( precondDestroyFunction( precond ) );

  //TODO: should we return performance feedback to have GEOSX pretty print details?:
  //      i.e. iterations to convergence, residual reduction, etc.
}

} // end geosx namespace
