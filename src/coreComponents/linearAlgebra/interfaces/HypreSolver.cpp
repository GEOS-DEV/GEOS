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

// BEGIN_RST_NARRATIVE TrilinosSolver.rst
// ==============================
// Hypre Solver
// ==============================
// ...

// Include the corresponding header file.
#include "HypreSolver.hpp"

#include "linearAlgebra/interfaces/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/HypreVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include "_hypre_utilities.h"
#include "_hypre_parcsr_ls.h"
#include "_hypre_IJ_mv.h"
#include "krylov.h"

// Put everything under the geosx namespace.
namespace geosx
{

typedef HYPRE_Int (*HYPRE_PtrToDestroyFcn)(HYPRE_Solver);

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

HypreSolver::HypreSolver( LinearSolverParameters const & parameters )
  :
  m_parameters( parameters )
{}


//// ----------------------------
//// Top-Level Solver
//// ----------------------------
//// We switch between different solverTypes here
//
void HypreSolver::solve( HypreMatrix & mat,
                         HypreVector & sol,
                         HypreVector & rhs )
{
  GEOSX_ASSERT_MSG( mat.isClosed(),
                   "Matrix has not been closed");
  GEOSX_ASSERT_MSG( sol.unwrappedPointer() != nullptr,
                   "Invalid solution vector");
  GEOSX_ASSERT_MSG( rhs.unwrappedPointer() != nullptr,
                   "Invalid right-hand side vector");

  if( m_parameters.solverType == "direct" )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs );
  }
}

// ----------------------------
// Direct solver
// ----------------------------

void HypreSolver::solve_direct( HypreMatrix & mat,
                                HypreVector & sol,
                                HypreVector & rhs )
{
  // Instantiate solver
  HYPRE_Solver solver;

  hypre_SLUDistSetup( &solver,
                      HYPRE_ParCSRMatrix(mat),
                      0 );
  hypre_SLUDistSolve( solver,
                      HYPRE_ParVector( rhs ),
                      HYPRE_ParVector( sol ) );
  hypre_SLUDistDestroy( solver );

}


// ----------------------------
// Iterative solver
// ----------------------------

void HypreSolver::solve_krylov( HypreMatrix & mat,
                                HypreVector & sol,
                                HypreVector & rhs )
{

  // Instantiate preconditioner and solver
  HYPRE_Solver precond = nullptr;
  HYPRE_Solver solver;

  // Get MPI communicator
  MPI_Comm comm = hypre_IJMatrixComm( HYPRE_IJMatrix(mat) );


  // Setup the preconditioner
  HYPRE_Int (*precondSetupFunction)( HYPRE_Solver,
                                     HYPRE_Matrix,
                                     HYPRE_Vector,
                                     HYPRE_Vector ) = nullptr;
  HYPRE_Int (*precondApplyFunction)( HYPRE_Solver,
                                     HYPRE_Matrix,
                                     HYPRE_Vector,
                                     HYPRE_Vector ) = nullptr;
  HYPRE_Int (*precondDestroyFunction)( HYPRE_Solver ) = nullptr;


  if( m_parameters.preconditionerType == "none" )
  {
    GEOSX_ERROR( "precond none: Not implemented yet" );
//    solver.SetAztecOption( AZ_precond, AZ_none );
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
//    solver.SetAztecOption( AZ_precond, AZ_Jacobi );
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
    HYPRE_ParCSRPilutCreate(comm, &precond);

    if (m_parameters.ilu.threshold >= 0 )
    {
      HYPRE_ParCSRPilutSetDropTolerance( precond,
                                         m_parameters.ilu.threshold );
    }
    if (m_parameters.ilu.fill >= 0 )
    {
      HYPRE_ParCSRPilutSetFactorRowSize( precond,
                                         integer_conversion<HYPRE_Int>(m_parameters.ilu.fill) );
    }

    precondApplyFunction = (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve;
    precondSetupFunction = (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup;
    precondDestroyFunction = (HYPRE_PtrToDestroyFcn) HYPRE_ParCSRPilutDestroy;

  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, -1); /* print amg solution info */
    HYPRE_BoomerAMGSetCoarsenType(precond, 6);
    HYPRE_BoomerAMGSetOldDefault(precond);
    HYPRE_BoomerAMGSetRelaxType(precond, 6);   /* Sym G.S./Jacobi hybrid */
    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, 0.0);       /* conv. tolerance zero */
    HYPRE_BoomerAMGSetMaxIter(precond, 1);     /* do only one iteration! */

    precondApplyFunction = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve;
    precondSetupFunction = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup;
    precondDestroyFunction = (HYPRE_PtrToDestroyFcn) HYPRE_BoomerAMGDestroy;
  }
  else
  {
    GEOSX_ERROR( "The requested preconditionerType doesn't seem to exist" );
  }

  // Choose the solver type - set parameters - solve
  if( m_parameters.solverType == "gmres" )
  {

    HYPRE_ParCSRGMRESCreate( comm, &solver);
    HYPRE_GMRESSetMaxIter(solver, m_parameters.krylov.maxIterations);
    HYPRE_GMRESSetKDim(solver, m_parameters.krylov.maxRestart);
    HYPRE_GMRESSetTol(solver, m_parameters.krylov.tolerance);

    // Default for now
    HYPRE_GMRESSetPrintLevel(solver, m_parameters.logLevel); /* prints out the iteration info */
    HYPRE_GMRESSetLogging(solver, 1); /* needed to get run info later */

    // Set the preconditioner
    HYPRE_GMRESSetPrecond( solver,
                           precondApplyFunction,
                           precondSetupFunction,
                           precond );

    // Setup
    HYPRE_ParCSRGMRESSetup( solver,
                          HYPRE_ParCSRMatrix(mat),
                          HYPRE_ParVector( rhs ),
                          HYPRE_ParVector( sol ));

    // Solve
    HYPRE_ParCSRGMRESSolve( solver,
                          HYPRE_ParCSRMatrix(mat),
                          HYPRE_ParVector( rhs ),
                          HYPRE_ParVector( sol ));


    /* Destroy solver and preconditioner */
    HYPRE_ParCSRGMRESDestroy(solver);
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    HYPRE_ParCSRBiCGSTABCreate( comm, &solver);
    HYPRE_BiCGSTABSetMaxIter(solver, m_parameters.krylov.maxIterations);
    HYPRE_BiCGSTABSetTol(solver, m_parameters.krylov.tolerance);

    // Default for now
    HYPRE_BiCGSTABSetPrintLevel(solver, 2); /* prints out the iteration info */
    HYPRE_BiCGSTABSetLogging(solver, 1);    /* needed to get run info later */

    // Set the preconditioner
    HYPRE_BiCGSTABSetPrecond( solver,
                              precondApplyFunction,
                              precondSetupFunction,
                                precond);

    // Setup
    HYPRE_ParCSRBiCGSTABSetup( solver,
                               HYPRE_ParCSRMatrix(mat),
                               HYPRE_ParVector( rhs ),
                               HYPRE_ParVector( sol ));

    // Solve
    HYPRE_ParCSRBiCGSTABSolve( solver,
                               HYPRE_ParCSRMatrix(mat),
                               HYPRE_ParVector( rhs ),
                               HYPRE_ParVector( sol ));


    /* Destroy solver and preconditioner */
    HYPRE_ParCSRBiCGSTABDestroy(solver);
  }
  else if( m_parameters.solverType == "cg" )
  {
    HYPRE_ParCSRPCGCreate( comm, &solver);
    HYPRE_PCGSetMaxIter(solver, m_parameters.krylov.maxIterations);
    HYPRE_PCGSetTol(solver, m_parameters.krylov.tolerance);

    // Default for now
    HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
    HYPRE_PCGSetLogging(solver, 1);    /* needed to get run info later */
    HYPRE_PCGSetTwoNorm(solver, 1);    /* use the two norm as the stopping criteria */

    // Set the preconditioner
    HYPRE_PCGSetPrecond( solver,
                         precondApplyFunction,
                         precondSetupFunction,
                         precond);

    // Setup
    HYPRE_ParCSRPCGSetup( solver,
                          HYPRE_ParCSRMatrix(mat),
                          HYPRE_ParVector( rhs ),
                          HYPRE_ParVector( sol ));

    // Solve
    HYPRE_ParCSRPCGSolve( solver,
                          HYPRE_ParCSRMatrix(mat),
                          HYPRE_ParVector( rhs ),
                          HYPRE_ParVector( sol ));


    /* Destroy solver and preconditioner */
    HYPRE_ParCSRPCGDestroy(solver);
  }
  else
  {
    GEOSX_ERROR( "The requested linear solverType doesn't seem to exist" );
  }

  // Destroy preconditioner
  precondDestroyFunction(precond);

//  //TODO: should we return performance feedback to have GEOSX pretty print details?:
//  //      i.e. iterations to convergence, residual reduction, etc.
}

} // end geosx namespace
