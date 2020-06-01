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

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"
#include "linearAlgebra/interfaces/hypre/HyprePreconditioner.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_IJ_mv.h>
#include <krylov.h>

namespace geosx
{

typedef HYPRE_Int (* HYPRE_PtrToSolverDestroyFcn)( HYPRE_Solver );

HypreSolver::HypreSolver( LinearSolverParameters parameters )
  :
  m_parameters( std::move( parameters ) )
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

  Stopwatch watch;
  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistSetup( &solver,
                                             mat.unwrapped(),
                                             0 ) );
  m_result.setupTime = watch.elapsedTime();
  watch.zero();
  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistSolve( solver,
                                             rhs.unwrapped(),
                                             sol.unwrapped() ) );
  m_result.solveTime = watch.elapsedTime();

  m_result.status = LinearSolverResult::Status::Success;
  m_result.numIterations = 1;
  m_result.residualReduction = NumericTraits< real64 >::eps;

  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistDestroy( solver ) );
}

namespace
{

void CreateHypreGMRES( LinearSolverParameters const & params,
                       MPI_Comm const comm,
                       HYPRE_Solver & solver,
                       HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetKDim( solver, params.krylov.maxRestart ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrintLevel( solver, params.logLevel ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetLogging( solver, 1 ) ); /* needed to get run info later */

  solverFuncs.setPrecond = HYPRE_ParCSRGMRESSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRGMRESSetup;
  solverFuncs.solve = HYPRE_ParCSRGMRESSolve;
  solverFuncs.getNumIter = HYPRE_GMRESGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_GMRESGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRGMRESDestroy;
}

void CreateHypreBiCGSTAB( LinearSolverParameters const & params,
                          MPI_Comm const comm,
                          HYPRE_Solver & solver,
                          HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetPrintLevel( solver, params.logLevel ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetLogging( solver, 1 ) ); // needed to get run info later

  solverFuncs.setPrecond = HYPRE_ParCSRBiCGSTABSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRBiCGSTABSetup;
  solverFuncs.solve = HYPRE_ParCSRBiCGSTABSolve;
  solverFuncs.getNumIter = HYPRE_BiCGSTABGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_BiCGSTABGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRBiCGSTABDestroy;
}

void CreateHypreCG( LinearSolverParameters const & params,
                    MPI_Comm const comm,
                    HYPRE_Solver & solver,
                    HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetPrintLevel( solver, params.logLevel ) ); /* print the iteration info */
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetLogging( solver, 1 ) );    /* needed to get run info later */
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTwoNorm( solver, 1 ) );    /* use the two norm as the stopping criteria */

  solverFuncs.setPrecond = HYPRE_ParCSRPCGSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRPCGSetup;
  solverFuncs.solve = HYPRE_ParCSRPCGSolve;
  solverFuncs.getNumIter = HYPRE_PCGGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_PCGGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRPCGDestroy;
}

void CreateHypreKrylovSolver( LinearSolverParameters const & params,
                              MPI_Comm const comm,
                              HYPRE_Solver & solver,
                              HypreSolverFuncs & solverFuncs )
{
  if( params.solverType == "gmres" )
  {
    CreateHypreGMRES( params, comm, solver, solverFuncs );
  }
  else if( params.solverType == "bicgstab" )
  {
    CreateHypreBiCGSTAB( params, comm, solver, solverFuncs );
  }
  else if( params.solverType == "cg" )
  {
    CreateHypreCG( params, comm, solver, solverFuncs );
  }
  else
  {
    GEOSX_ERROR( "Unsupported Hypre solver type: " << params.solverType );
  }
}

} // namespace

void HypreSolver::solve_krylov( HypreMatrix & mat,
                                HypreVector & sol,
                                HypreVector & rhs )
{
  Stopwatch watch;

  // Create the preconditioner, but don't compute (this is done by solver setup)
  HyprePreconditioner precond( m_parameters );

  // Deal with separate component approximation
  HypreMatrix separateComponentMatrix;
  if( m_parameters.amg.separateComponents )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat, separateComponentMatrix, m_parameters.dofsPerNode );
  }
  HypreMatrix & precondMat = m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // Instantiate the solver
  HYPRE_Solver solver{};
  HypreSolverFuncs solverFuncs;
  CreateHypreKrylovSolver( m_parameters, mat.getComm(), solver, solverFuncs );

  // Set the preconditioner
  GEOSX_LAI_CHECK_ERROR( solverFuncs.setPrecond( solver,
                                                 precond.unwrappedFuncs().apply,
                                                 precond.unwrappedFuncs().setup,
                                                 precond.unwrapped() ) );

  // Setup
  GEOSX_LAI_CHECK_ERROR( solverFuncs.setup( solver,
                                            precondMat.unwrapped(),
                                            rhs.unwrapped(),
                                            sol.unwrapped() ) );
  m_result.setupTime = watch.elapsedTime();

  // Solve
  watch.zero();
  HYPRE_Int const result = solverFuncs.solve( solver,
                                              mat.unwrapped(),
                                              rhs.unwrapped(),
                                              sol.unwrapped() );
  m_result.solveTime = watch.elapsedTime();

  // Set result status based on return value
  m_result.status = result ? LinearSolverResult::Status::NotConverged : LinearSolverResult::Status::Success;

  // Clear error code to avoid GEOSX from crashing if Krylov method did not converge
  GEOSX_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

  // Get final residual norm
  HYPRE_Real finalNorm;
  GEOSX_LAI_CHECK_ERROR( solverFuncs.getFinalNorm( solver, &finalNorm ) );
  m_result.residualReduction = finalNorm;

  // Get number of iterations
  HYPRE_Int numIter;
  GEOSX_LAI_CHECK_ERROR( solverFuncs.getNumIter( solver, &numIter ) );
  m_result.numIterations = numIter;

  // Destroy solver
  GEOSX_LAI_CHECK_ERROR( solverFuncs.destroy( solver ) );
}

} // end geosx namespace
