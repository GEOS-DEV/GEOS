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

#include "interfaces/hypre/HypreMatrix.hpp"
#include "interfaces/hypre/HypreVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

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

  // Extra scratch matrix, needed if the separate displacement component is requested
  HypreMatrix scratch; // default constructed, does nothing
  HYPRE_ParCSRMatrix precondParCSRMat = mat.unwrappedParCSR();

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
    precondSetupFunction = (HYPRE_PtrToParSolverFcn) HYPRE_ParCSRDiagScaleSetup;
    precondApplyFunction = (HYPRE_PtrToParSolverFcn) HYPRE_ParCSRDiagScale;
  }
  else if( m_parameters.preconditionerType == "ilu" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &precond ) );

    // Hypre's parameters to use ParCSR ILU as a preconditioner
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

    // Hypre's parameters to use ParCSR ILU as a preconditioner
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


    // Hypre's parameters to use BoomerAMG as a preconditioner
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( precond, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( precond,
                                                         LvArray::integerConversion< HYPRE_Int >( m_parameters.amg.logLevel ) ) );;

    // Set maximum number of multigrid levels (default 25)
    if( m_parameters.amg.maxLevels > 0 )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxLevels( precond,
                                                          LvArray::integerConversion< HYPRE_Int >( m_parameters.amg.maxLevels ) ) );
    }

    // Set type of cycle (1: V-cycle (default); 2: W-cycle)
    if( m_parameters.amg.cycleType == "V" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleType( precond, 1 ) );
    }
    else if( m_parameters.amg.cycleType == "W" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleType( precond, 2 ) );
    }

    // Set smoother to be used (other options available, see hypre's documentation)
    // (default "gaussSeidel", i.e. local symmetric Gauss-Seidel)
    if( m_parameters.amg.smootherType == "jacobi" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 0 ) );
    }
    else if( m_parameters.amg.smootherType == "hybridForwardGaussSeidel" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 3 ) );
    }
    else if( m_parameters.amg.smootherType == "hybridBackwardGaussSeidel" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 4 ) );
    }
    else if( m_parameters.amg.smootherType == "hybridSymmetricGaussSeidel" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 6 ) );
    }
    else if( m_parameters.amg.smootherType == "gaussSeidel" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 6 ) );
    }
    else if( m_parameters.amg.smootherType == "L1hybridSymmetricGaussSeidel" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 8 ) );
    }
    else if( m_parameters.amg.smootherType.substr( 0, 9 ) == "chebyshev" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 16 ) );
      // Set order for Chebyshev smoother valid options 1, 2 (default), 3, 4)
      if( m_parameters.amg.smootherType == "chebyshev1" )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( precond, 1 ) );
      }
      else if( m_parameters.amg.smootherType == "chebyshev3" )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( precond, 3 ) );
      }
      else if( m_parameters.amg.smootherType == "chebyshev4" )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( precond, 4 ) );
      }
    }
    else if( m_parameters.amg.smootherType == "L1jacobi" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond, 18 ) );
    }
    else if( m_parameters.amg.smootherType.substr( 0, 3 ) == "ilu" )
    {
      // - block Jacobi ILU0 (default)
      // - block Jacobi ILU1
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( precond, 5 ) );
      GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( precond, 0 ) );
      if( m_parameters.amg.smootherType == "ilu1" )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( precond, 1 ) );
      }
    }

    // Set coarsest level solver
    // (by default for coarsest grid size above 5,000 superlu_dist is used)
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetDSLUThreshold( precond, 5000 ) );
    if( m_parameters.amg.coarseType == "direct" )
    {
      GEOSX_LAI_CHECK_ERROR( hypre_BoomerAMGSetCycleRelaxType( precond, 9, 3 ) );
    }

    // Set the number of sweeps
    if( m_parameters.amg.preOrPostSmoothing == "both" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( precond,
                                                          LvArray::integerConversion< HYPRE_Int >( m_parameters.amg.numSweeps ) ) );
    }

    // Set the number of pre-smoothing sweeps
    if( m_parameters.amg.preOrPostSmoothing == "pre" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond,
                                                               LvArray::integerConversion< HYPRE_Int >( m_parameters.amg.numSweeps ),
                                                               1 ) );
    }

    // Set the number of post-smoothing sweeps
    if( m_parameters.amg.preOrPostSmoothing == "post" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond,
                                                               LvArray::integerConversion< HYPRE_Int >( m_parameters.amg.numSweeps ),
                                                               2 ) );
    }

    // Set strength of connection
    if( ( 0.0 < m_parameters.amg.strenghtOfConnection ) && ( m_parameters.amg.strenghtOfConnection < 1.0 ) )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( precond, m_parameters.amg.strenghtOfConnection ) );
    }

//    //TODO: add user-defined null space / rigid body mode support
//    if ( m_parameters.amg.nullSpaceType )
//    {
//      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVectors ( precond ,
//                                                               HYPRE_Int num_vectors ,
//                                                               HYPRE_ParVector *interp_vectors ) );

    // apply separate displacement component filter
    if( m_parameters.amg.separateComponents )
    {
      LAIHelperFunctions::SeparateComponentFilter< HypreInterface >( mat, scratch, m_parameters.dofsPerNode );
      precondParCSRMat = scratch.unwrappedParCSR();
    }

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
                                                   precondParCSRMat,
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
                                                      precondParCSRMat,
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
                                                 precondParCSRMat,
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
