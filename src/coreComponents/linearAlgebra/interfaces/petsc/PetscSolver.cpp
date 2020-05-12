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
 * @file PetscSolver.cpp
 */

#include "PetscSolver.hpp"

#include "linearAlgebra/interfaces/petsc/PetscMatrix.hpp"
#include "linearAlgebra/interfaces/petsc/PetscVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

// Put everything under the geosx namespace.
namespace geosx
{

PetscSolver::PetscSolver( LinearSolverParameters const & parameters )
  :
  m_parameters( parameters )
{}

void PetscSolver::solve( PetscMatrix & mat,
                         PetscVector & sol,
                         PetscVector & rhs )
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

void PetscSolver::solve_direct( PetscMatrix & mat,
                                PetscVector & sol,
                                PetscVector & rhs )
{
  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;
  GEOSX_LAI_CHECK_ERROR( KSPCreate( comm, &ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, mat.unwrapped(), mat.unwrapped() ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPPREONLY ) );

  // use direct solve preconditioner SUPERLU DIST
  PC prec;
  GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp, &prec ) );
  GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCLU ) );
  GEOSX_LAI_CHECK_ERROR( PCFactorSetMatSolverType( prec, MATSOLVERSUPERLU_DIST ) );

  // solve system
  GEOSX_LAI_CHECK_ERROR( KSPSetFromOptions( ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, rhs.unwrapped(), sol.unwrapped() ) );

  // destroy
  GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );
}

void PetscSolver::solve_krylov( PetscMatrix & mat,
                                PetscVector & sol,
                                PetscVector & rhs )
{
  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;

  // Extra scratch matrix,  needed if the separate displacement component is requested
  PetscMatrix scratch; // default constructed, does nothing

  GEOSX_LAI_CHECK_ERROR( KSPCreate( comm, &ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, mat.unwrapped(), mat.unwrapped() ) );
  GEOSX_LAI_CHECK_ERROR( KSPGMRESSetRestart( ksp, m_parameters.krylov.maxRestart ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetTolerances( ksp, m_parameters.krylov.tolerance, PETSC_DEFAULT,
                                           PETSC_DEFAULT, m_parameters.krylov.maxIterations ) );

  // pick the solver type
  if( m_parameters.solverType == "gmres" )
  {
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPGMRES ) );
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPBCGS ) );
  }
  else if( m_parameters.solverType == "cg" )
  {
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPCG ) );
  }
  else
  {
    GEOSX_ERROR( "The requested linear solverType doesn't seem to exist" );
  }

  // create a preconditioner and pick type
  PC prec;
  GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp, &prec ) );

  if( m_parameters.preconditionerType == "none" )
  {
    PCSetType( prec, PCNONE );
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
    PCSetType( prec, PCJACOBI );
  }
  else if( m_parameters.preconditionerType == "icc" )
  {
    PCSetType( prec, PCICC );
  }
  else if( m_parameters.preconditionerType == "iluk" )
  {
    // Set up additive Schwartz outer preconditioner
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCASM ) );
    GEOSX_LAI_CHECK_ERROR( PCASMSetOverlap( prec, m_parameters.dd.overlap ) );
    GEOSX_LAI_CHECK_ERROR( PCASMSetType( prec, PC_ASM_RESTRICT ) );
    GEOSX_LAI_CHECK_ERROR( PCSetUp( prec ) );

    // Get local preconditioning context
    KSP * ksp_local;
    PetscInt n_local, first_local;
    GEOSX_LAI_CHECK_ERROR( PCASMGetSubKSP( prec, &n_local, &first_local, &ksp_local ) );

    // Sanity checks
    GEOSX_LAI_ASSERT_EQ( n_local, 1 );
    GEOSX_LAI_ASSERT_EQ( first_local, MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) );

    // Set up local block ILU preconditioner
    PC prec_local;
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp_local[0], KSPPREONLY ) );
    GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp_local[0], &prec_local ) );
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec_local, PCILU ) );
    GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, m_parameters.ilu.fill ) );
    GEOSX_LAI_CHECK_ERROR( PCSetUpOnBlocks( prec ) );
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
    GEOSX_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    // apply separate displacement component filter
    if( m_parameters.amg.separateComponents )
    {
      LAIHelperFunctions::SeparateComponentFilter( mat, scratch, m_parameters.dofsPerNode );
      GEOSX_LAI_CHECK_ERROR( KSPGetOperators( ksp, &mat.unwrapped(), nullptr ) );
      GEOSX_LAI_CHECK_ERROR( PetscObjectReference((PetscObject)mat.unwrapped()); );
      GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, mat.unwrapped(), scratch.unwrapped() ) );
    }

    //  Default options only for the moment
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCGAMG ) );

#if 0
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCHMG ) );
    GEOSX_LAI_CHECK_ERROR( PCHMGSetInnerPCType( prec, PCGAMG ) );

    // Set maximum number of multigrid levels
    if( m_parameters.amg.maxLevels > 0 )
    {
      GEOSX_LAI_CHECK_ERROR( PCMGSetLevels( prec,
                                            LvArray::integerConversion< PetscInt >( m_parameters.amg.maxLevels ),
                                            nullptr ) );
    }

    // Set the number of sweeps
    if( m_parameters.amg.numSweeps > 1 )
    {
      GEOSX_LAI_CHECK_ERROR( PCMGSetNumberSmooth( prec,
                                                  LvArray::integerConversion< PetscInt >( m_parameters.amg.numSweeps ) ) );
    }

    // Set type of cycle (1: V-cycle (default); 2: W-cycle)
    if( m_parameters.amg.cycleType == "V" )
    {
      GEOSX_LAI_CHECK_ERROR( PCMGSetCycleType( prec, PC_MG_CYCLE_V ) );
    }
    else if( m_parameters.amg.cycleType == "W" )
    {
      GEOSX_LAI_CHECK_ERROR( PCMGSetCycleType( prec, PC_MG_CYCLE_W ) );
    }

    // Set smoother to be used (for all levels)
    PetscInt numLevels;
    PetscInt l;
    KSP smoother;
    PC smootherPC;
    GEOSX_LAI_CHECK_ERROR( PCMGGetLevels( prec, &numLevels ) );

    GEOSX_LOG_RANK_VAR( numLevels );

    for( l = 0; l < numLevels; ++l )
    {
      GEOSX_LAI_CHECK_ERROR( PCMGGetSmoother( prec, l, &smoother ) );
      GEOSX_LAI_CHECK_ERROR( KSPSetType( smoother, KSPRICHARDSON ) );
      GEOSX_LAI_CHECK_ERROR( KSPGetPC( smoother, &smootherPC ) );

      if( m_parameters.amg.smootherType == "jacobi" )
      {
        GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCJACOBI ) );
      }
      else if( m_parameters.amg.smootherType == "gaussSeidel" )
      {
        GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCSOR ) );
      }
      else if( m_parameters.amg.smootherType.substr( 0, 3 ) == "ilu" )
      {
        // Set up additive Schwartz preconditioner
        GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCASM ) );
        GEOSX_LAI_CHECK_ERROR( PCASMSetOverlap( smootherPC, 0 ) );
        GEOSX_LAI_CHECK_ERROR( PCASMSetType( smootherPC, PC_ASM_RESTRICT ) );
        // GEOSX_LAI_CHECK_ERROR( PCSetUp( smootherPC ) );

        // Get local preconditioning context
        KSP * ksp_local;
        PetscInt n_local, first_local;
        GEOSX_LAI_CHECK_ERROR( PCASMGetSubKSP( smootherPC, &n_local, &first_local, &ksp_local ) );

        // Sanity checks
        GEOSX_LAI_ASSERT_EQ( n_local, 1 );
        GEOSX_LAI_ASSERT_EQ( first_local, MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) );

        // Set up local block ILU preconditioner
        PC prec_local;
        GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp_local[0], KSPPREONLY ) );
        GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp_local[0], &prec_local ) );
        GEOSX_LAI_CHECK_ERROR( PCSetType( prec_local, PCILU ) );
        if( m_parameters.amg.smootherType == "ilu1" )
        {
          GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, 1 ) );
        }
        else
        {
          GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, 0 ) );
        }
        // GEOSX_LAI_CHECK_ERROR( PCSetUpOnBlocks( smootherPC ) );
      }
    }

    // Set coarsest level solver
    // TODO
    // m_parameters.amg.coarseType



    // Set aggretation threshold
    // TODO
    // m_parameters.amg.aggregationThreshold

    // TODO: add user-defined null space / rigid body mode support
    // if ( m_parameters.amg.nullSpaceType )
    // {
    //   ...
    // }
#endif
  }
  else
  {
    GEOSX_ERROR( "The requested preconditioner type isn't available in PETSc" );
  }

  // display output
  if( m_parameters.logLevel > 0 )
  {
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-ksp_monitor", nullptr ) );
  }

  // Actually solve
  GEOSX_LAI_CHECK_ERROR( KSPSetFromOptions( ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, rhs.unwrapped(), sol.unwrapped() ) );

  KSPConvergedReason result;
  GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( ksp, &result ) );
  GEOSX_WARNING_IF( result < 0, "PetscSolver: Krylov convergence not achieved" );

  // reset verbosity option
  GEOSX_LAI_CHECK_ERROR( PetscOptionsClearValue( nullptr, "-ksp_monitor" ) );

  // Destroy solver
  GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );

}

} // end geosx namespace
