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

#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

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
}

void PetscSolver::solve_krylov( PetscMatrix & mat,
                                PetscVector & sol,
                                PetscVector & rhs )
{
  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;
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
  else if( m_parameters.preconditionerType == "ilu" )
  {
    GEOSX_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "icc" )
  {
    GEOSX_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
#ifdef GEOSX_USE_MPI
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCHYPRE ) );
    GEOSX_LAI_CHECK_ERROR( PCHYPRESetType( prec, "pilut" ) );
#else
    GEOSX_ERROR( "Can't use HYPRE through PETSc in serial" );
#endif
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    std::map< std::string, std::string > translate; // maps GEOSX to PETSc syntax for Hyper options

    translate.insert( std::make_pair( "jacobi", "Jacobi" ));
    translate.insert( std::make_pair( "sequentialGaussSeidel", "sequential-Gauss-Seidel" ));
    translate.insert( std::make_pair( "seqboundaryGaussSeidel", "seqboundary-Gauss-Seidel" ));
    translate.insert( std::make_pair( "sorJacobi", "SOR/Jacobi" ));
    translate.insert( std::make_pair( "backwardSorJacobi", "backward-SOR/Jacobi" ));
    translate.insert( std::make_pair( "symmetricSorJacobi", "symmetric-SOR/Jacobi" ));
    translate.insert( std::make_pair( "l1scaledSorJacobi", "l1scaled-SOR/Jacobi" ));
    translate.insert( std::make_pair( "direct", "Gaussian-elimination" ));
    translate.insert( std::make_pair( "l1GaussSeidel", "l1-Gauss-Seidel" ));
    translate.insert( std::make_pair( "backwardl1GaussSeidel", "backward-l1-Gauss-Seidel" ));
    translate.insert( std::make_pair( "cg", "CG" ));
    translate.insert( std::make_pair( "chebyshev", "Chebyshev" ));
    translate.insert( std::make_pair( "fcfJacobi", "FCF-Jacobi" ));
    translate.insert( std::make_pair( "l1scaledJacobi", "l1scaled-Jacobi" ));

#ifdef GEOSX_USE_MPI
    GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCHYPRE ) );
    GEOSX_LAI_CHECK_ERROR( PCHYPRESetType( prec, "boomeramg" ) );
#else
    GEOSX_ERROR( "Can't use HYPRE through PETSc in serial" );
#endif
    // PETSc needs char[]
    char max_levels[10], cycle_type[10], num_sweeps[10], smoother_type[30], coarse_type[30];
    sprintf( max_levels, "%d", m_parameters.amg.maxLevels );
    sprintf( cycle_type, "%s", m_parameters.amg.cycleType.c_str() );
    sprintf( num_sweeps, "%d", m_parameters.amg.maxLevels );
    sprintf( smoother_type, "%s", translate[m_parameters.amg.smootherType].c_str() );
    sprintf( coarse_type, "%s", translate[m_parameters.amg.coarseType].c_str() );

    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_max_levels", max_levels ) );
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_cycle_type", cycle_type ) );
    // relaxation method
    // available in HYPRE: Jacobi, sequential-Gauss-Seidel, seqboundary-Gauss-Seidel, SOR/Jacobi backward-SOR/Jacobi,
    // symmetric-SOR/Jacobi
    //   l1scaled-SOR/Jacobi Gaussian-elimination, l1-Gauss-Seidel, backward-l1-Gauss-Seidel, CG, Chebyshev FCF-Jacobi,
    // l1scaled-Jacobi
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_relax_type_all", smoother_type ) ); // default:
                                                                                                                   // symmetric-SOR/Jacobi
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_relax_type_coarse", coarse_type ) ); // default:
                                                                                                                    // Gaussian-elimination
    // number of relaxation sweeps
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_grid_sweeps_all", num_sweeps ) );
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_grid_sweeps_coarse", num_sweeps ) ); // coarsest
                                                                                                                    // grid
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
}

} // end geosx namespace
