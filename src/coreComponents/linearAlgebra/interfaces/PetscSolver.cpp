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

// BEGIN_RST_NARRATIVE PetscSolver.rst
// ==============================
// Petsc Solver
// ==============================
// This class implements solvers from the PETSc library. 

// Include the corresponding header file.
#include "PetscSolver.hpp"

#include "linearAlgebra/interfaces/PetscVector.hpp"
#include "linearAlgebra/interfaces/PetscSparseMatrix.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

PetscSolver::PetscSolver( LinearSolverParameters const & parameters )
  :
  m_parameters( parameters )
{}

// ----------------------------
// Top-Level Solver
// ----------------------------
// We switch between different solverTypes here

void PetscSolver::solve( PetscSparseMatrix &mat,
                         PetscVector &sol,
                         PetscVector &rhs )
{
  if( m_parameters.solverType == "direct" )
    solve_direct( mat, sol, rhs );
  else
    solve_krylov( mat, sol, rhs );
}

// ----------------------------
// Direct solver
// ----------------------------

void PetscSolver::solve_direct( PetscSparseMatrix &mat,
                                PetscVector &sol,
                                PetscVector &rhs )
{
  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;
  KSPCreate( comm, &ksp );
  KSPSetOperators( ksp, mat.getMat(), mat.getMat() );
  KSPSetType( ksp, KSPPREONLY );

  // use direct solve preconditioner SUPERLU DIST
  PC prec;
  KSPGetPC( ksp, &prec );
  PCSetType( prec, PCLU );
  PCFactorSetMatSolverType( prec, MATSOLVERSUPERLU_DIST );

  // solve system
  KSPSetFromOptions( ksp );
  KSPSolve( ksp, rhs.getVec(), sol.getVec() );
}


// ----------------------------
// Iterative solver
// ----------------------------

void PetscSolver::solve_krylov( PetscSparseMatrix &mat,
                                PetscVector &sol,
                                PetscVector &rhs )
{
  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;
  KSPCreate( comm, &ksp );
  KSPSetOperators( ksp, mat.getMat(), mat.getMat() );
  KSPGMRESSetRestart( ksp, m_parameters.krylov.maxRestart );
  KSPSetTolerances( ksp, m_parameters.krylov.tolerance, PETSC_DEFAULT, PETSC_DEFAULT, m_parameters.krylov.maxIterations );

  // pick the solver type
  if( m_parameters.solverType == "gmres" )
  {
    KSPSetType( ksp, KSPGMRES );
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    KSPSetType( ksp, KSPBCGS );
  }
  else if( m_parameters.solverType == "cg" )
  {
    KSPSetType( ksp, KSPCG );
  }
  else
  {
    GEOSX_ERROR( "The requested linear solverType doesn't seem to exist" );
  }
  
  // create a preconditioner and pick type
  PC prec;
  KSPGetPC( ksp, &prec );

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
    PCSetType( prec, PCHYPRE );
    PCHYPRESetType( prec, "pilut" );
#else
    GEOSX_ERROR("Can't use HYPRE through PETSc in serial");
#endif
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    std::map<std::string, std::string> translate; // maps GEOSX to PETSc syntax for Hyper options

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
    PCSetType( prec, PCHYPRE );
    PCHYPRESetType( prec, "boomeramg" );
#else
    GEOSX_ERROR("Can't use HYPRE through PETSc in serial");
#endif
    // PETSc needs char[]
    char max_levels[10], cycle_type[10], num_sweeps[10], smoother_type[30], coarse_type[30];
    sprintf( max_levels, "%d", m_parameters.amg.maxLevels );
    sprintf( cycle_type, "%s", m_parameters.amg.cycleType.c_str() );
    sprintf( num_sweeps, "%d", m_parameters.amg.maxLevels );
    sprintf( smoother_type, "%s", translate[m_parameters.amg.smootherType].c_str() );
    sprintf( coarse_type, "%s", translate[m_parameters.amg.coarseType].c_str() );

    PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_max_levels", max_levels ); 
    PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_cycle_type", cycle_type ); 
    // relaxation method
    // available in HYPRE: Jacobi, sequential-Gauss-Seidel, seqboundary-Gauss-Seidel, SOR/Jacobi backward-SOR/Jacobi, symmetric-SOR/Jacobi  
    //   l1scaled-SOR/Jacobi Gaussian-elimination, l1-Gauss-Seidel, backward-l1-Gauss-Seidel, CG, Chebyshev FCF-Jacobi, l1scaled-Jacobi
    PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_relax_type_all", smoother_type ); // default: symmetric-SOR/Jacobi
    PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_relax_type_coarse", coarse_type ); // default: Gaussian-elimination
    // number of relaxation sweeps
    PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_grid_sweeps_all", num_sweeps ); 
    PetscOptionsSetValue( nullptr, "-pc_hypre_boomeramg_grid_sweeps_coarse", num_sweeps ); // coarsest grid
  }
  else
  {
    GEOSX_ERROR( "The requested preconditionerType isn't availbe in exist" );
  }

  // display output
  if ( m_parameters.logLevel > 0 )
  {
    PetscOptionsSetValue( nullptr, "-ksp_monitor", nullptr ); 
  }

  // Actually solve
  KSPSetFromOptions( ksp ); 
  KSPSolve( ksp, rhs.getVec(), sol.getVec() );
}

} // end geosx namespace

