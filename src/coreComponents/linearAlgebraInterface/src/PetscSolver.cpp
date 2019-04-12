/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * PetscSolver.cpp
 *
 *  Created on: Feb. 8, 2019
 *  Author: Hannah Morgan
 */

// BEGIN_RST_NARRATIVE PetscSolver.rst
// ==============================
// Petsc Solver
// ==============================
// This class implements solvers from the PETSc library. 

// Include the corresponding header file.
#include "PetscSolver.hpp"

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

void PetscSolver::solve( PetscMatrix &mat,
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

void PetscSolver::solve_direct( MPI_Comm const comm,
                                PetscMatrix &mat,
                                PetscVector &sol,
                                PetscVector &rhs )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);

  KSPSetOperators(ksp, mat.getMat(), mat.getMat());
  KSPSetType(ksp, KSPPREONLY);

  // use direct solve preconditioner
  PC prec;
  KSPGetPC(ksp, &prec);
  PCSetType(prec, PCLU);
  PCFactorSetMatSolverType(prec, MATSOLVERMUMPS);
  // more MUMPS parameters can go here

  // solve system
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}


// ----------------------------
// Iterative solver
// ----------------------------

void PetscSolver::solve_krylov( MPI_Comm const comm,
                                PetscMatrix &mat,
                                PetscVector &sol,
                                PetscVector &rhs )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);
  KSPSetOperators(ksp, mat.getMat(), mat.getMat());
  KSPGMRESSetRestart(ksp, m_parameters.krylov.maxRestart);
  KSPSetTolerances(ksp, m_parameters.krylov.tolerance, PETSC_DEFAULT, PETSC_DEFAULT, m_parameters.krylov.maxIterations);

  // Choose the solver type
  if( m_parameters.solverType == "gmres" )
  {
    KSPSetType(ksp, KSPGMRES);
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    KSPSetType(ksp, KSPBCGS);
  }
  else if( m_parameters.solverType == "cg" )
  {
    KSPSetType(ksp, KSPCG);
  }
  else
    GEOS_ERROR( "The requested linear solverType doesn't seem to exist" );

  // Create a preconditioner
  PC prec;
  KSPGetPC(ksp, &prec);

  // Choose the preconditioner type
  if( m_parameters.preconditionerType == "none" )
  {
    PCSetType(prec, PCNONE);
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
    PCSetType(prec, PCJACOBI);
  }
  else if( m_parameters.preconditionerType == "ilu" )
  {
    GEOS_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "icc" )
  {
    GEOS_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
    PCSetType(prec, PCHYPRE);
    PCHYPRESetType(prec, "pilut");
    // HYPRE Pilut Options
    // -pc_hypre_pilut_maxiter <-2>: Number of iterations (None)
    // -pc_hypre_pilut_tol <-2.>: Drop tolerance (None)
    // -pc_hypre_pilut_factorrowsize <-2>: FactorRowSize (None)
  }
  else if( m_parameters.preconditionerType == "amg" )
  {

    if(m_parameters.amg.solver== "Petsc")
    {
      PCSetType(prec, PCMG);
      PetscOptionsSetValue(NULL, "-pc_mg_levels", m_parameters.amg.maxLevels);
      if (m_parameters.amg.cycleType == "V"){
        PCMGSetCycleType(prec, PC_MG_CYCLE_V);
      } else {
        PCMGSetCycleType(prec, PC_MG_CYCLE_W);
      }
      PCMGSetNumberSmooth(prec, m_parameters.amg.numSweeps); 
      // by default, Chebyshev + SOR smoother
    } 
    else if(m_parameters.amg.solver == "Hypre")
    {
      PCSetType(prec, PCHYPRE);
      PCHYPRESetType(prec, "boomeramg")
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_max_levels", m_parameters.amg.maxLevels); 
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type", m_parameters.amg.cycleType); 
      // relaxation method on one level
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_down", m_parameters.amgHypre.smootherType_up); 
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_up", m_parameters.amgHypre.smootherType_down); 
      // move between grids
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_smooth_type", m_parameters.amgHypre.coarseType_up);
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", m_parameters.amgHypre.coarseType_down);
      // number of relaxation sweeps on one level
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_grid_sweeps_down", m_parameters.amg.numSweeps); 
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_grid_sweeps_up", m_parameters.amg.numSweeps); 
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_grid_sweeps_coarse", m_parameters.amg.numSweeps); // coarsest grid
    } 
    else
    {
      GEOS_ERROR( "The requested linear AMG solver isn't available in PETSc" );
    }
  }
  else
    GEOS_ERROR( "The requested preconditionerType doesn't seem to exist" );

  // HANNAH: in PETSc? Ask for a convergence normalized by the right hand side
  KSPSetNormType(prec, KSP_NORM_PRECONDITIONED);

  // Control output
  switch( m_parameters.verbosity )
  {
  case 1:
    PetscOptionsSetValue(NULL, "-ksp_monitor_short", "");
    break;
  case 2:
    PetscOptionsSetValue(NULL, "-ksp_monitor", "");
  }

  // Actually solve
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}

} // end geosx namespace

