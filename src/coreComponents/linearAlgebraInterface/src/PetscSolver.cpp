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

    PCSetType(prec, GAMG);

    // Options
    if (m_parameters.amg.cycleType == "v"){
      PCMGSetCycleType(prec, PC_MG_CYCLE_V);
    } else {
      PCMGSetCycleType(prec, PC_MG_CYCLE_W);
    }
  
    // GAMG options
    // -pc_gamg_type <agg>: Type of AMG method (one of) geo agg classical (PCGAMGSetType)
    // -pc_gamg_repartition: <FALSE> Repartion coarse grids (PCGAMGSetRepartition)
    // -pc_gamg_reuse_interpolation: <FALSE> Reuse prolongation operator (PCGAMGReuseInterpolation)
    // -pc_gamg_asm_use_agg: <FALSE> Use aggregation aggregates for ASM smoother (PCGAMGASMSetUseAggs)
    // -pc_gamg_use_parallel_coarse_grid_solver: <FALSE> Use parallel coarse grid solver (otherwise put last grid on one process) (PCGAMGSetUseParallelCoarseGridSolve)
    // -pc_gamg_process_eq_limit <50>: Limit (goal) on number of equations per process on coarse grids (PCGAMGSetProcEqLim)
    // -pc_gamg_coarse_eq_limit <50>: Limit on number of equations for the coarse grid (PCGAMGSetCoarseEqLim)
    // -pc_gamg_threshold_scale <1.>: Scaling of threshold for each level not specified (PCGAMGSetThresholdScale)
    // -pc_gamg_threshold <0.>: Relative threshold to use for dropping edges in aggregation graph (PCGAMGSetThreshold)
    // -pc_mg_levels <30>: Set number of MG levels (PCGAMGSetNlevels)
    // GAMG-AGG options
    // -pc_gamg_agg_nsmooths <1>: smoothing steps for smoothed aggregation, usually 1 (PCGAMGSetNSmooths)
    // -pc_gamg_sym_graph: <FALSE> Set for asymmetric matrices (PCGAMGSetSymGraph)
    // -pc_gamg_square_graph <1>: Number of levels to square graph for faster coarsening and lower coarse grid complexity (PCGAMGSetSquareGraph)
    // Matrix/graph coarsen (MatCoarsen) options -------------------------------------------------
    // -mat_coarsen_type <mis>: Type of aggregator (one of) mis hem (MatCoarsenSetType)

    // list.set( "ML output", m_parameters.verbosity );
    // list.set( "max levels", m_parameters.amg.maxLevels );
    // list.set( "aggregation: type", "Uncoupled" );
    // list.set( "PDE equations", m_parameters.dofsPerNode );
    // list.set( "smoother: sweeps", m_parameters.amg.numSweeps );
    // list.set( "prec type", translate[m_parameters.amg.cycleType] );
    // list.set( "smoother: type", translate[m_parameters.amg.smootherType] );
    // list.set( "coarse: type", translate[m_parameters.amg.coarseType] );
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

// HANNAH: Trilinos options for ILU, ICC, ILUT

// solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
          // Domain decomposition preconditioner (additive
          // Schwarz). That is, each processor augments its
          // submatrix according to options[AZ overlap] and
          // approximately “solves” the resulting subsystem
          // using the solver specified by
          // options[AZ subdomain solve].
          // Note: options[AZ reorder] determines whether
          // matrix equations are reordered (RCM) before
          // “solving” submatrix problem.
// solver.SetAztecOption( AZ_overlap, 0 );
          // Determines the submatrices factored with the domain
          // decomposition algorithms (see options[AZ precond]).
          // DEFAULT: 0.
// solver.SetAztecOption( AZ_subdomain_solve, AZ_icc );
          // Specifies the solver to use on each subdomain when
          // options[AZ precond] is set to AZ dom decomp DEFAULT: AZ ilut. 
          // Similar to AZ ilu using icc(k) instead of ilu(k)
// solver.SetAztecOption( AZ_graph_fill, m_parameters.ilu.fill );
          // The level of graph fill-in (k) for incomplete factorizations: ilu(k), icc(k), bilu(k). DEFAULT: 0
