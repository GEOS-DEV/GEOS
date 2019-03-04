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
    PCHYPRESetType(prec, "ilut");
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    Teuchos::ParameterList list;

    if( m_parameters.amg.isSymmetric )
      ML_Petsc::SetDefaults( "SA", list );
    else
      ML_Petsc::SetDefaults( "NSSA", list );

    std::map<string, string> translate; // maps GEOSX to ML syntax

    translate.insert( std::make_pair( "V", "MGV" ));
    translate.insert( std::make_pair( "W", "MGW" ));
    translate.insert( std::make_pair( "direct", "Amesos-KLU" ));
    translate.insert( std::make_pair( "jacobi", "Jacobi" ));
    translate.insert( std::make_pair( "blockJacobi", "block Jacobi" ));
    translate.insert( std::make_pair( "gaussSeidel", "Gauss-Seidel" ));
    translate.insert( std::make_pair( "blockGaussSeidel", "block Gauss-Seidel" ));
    translate.insert( std::make_pair( "chebyshev", "Chebyshev" ));
    translate.insert( std::make_pair( "ilu", "ILU" ));
    translate.insert( std::make_pair( "ilut", "ILUT" ));

    list.set( "ML output", m_parameters.verbosity );
    list.set( "max levels", m_parameters.amg.maxLevels );
    list.set( "aggregation: type", "Uncoupled" );
    list.set( "PDE equations", m_parameters.dofsPerNode );
    list.set( "smoother: sweeps", m_parameters.amg.numSweeps );
    list.set( "prec type", translate[m_parameters.amg.cycleType] );
    list.set( "smoother: type", translate[m_parameters.amg.smootherType] );
    list.set( "coarse: type", translate[m_parameters.amg.coarseType] );

    //TODO: add user-defined null space / rigid body mode support
    //list.set("null space: type","pre-computed");
    //list.set("null space: vectors",&rigid_body_modes[0]);
    //list.set("null space: dimension", n_rbm);

    ml_preconditioner.reset( new ML_Petsc::MultiLevelPreconditioner( *mat.unwrappedPointer(), list ));
    solver.SetPrecOperator( ml_preconditioner.get() );
  }
  else
    GEOS_ERROR( "The requested preconditionerType doesn't seem to exist" );

  // HANNAH: in PETSc? Ask for a convergence normalized by the right hand side
  KSPSetNormType(prec, KSP_NORM_PRECONDITIONED);

  // Control output
  switch( m_parameters.verbosity )
  {
  case 1:
    solver.SetAztecOption( AZ_output, AZ_summary );
    solver.SetAztecOption( AZ_diagnostics, AZ_all );
    break;
  case 2:
    solver.SetAztecOption( AZ_output, AZ_all );
    solver.SetAztecOption( AZ_diagnostics, AZ_all );
    break;
  default:
    solver.SetAztecOption( AZ_output, AZ_none );
  }

  // Actually solve
  KSPSolve(ksp, rhs.getVec(), sol.getVec());

  //TODO: should we return performance feedback to have GEOSX pretty print details?:
  //      i.e. iterations to convergence, residual reduction, etc.
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
