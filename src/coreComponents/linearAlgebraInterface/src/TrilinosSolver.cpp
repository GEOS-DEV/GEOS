/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * TrilinosSolver.cpp
 *
 *  Created on: Aug 9, 2018
 *  Author: Matthias Cremon
 */

// BEGIN_RST_NARRATIVE TrilinosSolver.rst
// ==============================
// Trilinos Solver
// ==============================
// This class implements solvers from the Trilinos library. Iterative solvers come from
// the AztecOO package, and direct solvers from the Amesos package.

// Include the corresponding header file.
#include "TrilinosSolver.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

TrilinosSolver::TrilinosSolver(LinearSolverParameters const & parameters)
  :
  m_parameters(parameters)
{}


// ----------------------------
// Top-Level Solver
// ----------------------------
// We switch between different solverTypes here

void TrilinosSolver::solve( EpetraMatrix &mat,
                            EpetraVector &sol,
                            EpetraVector &rhs)
{
  if(m_parameters.solverType == "direct")
    solve_direct(mat,sol,rhs);
  else
    solve_krylov(mat,sol,rhs);
}
 
// ----------------------------
// Direct solver
// ----------------------------

void TrilinosSolver::solve_direct( EpetraMatrix &mat,
                                   EpetraVector &sol,
                                   EpetraVector &rhs )
{
  // Create Epetra linear problem and instantiate solver.
  Epetra_LinearProblem problem( mat.unwrappedPointer(), sol.unwrappedPointer(), rhs.unwrappedPointer());

  // Instantiate the Amesos solver.
  Amesos_BaseSolver* solver;
  Amesos Factory;

  // Select KLU solver (only one available as of 9/20/2018)
  solver = Factory.Create( "Klu", problem );

  // Factorize the matrix
  solver->SymbolicFactorization();
  solver->NumericFactorization();

  // Solve the system
  solver->Solve();

  // Basic output
  if(m_parameters.verbosity > 0)
  {
   solver->PrintStatus();
   solver->PrintTiming();
  }
}


// ----------------------------
// Iterative solver
// ----------------------------

void TrilinosSolver::solve_krylov( EpetraMatrix &mat,
                                   EpetraVector &sol,
                                   EpetraVector &rhs)
{
  // Create Epetra linear problem.
  Epetra_LinearProblem problem( mat.unwrappedPointer(), sol.unwrappedPointer(), rhs.unwrappedPointer());

  // Instantiate the AztecOO solver.
  AztecOO solver( problem );

  // Choose the solver type
  if(m_parameters.solverType == "gmres")
  {
    solver.SetAztecOption( AZ_solver, AZ_gmres );
    solver.SetAztecOption( AZ_kspace, m_parameters.krylov.maxRestart );
  }
  else if(m_parameters.solverType == "bicgstab")
  {
    solver.SetAztecOption( AZ_solver, AZ_bicgstab );
  }
  else if(m_parameters.solverType == "cg")
  {
    solver.SetAztecOption( AZ_solver, AZ_cg );
  }
  else 
    GEOS_ERROR("The requested linear solverType doesn't seem to exist");

  // Choose the preconditioner type 
  if(m_parameters.preconditionerType == "none")
  {
    solver.SetAztecOption( AZ_precond, AZ_none );
  }
  else if(m_parameters.preconditionerType == "jacobi")
  {
    solver.SetAztecOption( AZ_precond, AZ_Jacobi );
  }
  else if(m_parameters.preconditionerType == "ilu")
  {
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_overlap, 0 );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_ilu );
    solver.SetAztecOption( AZ_graph_fill, m_parameters.ilu.fill );
  }
  else if(m_parameters.preconditionerType == "icc")
  {
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_overlap, 0 );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_icc );
    solver.SetAztecOption( AZ_graph_fill, m_parameters.ilu.fill );
  }
  else if(m_parameters.preconditionerType == "ilut")
  {
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_overlap, 0 );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_ilut );
    solver.SetAztecParam( AZ_ilut_fill, (m_parameters.ilu.fill>0?real64(m_parameters.ilu.fill):1.0)); 
  }
  else 
    GEOS_ERROR("The requested preconditionerType doesn't seem to exist");

  // Ask for a convergence normalized by the right hand side
  solver.SetAztecOption( AZ_conv, AZ_rhs );

  // Control output
  switch(m_parameters.verbosity)
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
  solver.Iterate( m_parameters.krylov.maxIterations, 
                  m_parameters.krylov.tolerance );
}

/**
 * @brief Solve system using an ML preconditioner. (Testing purposes mostly)
 *
 * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
 * This function is a very early design and should be further improved
 * before real usage.
 *
 */

// ----------------------------
// Iterative solver with ML
// ----------------------------
// Initial test of an ML-based preconditioner.
#if 0
void TrilinosSolver::ml_solve( EpetraMatrix &Mat,
                               EpetraVector &sol,
                               EpetraVector &rhs,
                               integer const max_iter,
                               real64 const newton_tol,
                               std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec )
{
  // Create Epetra linear problem.
  Epetra_LinearProblem problem( Mat.unwrappedPointer(), sol.unwrappedPointer(), rhs.unwrappedPointer());

  // Instantiate AztecOO solver.
  AztecOO solver( problem );

  // Use a default parameter list from Teuchos
  Teuchos::ParameterList MLList;
  ML_Epetra::SetDefaults( "SA", MLList );

  // Compute the preconditioner
  if( MLPrec == nullptr )
  {
    MLPrec = std::make_unique<ML_Epetra::MultiLevelPreconditioner>( *Mat.unwrappedPointer(), MLList );
  }
  // Set output (TODO verbosity manager?)
  solver.SetAztecOption( AZ_output, AZ_last );

  // Set the precondtioner
  solver.SetPrecOperator( MLPrec.get());

  // Solver the system
  solver.Iterate( max_iter, newton_tol );
}
#endif

} // end geosx namespace

