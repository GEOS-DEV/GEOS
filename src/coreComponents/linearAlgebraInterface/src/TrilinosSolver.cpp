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

#include "TrilinosSolver.hpp"

namespace geosx
{
/**
 * @brief Empty solver constructor.
 *
 * Create an empty (distributed) vector.
 */
TrilinosSolver::TrilinosSolver()
{}

/**
 * @brief Solve system.
 *
 * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector. Prec is an optional
 * pointer to a preconditioner.
 */
void TrilinosSolver::solve( EpetraSparseMatrix &Mat,
                            EpetraVector &sol,
                            EpetraVector &rhs,
                            integer max_iter,
                            real64 newton_tol,
                            std::unique_ptr<Epetra_Operator> Prec )
{
  // Create Epetra linear problem and instantiate solver.
  Epetra_LinearProblem problem( Mat.getPointer(), sol.getPointer(), rhs.getPointer());
  AztecOO solver( problem );
  if( Prec != nullptr )
  {
    solver.SetPrecOperator( Prec.get());
  }
  else
  {
    // Other parameters used to debug
//    solver.SetAztecOption( AZ_solver, AZ_cg );
//    solver.SetAztecOption( AZ_precond, AZ_Jacobi );
//    solver.SetAztecOption( AZ_conv, AZ_noscaled );

    // Choose the solver, preconditioner and options (HARD CODED)
    solver.SetAztecOption( AZ_solver, AZ_gmres );
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_ilut );
    solver.SetAztecParam( AZ_ilut_fill, 3.0 );

  }
  // Ask for a convergence normalized by the rhs
  solver.SetAztecOption( AZ_conv, AZ_rhs );
  // Set output (TODO verbosity manager?)
  solver.SetAztecOption( AZ_output, AZ_last );
  // Solve
  solver.Iterate( max_iter, newton_tol );
}

/**
 * @brief Solve system using an ML preconditioner.
 *
 * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
 * This function is a very early design and should be further improved
 * before real usage.
 *
 */
void TrilinosSolver::ml_solve( EpetraSparseMatrix &Mat,
                               EpetraVector &sol,
                               EpetraVector &rhs,
                               integer max_iter,
                               real64 newton_tol,
                               std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec )
{
  // Create Epetra linear problem and instantiate solver.
  Epetra_LinearProblem problem( Mat.getPointer(), sol.getPointer(), rhs.getPointer());
  AztecOO solver( problem );

  Teuchos::ParameterList MLList;
  ML_Epetra::SetDefaults( "SA", MLList );

  if( MLPrec == nullptr )
  {
    MLPrec = std::make_unique<ML_Epetra::MultiLevelPreconditioner>( *Mat.getPointer(), MLList );
  }
  solver.SetAztecOption( AZ_output, 0 );
  solver.SetPrecOperator( MLPrec.get());
  solver.Iterate( max_iter, newton_tol );
}

/**
 * @brief Solve system using a direct solver (sequential!).
 *
 * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVectors.
 *
 */

void TrilinosSolver::dsolve( EpetraSparseMatrix &Mat,
                             EpetraVector &sol,
                             EpetraVector &rhs )
{
  // Create Epetra linear problem and instantiate solver.
  Epetra_LinearProblem problem( Mat.getPointer(), sol.getPointer(), rhs.getPointer());
  Amesos_BaseSolver* solver;
  Amesos Factory;

  // Select KLU solver (only one available as of 9/20/2018)
  solver = Factory.Create( "Klu", problem );
  // Factorize the matrix
  solver->SymbolicFactorization();
  solver->NumericFactorization();
  // Solve the system
  solver->Solve();
  // Set output (TODO verbosity manager?)
  solver->PrintStatus();
  solver->PrintTiming();
}

}
