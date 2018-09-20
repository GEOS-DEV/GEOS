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
  Epetra_LinearProblem problem( Mat.getPointer(), sol.getPointer(), rhs.getPointer());
  AztecOO solver( problem );
  solver.SetAztecOption( AZ_solver, AZ_cg );
  if( Prec != nullptr )
  {
    solver.SetPrecOperator( Prec.get());
  }
  else
  {
    solver.SetAztecOption( AZ_precond, AZ_Jacobi );
    solver.SetAztecOption( AZ_conv, AZ_noscaled );
//    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
//    solver.SetAztecOption( AZ_subdomain_solve, AZ_ilut );
//    solver.SetAztecParam( AZ_ilut_fill, 5.0 );
  }
  solver.SetAztecOption( AZ_output, 0 );
  solver.Iterate( max_iter, newton_tol );
}

/**
 * @brief Solve system using an ML preconditioner.
 *
 * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
 */
void TrilinosSolver::ml_solve( EpetraSparseMatrix &Mat,
                               EpetraVector &sol,
                               EpetraVector &rhs,
                               integer max_iter,
                               real64 newton_tol,
                               std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec )
{
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
 * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
 */
void TrilinosSolver::dsolve( EpetraSparseMatrix &Mat,
                             EpetraVector &sol,
                             EpetraVector &rhs )
{
  Epetra_LinearProblem problem( Mat.getPointer(), sol.getPointer(), rhs.getPointer());
  Amesos_BaseSolver* solver;
  Amesos Factory;

  solver = Factory.Create( "Klu", problem );
  solver->SymbolicFactorization();
  solver->NumericFactorization();
  solver->Solve();
}

}
