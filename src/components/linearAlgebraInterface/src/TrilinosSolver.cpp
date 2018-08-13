/*
 * TrilinosSolver.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: Matthias
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
  {

  }

  /**
   * @brief Solve system.
   *
   * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
   */
  void TrilinosSolver::solve( EpetraSparseMatrix &Mat,
                              EpetraVector &sol,
                              EpetraVector &rhs,
                              integer max_iter,
                              real64 newton_tol )
  {
    problem(Mat,sol,rhs);
    AztecOO solver(problem);
    solver.SetAztecOption(AZ_solver,AZ_gmres);
    solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
    solver.SetAztecParam(AZ_ilut_fill,5.0);
    //solver.SetAztecOption(AZ_output,0);
    solver.Iterate(Mat.getPointer(),sol.getPointer(),rhs.getPointer(),max_iter,newton_tol);
  }

  /**
   * @brief Solve system using a direct solver.
   *
   * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
   */
  void TrilinosSolver::dsolve( EpetraSparseMatrix &Mat,
                              EpetraVector &sol,
                              EpetraVector &rhs )
  {
    problem(Mat,sol,rhs);
    Amesos_BaseSolver* solver;
    Amesos Factory;
    std::string solverName = "Klu";

    solver = Factory.Create(solverName,problem);
    solver->SymbolicFactorization();
    solver->NumericFactorization();
    solver->Solve();
  }

}


