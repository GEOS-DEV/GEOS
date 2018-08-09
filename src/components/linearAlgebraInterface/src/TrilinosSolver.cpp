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
    solver.SetAztecOption(AZ_solver,AZ_gmres);
    solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
    solver.SetAztecParam(AZ_ilut_fill,5.0);
    //solver.SetAztecOption(AZ_output,0);
    solver.Iterate(Mat.getPointer(),sol.getPointer(),rhs.getPointer(),max_iter,newton_tol);
  }

}


