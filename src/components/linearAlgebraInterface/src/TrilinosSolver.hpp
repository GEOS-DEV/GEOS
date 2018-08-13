/*
 * TrilinosSolver.hpp
 *
 *  Created on: Aug 9, 2018
 *      Author: Matthias
 */

#ifndef TRILINOSSOLVER_HPP_
#define TRILINOSSOLVER_HPP_

#include "EpetraSparseMatrix.hpp"
#include "EpetraVector.hpp"
#include <AztecOO.h>
#include <Amesos.h>

namespace geosx
{

class TrilinosSolver
{
public:

  /**
   * @brief Empty solver constructor.
   *
   * Create an empty (distributed) vector.
   */
  TrilinosSolver();

  /**
   * @brief Virtual destructor.
   */
  virtual ~TrilinosSolver() = default;

  /**
   * @brief Solve system.
   *
   * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
   */
  void solve( EpetraSparseMatrix &Mat,
              EpetraVector &rhs,
              EpetraVector &sol,
              integer max_iter,
              real64 newton_tol );

  /**
   * @brief Solve system using a direct solver.
   *
   * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
   */
  void dsolve( EpetraSparseMatrix &Mat,
               EpetraVector &rhs,
               EpetraVector &sol );

protected:
  Epetra_LinearProblem problem;
};

}

#endif /* TRILINOSSOLVER_HPP_ */
