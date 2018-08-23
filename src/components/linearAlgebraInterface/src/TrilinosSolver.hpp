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
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

namespace geosx
{

class TrilinosSolver
{
public:

  /**
   * @brief Empty solver constructor.
   *
   */
  TrilinosSolver();

  /**
   * @brief Copy constructor.
   *
   */
  TrilinosSolver(const TrilinosSolver &Solver);

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
              real64 newton_tol,
              std::unique_ptr<Epetra_Operator> Prec = nullptr );

  /**
   * @brief Solve system using the ml preconditioner.
   *
   * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
   */
  void ml_solve( EpetraSparseMatrix &Mat,
                 EpetraVector &rhs,
                 EpetraVector &sol,
                 integer max_iter,
                 real64 newton_tol,
                 std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec = nullptr );

  /**
   * @brief Solve system using a direct solver.
   *
   * Solve Ax=b with A an EpetraSparseMatrix, x and b EpetraVector.
   */
  void dsolve( EpetraSparseMatrix &Mat,
               EpetraVector &rhs,
               EpetraVector &sol );

protected:

};

}

#endif /* TRILINOSSOLVER_HPP_ */
