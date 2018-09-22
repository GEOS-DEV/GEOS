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

/**
 * \class TrilinosSolver
 * \brief This class creates and provides basic support for AztecOO, Amesos and ML libraries.
 */

class TrilinosSolver
{
public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty solver constructor.
   *
   */
  TrilinosSolver();

  /**
   * @brief Copy constructor.
   *
   */
  TrilinosSolver( const TrilinosSolver &Solver );

  /**
   * @brief Virtual destructor.
   *
   */
  virtual ~TrilinosSolver() = default;
  //@}

  //! @name Solvers
  //@{
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
  //@}

protected:

};

}

#endif /* TRILINOSSOLVER_HPP_ */
