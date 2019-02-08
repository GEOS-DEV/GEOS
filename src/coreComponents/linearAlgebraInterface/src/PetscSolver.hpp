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

/**
 * @file PetscSolver.hpp
 *
 *  Created on: Feb. 8, 2019
 *      Author: Hannah Morgan
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCSOLVER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCSOLVER_HPP_

#include "PetscSparseMatrix.hpp"
#include "PetscVector.hpp"
#include <petscksp.h>

namespace geosx
{

/**
 * \class PetscSolver
 * \brief This class creates and provides basic support for AztecOO, Amesos and ML libraries.
 */

class PetscSolver
{
public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty solver constructor.
   *
   */
  PetscSolver();

  /**
   * @brief Copy constructor.
   *
   */
  PetscSolver( PetscSolver const &Solver );

  /**
   * @brief Virtual destructor.
   *
   */
  virtual ~PetscSolver() = default;
  //@}

  //! @name Solvers
  //@{
  /**
   * @brief Solve system with an iterative solver (HARD CODED PARAMETERS, GMRES).
   *
   * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector.
   */
  void solve( MPI_Comm const comm,
              PetscSparseMatrix &Mat,
              PetscVector &sol,
              PetscVector &rhs,
              integer const max_iter,
              real64 const newton_tol,
              PC Prec = nullptr );

  /**
   * @brief Solve system using the ml preconditioner.
   *
   * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector.
   */
  void ml_solve( MPI_Comm const comm,
                 PetscSparseMatrix &Mat,
                 PetscVector &sol,
                 PetscVector &rhs,
                 integer const max_iter,
                 real64 const newton_tol,
                 PC MLPrec = nullptr );

  /**
   * @brief Solve system using a direct solver (KLU).
   *
   * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector.
   */
  void dsolve( MPI_Comm const comm,
               PetscSparseMatrix &Mat,
               PetscVector &sol,
               PetscVector &rhs );
  //@}

protected:

};

}

#endif /* PETSCSOLVER_HPP_ */
