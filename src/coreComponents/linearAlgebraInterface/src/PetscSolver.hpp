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

/**
 * @file PetscSolver.hpp
 *
 *  Created on: Feb. 8, 2019
 *      Author: Hannah Morgan
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCSOLVER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCSOLVER_HPP_

#include <petscksp.h>

#include "PetscSparseMatrix.hpp"
// #include "PetscVector.hpp"
#include "LinearSolverParameters.hpp"

namespace geosx
{

/**
 * \class PetscSolver
 * \brief This class creates and provides basic support for PETSc solvers.
 */

class PetscSolver
{
public:

  /**
   * @brief Solver constructor, with parameter list reference
   *
   */
  PetscSolver( LinearSolverParameters const & parameters );

  /**
   * @brief Virtual destructor.
   *
   */
  virtual ~PetscSolver() = default;

  /**
   * @brief Solve system with an iterative solver (HARD CODED PARAMETERS, GMRES).
   *
   * Solve Ax=b with A an PetscMatrix, x and b PetscVector.
   */

  void solve( PetscMatrix &mat,
              PetscVector &sol,
              PetscVector &rhs );

private:

  LinearSolverParameters const & m_parameters;

  void solve_direct( PetscMatrix &mat,
                     PetscVector &sol,
                     PetscVector &rhs );

  void solve_krylov( PetscMatrix &mat,
                     PetscVector &sol,
                     PetscVector &rhs );

};

} // end geosx namespace

#endif /* PETSCSOLVER_HPP_ */
