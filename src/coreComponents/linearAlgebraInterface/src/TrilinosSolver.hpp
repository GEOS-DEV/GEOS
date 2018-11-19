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
 * @file TrilinosSolver.hpp
 *
 *  Created on: Aug 9, 2018
 *      Author: Matthias Cremon
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSSOLVER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSSOLVER_HPP_

#include <AztecOO.h>
#include <Amesos.h>
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "EpetraMatrix.hpp"
#include "EpetraVector.hpp"

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
  TrilinosSolver( TrilinosSolver const &Solver );

  /**
   * @brief Virtual destructor.
   *
   */
  virtual ~TrilinosSolver() = default;
  //@}

  //! @name Solvers
  //@{
  /**
   * @brief Solve system with an iterative solver (HARD CODED PARAMETERS, GMRES).
   *
   * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
   */
  void solve( EpetraMatrix &Mat,
              EpetraVector &sol,
              EpetraVector &rhs,
              integer const max_iter,
              real64 const newton_tol,
              std::unique_ptr<Epetra_Operator> Prec = nullptr );

  /**
   * @brief Solve system using the ml preconditioner.
   *
   * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
   */
  void ml_solve( EpetraMatrix &Mat,
                 EpetraVector &sol,
                 EpetraVector &rhs,
                 integer const max_iter,
                 real64 const newton_tol,
                 std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec = nullptr );

  /**
   * @brief Solve system using a direct solver (KLU).
   *
   * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
   */
  void dsolve( EpetraMatrix &Mat,
               EpetraVector &sol,
               EpetraVector &rhs );
  //@}

protected:

};

}

#endif /* TRILINOSSOLVER_HPP_ */
