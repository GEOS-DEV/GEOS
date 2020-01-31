/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TrilinosSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSSOLVER_HPP_

// Forward declaration of Amesos_BaseSolver.
class Amesos_BaseSolver;

namespace geosx
{

class EpetraVector;
class EpetraMatrix;
class LinearSolverParameters;

/**
 * \class TrilinosSolver
 * \brief This class creates and provides basic support for AztecOO, Amesos and ML libraries.
 */

class TrilinosSolver
{
public:

  /**
   * @brief Solver constructor, with parameter list reference
   *
   */
  TrilinosSolver( LinearSolverParameters const & parameters );

  /**
   * @brief Virtual destructor.
   *
   */
  ~TrilinosSolver();

  /**
   * @brief Solve system with an iterative solver.
   *
   * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
   */

  void solve( EpetraMatrix & mat,
              EpetraVector & sol,
              EpetraVector & rhs );

private:

  LinearSolverParameters const & m_parameters;

  Amesos_BaseSolver * m_solver = nullptr;

  void solve_direct( EpetraMatrix & mat,
                     EpetraVector & sol,
                     EpetraVector & rhs );

  void solve_krylov( EpetraMatrix & mat,
                     EpetraVector & sol,
                     EpetraVector & rhs );

};

} // end geosx namespace

#endif /* TRILINOSSOLVER_HPP_ */
