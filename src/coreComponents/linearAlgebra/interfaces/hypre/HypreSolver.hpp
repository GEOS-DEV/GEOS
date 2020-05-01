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
 * @file HypreSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_HYPRESOLVER_HPP_
#define GEOSX_LINEARALGEBRA_HYPRESOLVER_HPP_

namespace geosx
{

class HypreVector;
class HypreMatrix;
class LinearSolverParameters;

/**
 * \class TrilinosSolver
 * \brief This class creates and provides basic support for ... .
 */

class HypreSolver
{
public:

  /**
   * @brief Solver constructor, with parameter list reference
   *
   */
  HypreSolver( LinearSolverParameters const & parameters );

  /**
   * @brief Virtual destructor.
   *
   */
  virtual ~HypreSolver() = default;

  /**
   * @brief Solve system with an iterative solver (HARD CODED PARAMETERS, GMRES).
   *
   * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
   */

  void solve( HypreMatrix & mat,
              HypreVector & sol,
              HypreVector & rhs );


private:

  LinearSolverParameters const & m_parameters;

  void solve_direct( HypreMatrix & mat,
                     HypreVector & sol,
                     HypreVector & rhs );

  void solve_krylov( HypreMatrix & mat,
                     HypreVector & sol,
                     HypreVector & rhs );

};

} // end geosx namespace

#endif /* HYPRESOLVER_HPP_ */
