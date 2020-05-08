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
 * @file TrilinosTpetraInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRAINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRAINTERFACE_HPP_

#include "linearAlgebra/interfaces/trilinos/TpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/TpetraMatrix.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosTpetraSolver.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geosx
{

/**
 * @brief This class holds aliases based on the Trilinos library.
 */
struct TrilinosTpetraInterface
{
/**
 * @brief Initializes the Trilinos library
 *
 * @param[in] argc standard argc as in any C main
 * @param[in] argv standard argv as in any C main
 */
  static void initialize( int & argc, char * * & argv );

  /**
   * @brief Finalizes the Trilinos library
   */
  static void finalize();

  /**
   * @brief Create a Trilinos-based preconditioner object.
   * @param params the preconditioner parameters
   * @return an owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< TrilinosTpetraInterface > >
  createPreconditioner( LinearSolverParameters params );

  /// Alias for TpetraMatrix
  using ParallelMatrix = TpetraMatrix;
  /// Alias for TpetraVector
  using ParallelVector = TpetraVector;
  /// Alias for TrilinosTpetraSolver
  using LinearSolver   = TrilinosTpetraSolver;
};

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRAINTERFACE_HPP_
