/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TrilinosInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_

#include "linearAlgebra/interfaces/trilinos/EpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/common/LinearSolverBase.hpp"

#include <memory>

namespace geosx
{

/**
 * @class TrilinosInterface
 * @brief This class holds aliases based on the Trilinos library.
 */
struct TrilinosInterface
{
  /**
   * @brief Initializes the Trilinos library.
   */
  static void initialize();

  /**
   * @brief Finalizes the Trilinos library.
   */
  static void finalize();

  /**
   * @brief Create a petsc-based solver object.
   * @param params the preconditioner parameters
   * @return owning pointer to the newly created solver
   */
  static std::unique_ptr< LinearSolverBase< TrilinosInterface > >
  createSolver( LinearSolverParameters params );

  /**
   * @brief Create a Trilinos-based preconditioner object.
   * @param params the preconditioner parameters
   * @return an owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< TrilinosInterface > >
  createPreconditioner( LinearSolverParameters params );

  /**
   * @brief Create a Trilinos-based preconditioner object.
   * @param params the preconditioner parameters
   * @param nearNullKernel the user-provided near null kernel
   * @return an owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< TrilinosInterface > >
  createPreconditioner( LinearSolverParameters params,
                        array1d< EpetraVector > const & nearNullKernel );

  /// Alias for EpetraMatrix
  using ParallelMatrix = EpetraMatrix;
  /// Alias for EpetraVector
  using ParallelVector = EpetraVector;
};

} /* namespace geosx */

#endif /* GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_ */
