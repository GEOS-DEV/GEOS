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
 * @file HypreInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREINTERFACE_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/common/LinearSolverBase.hpp"

#include <memory>

namespace geosx
{

/**
 * @class HypreInterface
 * @brief This class holds aliases based on the Hypre library.
 */
struct HypreInterface
{
  /**
   * @brief Initializes the Hypre library.
   */
  static void initialize();

  /**
   * @brief Finalizes the Hypre library.
   */
  static void finalize();

  /**
   * @brief Create a hypre-based solver object.
   * @param params the preconditioner parameters
   * @return owning pointer to the newly created solver
   */
  static std::unique_ptr< LinearSolverBase< HypreInterface > >
  createSolver( LinearSolverParameters params );

  /**
   * @brief Create a hypre-based preconditioner object.
   * @param params the preconditioner parameters
   * @return owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< HypreInterface > >
  createPreconditioner( LinearSolverParameters params );

  /**
   * @brief Create a hypre-based preconditioner object.
   * @param params the preconditioner parameters
   * @param nearNullKernel the user-provided near null kernel
   * @return owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< HypreInterface > >
  createPreconditioner( LinearSolverParameters params,
                        array1d< HypreVector > const & nearNullKernel );

  /// Alias for HypreMatrix
  using ParallelMatrix = HypreMatrix;

  /// Alias for HypreVector
  using ParallelVector = HypreVector;

};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREINTERFACE_HPP_*/
