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
 * @file TrilinosInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_

#include "linearAlgebra/interfaces/trilinos/EpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "TrilinosSolver.hpp"

namespace geosx
{

/**
 * @class TrilinosInterface
 * @brief This class holds aliases based on the Trilinos library.
 */
struct TrilinosInterface
{
  /**
   * @brief Initializes the MPI environment for the Trilinos library
   *
   * @param[in] argc standard argc as in any C main
   * @param[in] argv standard argv as in any C main
   */
  static void initialize( int & argc, char * * & argv )
  {
    GEOSX_UNUSED_VAR( argc );
    GEOSX_UNUSED_VAR( argv );
  }

  /**
   * @brief Finalizes the MPI environment for the Trilinos library
   */
  static void finalize() {}

  /// Alias for EpetraMatrix
  using ParallelMatrix = EpetraMatrix;
  /// Alias for EpetraVector
  using ParallelVector = EpetraVector;
  /// Alias for TrilinosSolver
  using LinearSolver   = TrilinosSolver;
};

} /* namespace geosx */

#endif /* GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_ */
