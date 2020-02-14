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

#include "TrilinosSolver.hpp"
#include "EpetraMatrix.hpp"
#include "EpetraVector.hpp"

namespace geosx
{

/**
 * \class TrilinosInterface
 * \brief This class holds aliases based on the Trilinos library.
 */

struct TrilinosInterface
{
  static void initialize( int & GEOSX_UNUSED_PARAM( argc ), char ** & GEOSX_UNUSED_PARAM( argv ) ) {}

  static void finalize() {}

  using ParallelMatrix = EpetraMatrix;
  using ParallelVector = EpetraVector;
  using LinearSolver   = TrilinosSolver;
};

} /* namespace geosx */

#endif /* GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSINTERFACE_HPP_ */
