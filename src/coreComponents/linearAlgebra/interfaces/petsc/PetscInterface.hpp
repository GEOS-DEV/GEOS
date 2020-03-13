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
 * @file PetscInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_

#include "PetscSolver.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"

namespace geosx
{

/**
 * \class PetscInterface
 * \brief This class holds aliases based on the Petsc library.
 */

struct PetscInterface
{
  static void initialize( int & argc, char * * & argv );

  static void finalize();

  using ParallelMatrix = PetscMatrix;
  using ParallelVector = PetscVector;
  using LinearSolver   = PetscSolver;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_*/
