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

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCINTERFACE_HPP_

#include "PetscSolver.hpp"
#include "PetscVector.hpp"
#include "PetscSparseMatrix.hpp"

namespace geosx
{

/**
 * \class PetscInterface
 * \brief This class holds aliases based on the Petsc library.
 */

class PetscInterface
{
public:

  // using lid = petscTypes::lid; // no longer necessary
  // using gid = petscTypes::gid; // no longer necessary

  // Petsc matrix and vector wrappers
  using ParallelMatrix = PetscSparseMatrix;
  using ParallelVector = PetscVector;
  using LinearSolver   = PetscSolver;

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   */
  PetscInterface() = default;

  /**
   * @brief Destructor.
   *
   */
  ~PetscInterface() = default;
  //@}

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_PETSCINTERFACE_HPP_ */
