/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
