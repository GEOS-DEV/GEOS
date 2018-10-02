/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file TrilinosInterface.hpp
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_

#include "TrilinosSolver.hpp"
#include "TrilinosSparseMatrix.hpp"
#include "TrilinosVector.hpp"

namespace geosx
{

/**
 * \class TrilinosInterface
 * \brief This class holds aliases based on the Trilinos library.
 */

class TrilinosInterface
{
public:

  using laiLID = trilinosTypes::lid;
  using laiGID = trilinosTypes::gid;

  // Epetra matrix and vector wrappers
  using ParallelMatrix = EpetraSparseMatrix;
  using ParallelVector = EpetraVector;

  using LinearSolver = TrilinosSolver;

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   */
  TrilinosInterface() = default;

  /**
   * @brief Destructor.
   *
   */
  ~TrilinosInterface() = default;
  //@}

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_ */
