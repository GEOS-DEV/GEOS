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
 * @file HypreInterface.hpp
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_

//#include "HypreSparseMatrix.hpp"
//#include "HypreVector.hpp"
//#include "HypreSolver.hpp"

namespace geosx
{

/**
 * \class HypreInterface
 * \brief This class holds aliases based on the Hypre library.
 */

class HypreInterface
{
public:

  using laiLID = hypreTypes::lid;
  using laiGID = hypreTypes::gid;

  // Epetra matrix and vector wrappers
  using ParallelMatrix = HypreSparseMatrix;
  using ParallelVector = HypreVector;

  using LinearSolver = HypreSolver;

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   */
  HypreInterface() = default;
  /**
   * @brief Destructor.
   */
  ~HypreInterface() = default;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_ */
