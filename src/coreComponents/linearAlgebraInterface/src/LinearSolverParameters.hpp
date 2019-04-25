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
 * @file LinearSolverParameters.hpp
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LINEARSOLVERPARAMETERS_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LINEARSOLVERPARAMETERS_HPP_

// #include "../../common/DataTypes.hpp"

namespace geosx
{

/**
 * \class LinearSolverParameters
 * This class holds a simple tree of linear solver options.  They are set to
 * default values, but can be overwritten as needed.
 */

class LinearSolverParameters
{
public:

  integer verbosity = 0;               //!< Output level [0=none, 1=basic, 2=everything]
  string  solverType = "cg";           //!< Solver type [direct, cg, gmres, bicgstab]
  string  preconditionerType = "ilut"; //!< Preconditioner type [none, ilu, ilut, icc, amg]
  integer dofsPerNode = 1;             //!< Can be used to enable dense-block algorithms if available

  struct
  {
    real64  tolerance = 1e-6;
    integer maxIterations = 200;
    integer maxRestart = 200;
  }
  krylov;

  struct
  {
    bool useRowScaling = false;
    bool useRowColScaling = false;
  }
  scaling; //TODO: not implemented

  struct
  {
    integer maxLevels = 20;
    string  cycleType = "V";
    string  smootherType = "gaussSeidel";
    string  coarseType = "direct";
    integer numSweeps = 2;
    bool    isSymmetric = true;
    string  nullSpaceType = "constantModes";
    solver = "Petsc";
  }
  amg;
  
  struct
  {
    integer fill = 0;
    real64  threshold = 0.0;
  }
  ilu;

  /**
   * @brief Constructor.
   */
  LinearSolverParameters() = default;

  /**
   * @brief Destructor.
   *
   */
  ~LinearSolverParameters() = default;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_LINEARSOLVERPARAMETERS_HPP_ */
