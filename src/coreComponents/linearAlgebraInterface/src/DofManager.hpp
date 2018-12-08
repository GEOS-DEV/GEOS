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
 * @file DofManager.hpp
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

class DofManager
{
public:

  /**
   * Constructor.
   */
  DofManager() = default;

  /**
   * Destructor.
   */
  ~DofManager() = default;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_ */

