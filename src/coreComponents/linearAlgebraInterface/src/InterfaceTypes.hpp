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
 * @file InterfaceTypes.hpp
 *
 *  Created on: Sep 14, 2018
 *      Author: settgast1
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_INTERFACETYPES_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_INTERFACETYPES_HPP_

#include "common/GeosxConfig.hpp"

#include "TrilinosInterface.hpp"

#define CONCAT_( A, B ) A##B
#define CONCAT( A, B ) CONCAT_( A, B )

namespace geosx
{

using LAInterface = CONCAT( GEOSX_LA_INTERFACE, Interface );

// The following aliases are added into geosx namespace for global use
using ParallelMatrix = LAInterface::ParallelMatrix;
using ParallelVector = LAInterface::ParallelVector;
using LinearSolver   = LAInterface::LinearSolver;

}



#endif /* SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_INTERFACETYPES_HPP_ */
