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
