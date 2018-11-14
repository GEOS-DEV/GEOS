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
 * @file InterfaceTypes.hpp
 *
 *  Created on: Sep 14, 2018
 *      Author: settgast1
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_INTERFACETYPES_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_INTERFACETYPES_HPP_

namespace geosx
{

namespace trilinosTypes
{
using gid = int;
using lid = int;
}

namespace hypreTypes
{
using gid = long long;
using lid = int;
}

namespace petscTypes
{
using gid = long long;
using lid = int;
}


}



#endif /* SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_INTERFACETYPES_HPP_ */
