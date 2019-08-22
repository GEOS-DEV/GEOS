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

/*
 * DataTypes.cpp
 *
 *  Created on: Mar 23, 2017
 *      Author: rrsettgast
 */


#include "DataTypes.hpp"
#include "Logger.hpp"
#include "StringUtilities.hpp"

namespace geosx
{
#ifdef GEOSX_USE_MPI
MPI_Comm MPI_COMM_GEOSX;
#endif

void printTypeSummary()
{
  GEOS_LOG_RANK_0( "real64 is alias of " <<cxx_utilities::demangle( typeid(real64).name() ) );
  GEOS_LOG_RANK_0( "localIndex is alias of " <<cxx_utilities::demangle( typeid(localIndex).name() ) );
  GEOS_LOG_RANK_0( "globalIndex is alias of "<<cxx_utilities::demangle( typeid(globalIndex).name()) );
}


}
