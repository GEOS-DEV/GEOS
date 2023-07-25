/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DataTypes.cpp
 */


#include "DataTypes.hpp"
#include "Logger.hpp"
#include "LvArray/src/system.hpp"

namespace geos
{
#ifdef GEOSX_USE_MPI
MPI_Comm MPI_COMM_GEOSX;
#else
int MPI_COMM_GEOSX = 0;
#endif

void printTypeSummary()
{
  logger.rank0Log( "real64 is alias of ", LvArray::system::demangle( typeid(real64).name() ) );
  logger.rank0Log( "localIndex is alias of ", LvArray::system::demangle( typeid(localIndex).name() ) );
  logger.rank0Log( "globalIndex is alias of ", LvArray::system::demangle( typeid(globalIndex).name()) );
}


}
