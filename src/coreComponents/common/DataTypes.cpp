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
 * @file DataTypes.cpp
 */


#include "DataTypes.hpp"
#include "common/Logger.hpp"
#include "cxx-utilities/src/StringUtilities.hpp"

namespace geosx
{
#ifdef GEOSX_USE_MPI
MPI_Comm MPI_COMM_GEOSX;
#endif

void printTypeSummary()
{
  GEOSX_LOG_RANK_0( "real64 is alias of " <<cxx_utilities::demangle( typeid(real64).name() ) );
  GEOSX_LOG_RANK_0( "localIndex is alias of " <<cxx_utilities::demangle( typeid(localIndex).name() ) );
  GEOSX_LOG_RANK_0( "globalIndex is alias of "<<cxx_utilities::demangle( typeid(globalIndex).name()) );
}


}
