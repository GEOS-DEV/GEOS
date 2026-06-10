/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StdContainerWrappers.cpp
 */

#include "StdContainerWrappers.hpp"

#include "common/logger/Logger.hpp"
#include <sstream>
#include <stdexcept>

namespace geos
{

namespace internal
{


[[noreturn]] void rethrowStdOutOfRange( const char * where,
                                        std::out_of_range const & cause )
{
  std::ostringstream msg;
  msg << "Out of range access in " << where << ": " << cause.what();

  GEOS_THROW( msg.str(), OutOfRangeError );

  // This line should never be reached, but helps the compiler understand the noreturn attribute
  std::terminate();
}


} // namespace geos

} // namespace geos
