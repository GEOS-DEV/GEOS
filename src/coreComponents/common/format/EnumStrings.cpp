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
 * @file EnumStrings.cpp
 */

 #include "common/format/EnumStrings.hpp"
 #include "common/logger/Logger.hpp"

namespace geos
{

void internal::EnumErrorMessageToString( size_t index,
                                         string const & typeName,
                                         std::size_t size )
{
  GEOS_THROW( "Invalid value " << index << " of type " << typeName<< ". Valid range is 0.." << size - 1,
              InputError );
}

void internal::EnumErrorMessageFromString( string const & s,
                                           string const & typeName,
                                           string const & concat )
{
  GEOS_THROW( "Invalid value '" << s << "' of type " << typeName << ". Valid options are: " << concat,
              InputError );
}

};
