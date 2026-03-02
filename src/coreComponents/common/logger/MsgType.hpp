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
 * @file MsgType.hpp
 */

#ifndef INITIALIZATION_MSG_TYPE_HPP
#define INITIALIZATION_MSG_TYPE_HPP

#include "common/format/EnumStrings.hpp"

namespace geos
{

/**
 * @enum MsgType
 * Enum listing the different types of possible errors
 */
enum class MsgType
{
  Error,
  ExternalError,
  Warning,
  Exception,
  Undefined
};

/// Declare strings associated with output MsgType values.
ENUM_STRINGS( MsgType,
              "Error",
              "ExternalError",
              "Warning",
              "Exception",
              "Undefined" );

}

#endif
