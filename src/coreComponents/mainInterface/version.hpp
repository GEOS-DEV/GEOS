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
 * @file version.hpp
 *
 * Utilities for getting and printing version info of GEOS,
 * its build environment and dependencies.
 */

#include <string>

namespace geos
{
/**
 * @brief Get GEOSX version.
 * @return The full version string.
 */
std::string getVersion();

/**
 * @brief output version info for dependencies to log
 */
void outputVersionInfo();

}
