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
 */
#include "common/DataTypes.hpp"

namespace geos
{
/**
 * @brief Get GEOSX version.
 * @return The full version string.
 */
string getVersion();

/**
 * @brief output version info for dependencies to log
 */
void outputVersionInfo();

/**
 * @brief Get a ID and version of the c++ compiler
 *
 * @return A string containing the c++ compiler ID and version
 */
string getCppCompilerIdString();

/**
 * @brief Get the Gpu Compiler type and version
 *
 * @return A string containing the GPU compiler ID and version
 */
string getGpuCompilerIdString();

}
