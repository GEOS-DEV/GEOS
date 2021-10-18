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

#ifndef GEOSX_MAININTERFACE_INITIALIZATION_HPP_
#define GEOSX_MAININTERFACE_INITIALIZATION_HPP_

// Source includes
#include "common/initializeEnvironment.hpp"

namespace geosx
{

/**
 * @brief Parse the command line options and populate @p commandLineOptions with the results.
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return The command line options.
 */
std::unique_ptr< CommandLineOptions > parseCommandLineOptions( int argc, char * * argv );

/**
 * @brief Perform the basic GEOSX initialization and optionally parse the command line input.
 * @param [in] argc The number of command line arguments.
 * @param [in,out] argv The command line arguments.
 * @param [in] parseCommandLine True iff the command line options should be parsed.
 * @return The command line options, if @c parseCommandLine is @c false then the returned value
 *   is default constructed (empty).
 */
std::unique_ptr< CommandLineOptions > basicSetup( int argc, char * argv[], bool const parseCommandLine=false );

/**
 * @brief Perform the basic GEOSX cleanup.
 */
void basicCleanup();

/**
 * @brief Get GEOSX version.
 * @return The full version string.
 */
string getVersion();

} // namespace geosx

#endif // GEOSX_MAININTERFACE_INITIALIZATION_HPP_
