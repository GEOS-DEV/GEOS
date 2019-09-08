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

#ifndef _COMMON_INITIALIZATION_HPP_
#define _COMMON_INITIALIZATION_HPP_

namespace geosx
{

/**
 * @brief Perform the basic GEOSX initialization.
 * @param [in] argc the number of command line arguments.
 * @param [in/out] argv the command line arguments.
 */
void basicSetup( int argc, char * argv[] );

/**
 * @brief Perform the basic GEOSX cleanup.
 */
void basicCleanup();

/**
 * @brief Setup the cxx-utilities library. This initializes the logger,
 * signal handling and the floating point environment.
 */
void setupCXXUtils();

/**
 * @brief Finalize the cxx-utilities library. This finalizes the logger and chai.
 */
void finalizeCXXUtils();

/**
 * @brief We link to MKL with the single dynamic library approach, so this sets some MKL parameters.
 */
void setupMKL();

/**
 * @brief Setup OpenMP.
 */
void setupOpenMP();

/**
 * @brief Setup MPI.
 * @param [in] argc the number of command line arguments.
 * @param [in/out] argv the command line arguments.
 */
void setupMPI( int argc, char * argv[] );

/**
 * @brief Finalize MPI.
 */
void finalizeMPI();

} // namespace geosx

#endif //_COMMON_INITIALIZATION_HPP_
