/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/**
 * @brief Setup PETSc.
 * @param [in] argc the number of command line arguments.
 * @param [in/out] argv the command line arguments.
 */
void setupPetsc( int argc, char * argv[] );

/**
 * @brief Finalize PETSc.
 */
void finalizePetsc();

} // namespace geosx

#endif //_COMMON_INITIALIZATION_HPP_
