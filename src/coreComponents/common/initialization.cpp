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

#include "initialization.hpp"
#include "common/DataTypes.hpp"
#include "SetFPE.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include <petscvec.h>


#ifdef GEOSX_USE_MKL
#include <mkl.h>
#endif

#ifdef GEOSX_USE_OPENMP
#include <omp.h>
#endif

namespace geosx
{

void setupMKL()
{
#ifdef GEOSX_USE_MKL
  GEOS_LOG_RANK_0( "MKL max threads: " << mkl_get_max_threads());
#endif
}

void setupOpenMP()
{
#ifdef GEOSX_USE_OPENMP
  GEOS_LOG_RANK_0( "Max threads: " << omp_get_max_threads());
#endif
}

void setupMPI( int argc, char * argv[] )
{
#ifdef GEOSX_USE_MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
#endif
}

void finalizeMPI()
{
#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif
}

void setupCXXUtils()
{
#ifdef GEOSX_USE_MPI
  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger() :
#endif

  cxx_utilities::setSignalHandling( cxx_utilities::handler1 );
  cxx_utilities::SetFPE();
}

void finalizeCXXUtils()
{
#ifdef GEOSX_USE_CHAI
  chai::ArrayManager::finalize();
#endif

  logger::FinalizeLogger();
}

void basicSetup( int argc, char * argv[] )
{
  setupMPI( argc, argv );
  setupOpenMP();
  setupMKL();
  setupCXXUtils();
}

void basicCleanup()
{
  FieldSpecificationManager::finalize();
  NewFunctionManager::finalize();

  finalizeCXXUtils();
  finalizeMPI();
}

void setupPetsc( int argc, char * argv[] )
{
  char help[] = "Setup PETSc.\n";
  PetscOptionsSetValue( nullptr, "-log_view", "" );
  PetscOptionsSetValue( nullptr, "-ksp_monitor", nullptr );
  PetscInitialize( &argc, &argv, nullptr, help );
}

void finalizePetsc()
{
  PetscFinalize();
}

} // namespace geosx
