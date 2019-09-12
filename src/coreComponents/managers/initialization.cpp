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

#include "initialization.hpp"

#include "common/DataTypes.hpp"
#include "SetFPE.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

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

void setupMPI( int MPI_PARAM(argc), char * MPI_PARAM(argv[]) )
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
  logger::InitializeLogger();
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
  setupLAI( argc, argv );
}

void basicCleanup()
{
  FieldSpecificationManager::finalize();
  NewFunctionManager::finalize();

  finalizeLAI();
  finalizeCXXUtils();
  finalizeMPI();
}

} // namespace geosx
