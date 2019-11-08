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
#include "cxx-utilities/src/SetFPE.hpp"
#include "cxx-utilities/src/SetSignalHandling.hpp"
#include "cxx-utilities/src/stackTrace.hpp"
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
  GEOSX_LOG_RANK_0( "MKL max threads: " << mkl_get_max_threads() );
#endif
}

void setupOpenMP()
{
#ifdef GEOSX_USE_OPENMP
  GEOSX_LOG_RANK_0( "Max threads: " << omp_get_max_threads() );
#endif
}

void setupMPI( int argc, char * argv[] )
{
  MpiWrapper::Init( &argc, &argv );
#ifdef GEOSX_USE_MPI
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

void setupLogger()
{
#ifdef GEOSX_USE_MPI
  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger();
#endif
}

void setupCXXUtils()
{
  cxx_utilities::setSignalHandling( cxx_utilities::handler1 );
  cxx_utilities::SetFPE();
}

void finalizeLogger()
{
  logger::FinalizeLogger();
}

void basicSetup( int argc, char * argv[] )
{
  setupMPI( argc, argv );
  setupLogger();
  setupCXXUtils();
  setupOpenMP();
  setupMKL();
  setupLAI( argc, argv );
}

void basicCleanup()
{
  finalizeLAI();
  finalizeLogger();
  finalizeMPI();
}

} // namespace geosx
