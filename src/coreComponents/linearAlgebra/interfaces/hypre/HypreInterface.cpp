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
 * @file HypreInterface.cpp
 */

#include "HypreInterface.hpp"

#include "common/GeosxConfig.hpp"
#include "linearAlgebra/interfaces/direct/SuiteSparse.hpp"
#ifdef GEOS_USE_HYPREDRV
#include "linearAlgebra/interfaces/hypre/hypredrive.hpp"
#endif
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HyprePreconditioner.hpp"
#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#if defined(GEOS_USE_SUPERLU_DIST)
#include "linearAlgebra/interfaces/direct/SuperLUDist.hpp"
#endif

#include "HYPRE_utilities.h"
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
#include "_hypre_utilities.h"
#include "_hypre_utilities.hpp"
#endif

namespace geos
{

void HypreInterface::initialize()
{
#if defined(GEOS_USE_OPENMP) && defined(HYPRE_USING_OPENMP)
  GEOS_LOG_RANK_0_IF( omp_get_max_threads() > 1,
                      "\n"
                      "********************************************************************\n"
                      "*                                                                  *\n"
                      "*    WARNING: OMP_NUM_THREADS > 1 MAY NOT BE OPTIMAL FOR CERTAIN   *\n"
                      "*             HYPRE PRECONDITIONING OPTIONS!                       *\n"
                      "*                                                                  *\n"
                      "********************************************************************\n"
                      );
#endif

#ifdef GEOS_USE_HYPREDRV
  hypre::hypreDrive::initializeRuntime();
#else
  HYPRE_Initialize();
#endif
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  HYPRE_SetExecutionPolicy( HYPRE_EXEC_DEVICE );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA
  HYPRE_SetSpGemmUseVendor( 0 );
#else
  HYPRE_SetSpGemmUseVendor( 1 );
#endif
#if !defined(GEOS_USE_HYPREDRV) || !HYPRE_CHECK_MIN_VERSION(23100, 0)
  HYPRE_DeviceInitialize();
#endif
#endif
  HYPRE_SetMemoryLocation( hypre::memoryLocation );
  HYPRE_SetPrintErrorMode( 1 );

#if defined(HYPRE_USING_UMPIRE)
  HYPRE_SetUmpireUMPoolName( "HYPRE_UM" );
  HYPRE_SetUmpireHostPoolName( "HYPRE_HOST" );
  HYPRE_SetUmpireDevicePoolName( "HYPRE_DEVICE" );
  HYPRE_SetUmpirePinnedPoolName( "HYPRE_PINNED" );
#endif

  HYPRE_SetLogLevel( getenv( "HYPRE_LOG_LEVEL" ) ? atoi( getenv( "HYPRE_LOG_LEVEL" ) ) : 0 );
}

void HypreInterface::finalize()
{
#ifdef GEOS_USE_HYPREDRV
  hypre::hypreDrive::finalizeRuntime();
#else
  HYPRE_Finalize();
#endif
}

std::unique_ptr< LinearSolverBase< HypreInterface > >
HypreInterface::createSolver( LinearSolverParameters params )
{
#ifdef GEOS_USE_HYPREDRV
  if( hypre::hypreDrive::shouldUse( params ) )
  {
    return std::make_unique< HypreDriveSolver >( std::move( params ) );
  }
#endif

  if( params.solverType == LinearSolverParameters::SolverType::direct )
  {
    if( params.direct.parallel )
    {
#if defined(GEOS_USE_SUPERLU_DIST)
      return std::make_unique< SuperLUDist< HypreInterface > >( std::move( params ) );
#else
      GEOS_ERROR( "GEOS is configured without support for SuperLU_dist." );
      return std::unique_ptr< LinearSolverBase< HypreInterface > >( nullptr );
#endif
    }
    else
    {
      return std::make_unique< SuiteSparse< HypreInterface > >( std::move( params ) );
    }
  }
  else
  {
    return std::make_unique< HypreSolver >( std::move( params ) );
  }
}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geos::HypreInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< HyprePreconditioner >( std::move( params ) );
}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geos::HypreInterface::createPreconditioner( LinearSolverParameters params,
                                            arrayView1d< HypreVector const > nearNullKernel )
{
  return std::make_unique< HyprePreconditioner >( std::move( params ), nearNullKernel );
}

}
