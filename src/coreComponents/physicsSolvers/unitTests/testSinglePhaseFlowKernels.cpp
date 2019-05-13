/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "Logger.hpp"

#include "physicsSolvers/FiniteVolume/SinglePhaseFlowKernels.hpp"

using namespace geosx;
using namespace geosx::SinglePhaseFlowKernels;

TEST( testSinglePhaseFlowKernels, mobility )
{
  int constexpr N = 3;

  real64 const dens[N]        = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[N] = { 1e-5,  1e-10,  0.0    };
  real64 const visc[N]        = { 5.0,   2.0,    1.0    };
  real64 const dVisc_dPres[N] = { 1e-7,  0.0,    0.0    };

  for (int i = 0; i < N; ++i)
  {
    SCOPED_TRACE( "Input # " + std::to_string(i) );

    real64 mob;
    real64 dMob_dPres;

    MobilityKernel::Compute( dens[i], dDens_dPres[i], visc[i], dVisc_dPres[i], mob, dMob_dPres );

    real64 const mob_et = dens[i] / visc[i];
    real64 const dMob_dPres_et = mob_et * (dDens_dPres[i] / dens[i] - dVisc_dPres[i] / visc[i]);

    EXPECT_DOUBLE_EQ( mob, mob_et );
    EXPECT_DOUBLE_EQ( dMob_dPres, dMob_dPres_et );
  }
}

int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

#ifdef GEOSX_USE_MPI
  int rank = 0;
  int nranks = 1;

  MPI_Init( &argc, &argv );
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
  MPI_Comm_size( MPI_COMM_GEOSX, &nranks );

  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger():
#endif

  cxx_utilities::setSignalHandling( cxx_utilities::handler1 );

  int const result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}

