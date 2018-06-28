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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @brief This file contains the "main" function for the GEOS code.
 * @file main.cpp
 * @author Randolph Settgast
 */

#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <signal.h>
#include "ObjectManagers/ProblemManagerT.h"
#include "Utilities/Utilities.h"
#include "fenv.h"

#ifndef __xlC__
#include <xmmintrin.h>
#endif

#if GPAC_MPI
#include <mpi.h>
#endif


int term_flag = 0;
extern "C"
{
void sigproc(int s);   //catch function
}


/// returns the amount of cpu time use for this process
realT getcputime(void);

#ifdef SRC_INTERNAL
/// Coupling to other codes Comm
extern MPI_Comm MPI_DOMAIN_COMM;
#endif

#include "DataStructures/EncapsulatedObjects/EncapsulatedObjectBase.h"

/**
 * @author Randolph Settgast
 * @param argc number of command line arguments including executable
 * @param argv array of command line arguments
 * @return 0
 */
int main( int argc, char* argv[] )
{

  timeval tim;
  gettimeofday(&tim, NULL);
  realT t_start=tim.tv_sec+(tim.tv_usec/1000000.0);
  realT t_initialize, t_run;

  signal(SIGINT, &sigproc);   //set handler for interrupt
  signal(SIGTERM, &sigproc);   //set handler for interrupt


  int rank = 0, size = 1, returnValue = 0;

#if GPAC_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#if SRC_INTERNAL
  MPI_DOMAIN_COMM=MPI_COMM_WORLD;
#endif
#endif


#ifdef  __IBMC__

#endif

#ifdef __APPLE__// && __MACH__
//  _MM_SET_EXCEPTION_MASK(  _MM_EXCEPT_DIV_ZERO);
//  _MM_SET_EXCEPTION_MASK(_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO |
// _MM_EXCEPT_OVERFLOW);
//  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
//fesetenv(FE_ALL_EXCEPT);
#endif

#ifdef __INTEL_COMPILER
  feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif

  {
    //should destruct all before MPI finish
    ProblemManagerT ProblemManager;
    try
    {
      ProblemManager.ParseCommandLineInput(argc, argv);
      ProblemManager.DisplaySplash();
      ProblemManager.ProblemSetup();
      ProblemManager.VerifySolvers();
      std::remove("geos.visit");


      MPI_Barrier(MPI_COMM_WORLD);
      realT cputime0 = getcputime();

      gettimeofday(&tim, NULL);
      t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

      double outputTime = 0.0;

      ProblemManager.Run(outputTime, t_start);

      gettimeofday(&tim, NULL);
      t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

      realT cputime1 = getcputime();
      realT cputime = cputime1 - cputime0;
      realT t_total = t_run - t_start;
      if (rank == 0)
      {
        printf(
          "   RANK   SETUP TIME     RUN TIME   TOTAL TIME     CPU TIME  OUTPUT TIME    WAIT TIME       %%WAIT\n");
        printf(
          " ------ ------------ ------------ ------------ ------------ ------------ ------------ ------------\n");
      }

      MPI_Barrier(MPI_COMM_WORLD);
      printf(" %6i %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f \n", rank,
             t_initialize - t_start, t_run - t_initialize, t_total, cputime, outputTime,
             t_total - cputime, (t_total - cputime) / cputime * 100);
      //  printf( "Scaling Data: mpiSize, rank, runtime:  %12.3f %12.3f %12.3f
      // \n", size, rank, t_run-t_initialize );

      realT runtime = t_run - t_initialize - outputTime;
      realT maxruntime;
      MPI_Reduce(&runtime, &maxruntime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (rank == 0)
        printf("Scaling Data: mpiSize, numCycles, runtime:  %6i %6i %12.3f \n", size,
               ProblemManager.m_cycleNumber, maxruntime);
    }
    catch (GPException& e)
    {
      std::cout << "GP Exception: " << e.what() << std::endl;
      returnValue = 1;
      ProblemManager.WriteSilo(false);
    }
    catch (std::exception& e)
    {
      std::cout << "Standard Exception: " << e.what() << std::endl;
      returnValue = 1;
      ProblemManager.WriteSilo(false);
    }
  }
#if GPAC_MPI
  MPI_Finalize();
#endif

  return returnValue;
}

void sigproc(int s )
{
  std::cout << "Received interrupt" << std::endl;
  throw GPException("Received interrupt.\n");
}
