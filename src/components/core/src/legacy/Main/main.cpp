//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  void sigproc(int s); //catch function
}


/// returns the amount of cpu time use for this process
realT getcputime(void)  ;

#ifdef SRC_INTERNAL
/// Coupling to other codes Comm
extern MPI_Comm MPI_DOMAIN_COMM ;
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
//  _MM_SET_EXCEPTION_MASK(_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO | _MM_EXCEPT_OVERFLOW);
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
      //  printf( "Scaling Data: mpiSize, rank, runtime:  %12.3f %12.3f %12.3f \n", size, rank, t_run-t_initialize );

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
