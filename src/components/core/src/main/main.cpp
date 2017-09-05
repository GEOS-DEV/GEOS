
#if USE_CALIPER
#include "caliper/Annotation.h"
#endif

#include "common/Logger.hpp"

#include <mpi.h>
#include <iostream>
#include <sys/time.h>
#include "../dataRepository/ManagedGroup.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"

using namespace geosx;


int main( int argc, char *argv[] )
{
  timeval tim;
  gettimeofday(&tim, NULL);
  real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);
  real64 t_initialize, t_run;

#if USE_MPI == 1
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  std::cout<<"starting main"<<std::endl;

#if ATK_FOUND
  axom::slic::initialize();
  std::string format =  std::string( "***********************************\n" )+
                        std::string( "* <TIMESTAMP>\n\n" ) +
                        std::string( "* LEVEL=<LEVEL>\n" ) +
                        std::string( "* MESSAGE=<MESSAGE>\n" ) +
                        std::string( "* FILE=<FILE>\n" ) +
                        std::string( "* LINE=<LINE>\n" ) +
                        std::string( "***********************************\n" );
  axom::slic::setLoggingMsgLevel( axom::slic::message::Debug );
  axom::slic::addStreamToAllMsgLevels( new axom::slic::GenericOutputStream( &std::cout, format ) );

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);
#endif



  // Mark begin of "initialization" phase
#if USE_CALIPER
  cali::Annotation init_ann = cali::Annotation("initialization").begin();
#endif

  ProblemManager problemManager( "ProblemManager", nullptr );

  problemManager.BuildDataStructure(nullptr);
  problemManager.SetDocumentationNodes( &problemManager );

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );
  problemManager.ParseInputFile();


  
  problemManager.Initialize( &problemManager );
#if USE_CALIPER
  init_ann.end();
#endif
  gettimeofday(&tim, NULL);
  t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

  problemManager.ApplyInitialConditions();
  std::cout << std::endl << "Running simulation:" << std::endl;

#if USE_CALIPER
  cali::Annotation simulation_ann = cali::Annotation("RunSimulation").begin();
#endif

  problemManager.RunSimulation();

#if USE_CALIPER
  simulation_ann.end();
#endif

  problemManager.ClosePythonInterpreter();

#if ATK_FOUND
  axom::slic::finalize();
#endif
  
  gettimeofday(&tim, NULL);
  t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

  printf("Done!\n\nScaling Data: initTime = %1.2fs, runTime = %1.2fs\n", t_initialize - t_start,  t_run - t_initialize );

#if USE_MPI
  MPI_Finalize();
#endif

  
  return 0;
}
