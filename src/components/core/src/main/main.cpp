
#if USE_CALIPER==1
#include "caliper/Annotation.h"
#endif
#include <slic/slic.hpp>
#include <slic/GenericOutputStream.hpp>
#include <iostream>
#include <sys/time.h>
//#include "ManagedArray.hpp"
#include "../dataRepository/ManagedGroup.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"
#include <mpi.h>

using namespace geosx;
using namespace axom;


int main( int argc, char *argv[] )
{
  timeval tim;
  gettimeofday(&tim, NULL);
  real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);
  real64 t_initialize, t_run;

#if USE_MPI

  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::cout<<"rank = "<<rank<<std::endl;
#endif


  std::cout<<"starting main"<<std::endl;


  slic::initialize();


  std::string format =  std::string( "***********************************\n" )+
                        std::string( "* <TIMESTAMP>\n\n" ) +
                        std::string( "* LEVEL=<LEVEL>\n" ) +
                        std::string( "* MESSAGE=<MESSAGE>\n" ) +
                        std::string( "* FILE=<FILE>\n" ) +
                        std::string( "* LINE=<LINE>\n" ) +
                        std::string( "***********************************\n" );
  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::addStreamToAllMsgLevels( new slic::GenericOutputStream( &std::cout, format ) );

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);



  // Mark begin of "initialization" phase
#if USE_CALIPER==1
  cali::Annotation init_ann = cali::Annotation("initialization").begin();
#endif

  ProblemManager problemManager( "ProblemManager", nullptr );

  problemManager.BuildDataStructure(nullptr);
  problemManager.SetDocumentationNodes( &problemManager );

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );
  problemManager.ParseInputFile();


  
  problemManager.Initialize( &problemManager );
#if USE_CALIPER==1
  init_ann.end();
#endif
  gettimeofday(&tim, NULL);
  t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

  problemManager.ApplyInitialConditions();
  std::cout << std::endl << "Running simulation:" << std::endl;
  problemManager.RunSimulation();
  

  problemManager.ClosePythonInterpreter();

  slic::finalize();

  gettimeofday(&tim, NULL);
  t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

  printf("Done!\n\nScaling Data: initTime = %1.2fs, runTime = %1.2fs\n", t_initialize - t_start,  t_run - t_initialize );
  
  return 0;
}
