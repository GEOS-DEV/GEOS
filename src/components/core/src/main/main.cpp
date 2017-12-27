

#include "common/Logger.hpp"
#include "common/TimingMacros.hpp"

#include <mpi.h>
#include <iostream>
#include <sys/time.h>
#include "../dataRepository/ManagedGroup.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"

using namespace geosx;

#ifdef USE_ATK
using namespace axom;
#endif

int main( int argc, char *argv[] )
{
  timeval tim;
  gettimeofday(&tim, NULL);
  real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);
  real64 t_initialize, t_run;

#ifdef USE_MPI
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  std::cout<<"starting main"<<std::endl;

#ifdef USE_ATK
  slic::initialize();
  std::string format =  std::string( "***********************************\n" )+
                       std::string( "* <TIMESTAMP>\n\n" ) +
                       std::string( "* LEVEL=<LEVEL>\n" ) +
                       std::string( "* MESSAGE=<MESSAGE>\n" ) +
                       std::string( "* FILE=<FILE>\n" ) +
                       std::string( "* LINE=<LINE>\n" ) +
                       std::string( "***********************************\n" );
  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::GenericOutputStream stream(&std::cout, format );
  slic::addStreamToAllMsgLevels( &stream );

#endif
  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  // Mark begin of "initialization" phase
  GEOS_MARK_BEGIN("Initialization");

  ProblemManager problemManager( "ProblemManager", nullptr );
  problemManager.SetDocumentationNodes();
  problemManager.RegisterDocumentationNodes();

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );
  problemManager.ParseInputFile();



  problemManager.Initialize( &problemManager );

  GEOS_MARK_END("Initialization");

  gettimeofday(&tim, NULL);
  t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

  problemManager.ApplyInitialConditions();
  std::cout << std::endl << "Running simulation:" << std::endl;

  GEOS_MARK_BEGIN("RunSimulation");

  problemManager.RunSimulation();

  GEOS_MARK_END("RunSimulation");

  problemManager.ClosePythonInterpreter();

#ifdef USE_ATK
  axom::slic::finalize();
#endif

  gettimeofday(&tim, NULL);
  t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

  printf("Done!\n\nScaling Data: initTime = %1.2fs, runTime = %1.2fs\n", t_initialize - t_start,  t_run - t_initialize );

#ifdef USE_MPI
  MPI_Finalize();
#endif


  return 0;
}
