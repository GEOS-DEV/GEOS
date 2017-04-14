#include <slic/slic.hpp>
#include <slic/GenericOutputStream.hpp>
#include <iostream>
//#include "ManagedArray.hpp"
#include "../dataRepository/ManagedGroup.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"
#include <mpi.h>

using namespace geosx;
using namespace asctoolkit;


int main( int argc, char *argv[] )
{

#if USE_MPI

  MPI_Init(&argc,&argv);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
#if SRC_INTERNAL
  MPI_DOMAIN_COMM=MPI_COMM_WORLD;
#endif
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




  ProblemManager problemManager( "ProblemManager", nullptr );

  problemManager.BuildDataStructure(nullptr);
  problemManager.SetDocumentationNodes( &problemManager );

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );
  problemManager.ParseInputFile();


  
  problemManager.InitializeObjects();
  problemManager.RunSimulation();
  

  problemManager.ClosePythonInterpreter();

  slic::finalize();


  std::cout<<"exiting main"<<std::endl;

  return 0;
}
