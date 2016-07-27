#include <slic/slic.hpp>
#include <slic/GenericOutputStream.hpp>
#include <iostream>

#include "dataRepository/WrapperCollection.hpp"
#include "PhysicsSolvers/SolidMechanicsLagrangianFEM.hpp"
#include "codingUtilities/SetSignalHandling.hpp"
#include "codingUtilities/stackTrace.hpp"
#include "managers/ProblemManager.hpp"

using namespace geosx;
using namespace asctoolkit;


int main( int argc, char *argv[] )
{
  std::cout<<"starting main"<<std::endl;
  geosx::setSignalHandling(geosx::stacktrace::handler1);


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

  ProblemManager problemManager;

  problemManager.ParseCommandLineInput( argc, argv );


#ifdef NDEBUG
  std::cout<<"NDEBUG is defined"<<std::endl;
#else
  std::cout<<"NDEBUG is NOT defined"<<std::endl;
#endif


  problemManager.ParseInputFile();





  slic::finalize();


  std::cout<<"exiting main"<<std::endl;

  return 0;
}
