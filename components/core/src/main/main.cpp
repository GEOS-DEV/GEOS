#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"
#include <iostream>

#include "../dataRepository/intrinsic/WrapperCollection.hpp"
#include "PhysicsSolvers/NewtonianMechanics.hpp"
#include "codingUtilities/SetSignalHandling.hpp"
#include "codingUtilities/stackTrace.hpp"

using namespace asctoolkit;


int main()
{
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



  geosx::dataRepository::WrapperCollection domain("domain") ;

  std::string newName("new solver");
//  auto solver = geosx::SolverBase::SolverFactory::Factory(geosx::NewtonianMechanics::CatalogueName(),newName);
  auto solver = geosx::SolverBase::CatalogueEntryBase::Factory(geosx::NewtonianMechanics::CatalogueName(),newName);

  std::cout<<"breakpoint 1: "<<LOCATION<<std::endl;

  solver->RegisterDataObjects( domain );



//  geosx::DataObjectManager& nodes = domain.GetDataObjectManager<geosx::DataObjectManager>("FEM_Nodes");


  solver->TimeStep( 0.0, 0.1, 0, domain );




  slic::finalize();



  return 0;
}
