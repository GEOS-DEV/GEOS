/*
 * ProblemManager.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#include "ProblemManager.hpp"
#include "PhysicsSolvers/SolverBase.hpp"

namespace geosx
{

ProblemManager::ProblemManager():
    WrapperCollection("ProblemManager",nullptr)
{
  RegisterChildWrapperCollection<WrapperCollection>("domain");
  RegisterChildWrapperCollection<WrapperCollection>("solvers");
}

ProblemManager::~ProblemManager()
{}

void ProblemManager::ParseCommandLineInput( int const& argc, char* const argv[])
{}

void ProblemManager::ParseInputFile()
{
  geosx::dataRepository::WrapperCollection& solvers = GetChildDataObjectManager<geosx::dataRepository::WrapperCollection>("solvers");
  geosx::dataRepository::WrapperCollection& domain  = GetChildDataObjectManager<geosx::dataRepository::WrapperCollection>("domain");

  std::string newName("new solver");
  std::string newName2("new solver2");

  std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory("SolidMechanics_LagrangianFEM", newName, &domain );
  auto& solver1 = solvers.RegisterChildWrapperCollection<SolverBase>( newName, std::move(solver) );
  solver = SolverBase::CatalogInterface::Factory( "NewComponent", newName2, &domain );
  auto& solver2 = solvers.RegisterChildWrapperCollection<SolverBase>( newName2, std::move(solver) );

  solver1.Registration( domain );


  double time = 0.0;
  double dt = 5.0e-5;
  for( int i=0 ; i<10 ; ++i )
  {
    solver1.TimeStep( time, dt, i, domain );
    time += dt;
    std::cout<<std::endl;
  }

}

void ProblemManager::ApplySchedulerEvent()
{

}


} /* namespace geosx */
